library(matrixStats)
library(dplyr)
library(caret)
library(glmnet)
library(pROC)
library(doMC)
library(future.apply)
library(R.utils)
library(sva)
## batch correct TCGA data according to cancer type
# exprs_batch <- t(ComBat(t(exprs), batch = samples[rownames(exprs),'DISEASE']))

fit_gene_classifier <- function(
    target_gene,
    gene_class_vec,
    expr,
    mutation,
    CNV_loss,
    CNV_gain,
    samples,
    mut_burden,
    n_features = 8000,
    usedummy = F,
    min_positive = 15,
    min_prevalence = 0.1,
    min_samples = 199,
    min_samples_per_ct = 100,
    train_frac = 0.8,
    useweights = F,
    useparallel = F,
    cores = 32,
    seed = 1
) {
  ## target_gene: string or vector of gene symbols to predict mutation status for
  ## gene_class_vec: named vector of gene classes (TSG, Oncogene, Other) for each gene symbol
  ## expr: expression dataframe, rownames = sample IDs, colnames = gene symbols
  ## CNV_gain: binary CNV gain dataframe, rownames = sample IDs, colnames = gene symbols
  ## CNV_loss: binary CNV loss dataframe, rownames = sample IDs, colnames = gene symbols
  ## mutation: binary mutation dataframe, rownames = sample IDs, colnames = gene symbols
  ## samples: sample ID's you want to include (rownames of expr and mutation)
  ## n_features: either a number (e.g. 8000) or a vector of gene symbols to use as features
  ## usedummy: whether to include dummy variables for cancer type in the model
  ## min_positive: minimum number of positive samples per cancer type to include in the model
  ## min_prevalence: minimum prevalence of positive samples per cancer type to include in the model
  ## min_samples: minimum number of samples to include in the model
  ## min_samples_per_ct: minimum number of samples per cancer type to include in the model
  ## train_frac: fraction of samples to use for training
  ## useweights: whether to use weights to account for sample class imbalance during model fitting
  ## useparallel: whether to use parallel processing for cross-validation
  ## cores: number of cores to use for parallel processing
  ## seed: random seed for reproducibility
  
  if (!all(target_gene %in% colnames(mutation))) {
    warning('not target genes are in the mutations dataframe')
    return(NULL)
  }
  mutation <- mutation[rownames(expr),]
  samples <- samples[rownames(expr),]
  CNV_gain <- CNV_gain[rownames(expr),]
  CNV_loss <- CNV_loss[rownames(expr),]
  message(paste(target_gene, collapse = ' '))
  
  ## Remove target gene from predictors
  expr_sub <- expr[, setdiff(colnames(expr), target_gene), drop = FALSE]
  
  ## Select most variable genes
  if(is.vector(n_features)){
    top_genes <- intersect(x = n_features, colnames(expr_sub))
  } else if (is.numeric(n_features)) {
    mad_scores <- colMads(as.matrix(expr_sub))
    top_genes <- names(sort(mad_scores, decreasing = TRUE))[1:n_features]
  }
  expr_sub <- expr_sub[, top_genes, drop = FALSE]
  ## Mutation status
  if (length(target_gene) == 1) {
    gene_status_mut <- mutation[, target_gene] == 1
  } else {
    gene_status_mut <- rowSums(mutation[, target_gene, drop = FALSE]) >= 1
  }

  tsg_genes <- target_gene[gene_class_vec[target_gene] == "TSG"]
  onc_genes <- target_gene[gene_class_vec[target_gene] == "Oncogene"]

  gene_status_tsg <- NULL
  if (length(tsg_genes) > 0 && !is.null(CNV_loss) && target_gene %in% colnames(CNV_loss)) {
    gene_status_tsg <- rowSums(
      CNV_loss[, tsg_genes, drop = FALSE] == 1
    ) >= 1
  }

  gene_status_onc <- NULL
  if (length(onc_genes) > 0 && !is.null(CNV_gain) && target_gene %in% colnames(CNV_gain)) {
    gene_status_onc <- rowSums(
      CNV_gain[, onc_genes, drop = FALSE] == 1
    ) >= 1
  }
  
  gene_status <- gene_status_mut
  
  if (!is.null(gene_status_tsg)) {
    gene_status <- gene_status | gene_status_tsg
  } else if (!is.null(CNV_loss)){
    warning("No CNV loss information included")
  }
  
  if (!is.null(gene_status_onc)) {
    gene_status <- gene_status | gene_status_onc
  } else if (!is.null(CNV_gain)){
    warning("No CNV gain information included")
  }
  gene_status <- as.integer(gene_status)

  tmp <- data.frame(
    SAMPLE = rownames(expr_sub),
    DISEASE = samples$DISEASE,
    STATUS = gene_status
  )
  
  keep_ct <-
    tmp %>%
    group_by(DISEASE) %>%
    summarise(
      n_samples = n(),
      positives = sum(STATUS),
      prevalence = mean(STATUS),
      .groups = "drop"
    ) %>%
    filter(
      n_samples > min_samples_per_ct,
      positives >= min_positive,
      prevalence >= min_prevalence
    ) 
  if (nrow(keep_ct)>1){
    keep_ct <- pull(keep_ct, DISEASE)
  } else {
    return(NULL)
  }
  
  keep <- samples$DISEASE %in% keep_ct
  if (sum(keep) < min_samples) {
    return(NULL)
  }
  
  expr_sub <- expr_sub[keep, , drop = FALSE]
  gene_status <- gene_status[keep]
  mut_burden_sub <- mut_burden[keep, , drop = FALSE]
  samples_sub <- samples[keep, , drop = FALSE]
  if (usedummy){
    X <- cbind(
      expr_sub,
      model.matrix(~ DISEASE - 1, data = samples_sub)
    )
  } else {
    X <- expr_sub
  }
  
  Y <- gene_status
  
  strata <- interaction(samples_sub$DISEASE, Y)
  set.seed(seed)
  train_idx <- createDataPartition(
    strata,
    p = train_frac,
    list = FALSE
  )
  
  X_train <- X[train_idx, , drop = FALSE]
  X_test <- X[-train_idx, , drop = FALSE]
  
  Y_train <- Y[train_idx]
  Y_test <- Y[-train_idx]
  
  train_df <- data.frame(
    X_train,
    STATUS = Y_train,
    DISEASE = samples_sub$DISEASE[train_idx]
  )
  
  test_df <- data.frame(
    X_test,
    STATUS = Y_test,
    DISEASE = samples_sub$DISEASE[-train_idx]
  )
  
  alpha_grid <- seq(0, 1, by = 0.1)
  weights <- NULL
  if (useparallel) registerDoMC(cores = cores)
  if (useweights){
    group_counts <- train_df %>%
      count(DISEASE, STATUS)
    
    train_df <- train_df %>%
      left_join(group_counts, by = c("DISEASE", "STATUS"))
    rownames(train_df) <- rownames(X)[train_idx]
    train_df$w <- 1 / train_df$n
    train_df$w <- train_df$w / mean(train_df$w)
    X_bal <- as.matrix(train_df[, !(names(train_df) %in% c("STATUS", "DISEASE", "n", "w"))])
    Y_bal <- train_df$STATUS
    weights <- train_df$w
    cvfits <- lapply(alpha_grid, function(a) {
      
      cv.glmnet(
        x = X_bal,
        y = Y_bal,
        family = "binomial",
        alpha = a,
        type.measure = "auc",
        nfolds = 5,
        weights = weights,
        parallel = useparallel
      )
      
    })
  } else {
    group_sizes <- train_df %>%
      count(STATUS)
    
    n_min <- min(group_sizes$n)
    
    balanced_train <-
      train_df %>%
      group_by(STATUS) %>%
      slice_sample(n = n_min) %>%
      ungroup()
    
    X_bal <- as.matrix(
      balanced_train[, !(names(balanced_train) %in% c("STATUS","DISEASE"))]
    )
    
    Y_bal <- balanced_train$STATUS
    
    ## Search over alpha
    cvfits <- lapply(alpha_grid, function(a) {
      
      cv.glmnet(
        x = X_bal,
        y = Y_bal,
        weights = weights,
        family = "binomial",
        alpha = a,
        type.measure = "auc",
        nfolds = 5,
        parallel = useparallel,
      )
      
    })
  }
  ## Best alpha
  cv_auc <- sapply(cvfits, function(f) max(f$cvm))
  
  best_idx <- which.max(cv_auc)
  
  best_alpha <- alpha_grid[best_idx]
  
  best_cv <- cvfits[[best_idx]]
  
  best_lambda <- best_cv$lambda.min
  
  final_model <- glmnet(
    x = as.matrix(X_train),
    y = Y_train,
    weights = weights,
    family = "binomial",
    alpha = best_alpha,
    lambda = best_lambda
  )
  
  prob_train <- predict(
    final_model,
    as.matrix(X_train),
    type = "response"
  )
  
  prob_test <- predict(
    final_model,
    as.matrix(X_test),
    type = "response"
  )
  
  auc_train <- roc(Y_train, as.numeric(prob_train))$auc
  auc_test <- roc(Y_test, as.numeric(prob_test))$auc

  
  cutoff <- as.numeric(
    coords(
      roc(Y_train, c(prob_train)),
      x = "best",
      best.method = "youden",
      ret = "threshold"
    )
  )
  
  test_class <- ifelse(prob_test >= cutoff, 1, 0)
  
  cm <- confusionMatrix(
    factor(test_class),
    factor(Y_test),
    positive = "1"
  )
  # Balanced Accuracy
  balanced_accuracy <- cm$byClass["Balanced Accuracy"]
  
  save_train <- paste0(paste(target_gene, collapse = '_'),"_train_combat_weight_nocnv.csv")
  write.table(train_df, file = save_train, sep = ',', quote = F)
  gzip(save_train, overwrite = T)
  
  save_test <- paste0(paste(target_gene, collapse = '_'),"_test_combat_weight_nocnv.csv")
  write.table(test_df, file = save_test, sep = ',', quote = F)
  gzip(save_test, overwrite = T)
  
  return(list(
    gene = target_gene,
    model = final_model,
    auc_train = auc_train,
    auc_test = auc_test,
    cutoff = cutoff,
    confusion = cm,
    Y_train = Y_train,
    Y_test = Y_test,
    Y_prob_test = as.numeric(prob_test),
    Y_prob_train = as.numeric(prob_train)
  ))
}

plan(multisession, workers = 64)
options(future.globals.maxSize = 10 * 1024^3)
cancer_gene_census <- read.csv("data/cancer_gene_census.csv")
GLM_Performance_C <- future_lapply(seq_len(nrow(cancer_gene_census)), function(i) {
  
  gene <- cancer_gene_census$Gene.Symbol[i]

  fit_gene_classifier(
    gene,
    NULL,
    combat_X, mutation, NULL, NULL,
    samples = samples,
    mut_burden = mut_burden,
    n_features = colnames(combat_X),
    useparallel = F,
    useweights = T,
    usedummy = T,
    seed = which(cancer_gene_census$Gene.Symbol == gene)
  )
  
}, future.seed = TRUE)
