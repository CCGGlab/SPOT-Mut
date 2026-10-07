#### A ####
library(fmsb)
source('functions/custom_radar.R')


data <- as.data.frame(matrix(c(0.9651908437576288, 0.9740222122087887, 0.9772197104755769, 0.8937474154697452, 0.9456312616259845, 0.9770272930163703, 0.9397493258435817, 0.9894109925614568, 0.9552065054425999, 0.9746483105615595, 0.9353009450252984, 0.9779534889511902, 0.962683531576758,
                               0.9708425027851746, 0.9460701757441708, 0.9598458317123182, 0.9901711931872788, 0.9793313242193801, 0.9621995241546314, 0.988467165008557, 0.8555085408184091, 0.9482046256205198, 0.9584419716296796, 0.955056590326736, 0.9787163232103004, 0.9827061965079312),
                               #0.9679745588491163, 0.9598041784955216, 0.9684209974417085, 0.9393876867838192, 0.9621494267462671, 0.9695212839838021, 0.9634595233483232, 0.9174613665871404, 0.9516421316293137, 0.9664425432154737, 0.945009693949516, 0.9783143020045858, 0.9725748484954075), 
                             nrow=2, byrow = TRUE))
colnames(data) <- c("Pan-Cancer", "BLCA" , "HNSC" , "STAD" , "COAD" , "BRCA","LIHC","SKCM","PRAD","LUAD","LUSC","LGG","UCEC")
rownames(data) <- c('Precision', 'Recall')

# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each variable to show on the plot!
data <- rbind(rep(1,13) , rep(0,13) , data)

# Color vector
colors_border=c( "#ED7D31", "#A6183A" )

# plot with default options:
radarchart2( data  , axistype=1, 
             #custom polygon
             pcol=colors_border , pfcol=NA , plwd=4 , plty=1,
             #custom the grid
             cglcol="darkgrey", cglty=1, axislabcol="black", caxislabels=c(0,0.25,0.5,0.75,1), cglwd=0.8,
             #custom labels
             vlcex=0.8, calcex = 0.65, vlabcol = rep('black',13)
)

# Add a legend
legend(x=1, y=1.3, legend = rownames(data[-c(1,2),]), bty = "n", pch=16 , col=colors_border , cex=1, pt.cex=1)
#####

## 
data <- as.data.frame(matrix(c(0.9347,0.9258,0.9136,0.9685,0.8589,0.9271,0.9372,0.9446,0.8958,0.9551,0.4023,0.9519,
                               0.9733,0.9684,0.9757,0.9709,0.9903,0.9879,0.9782,0.9515,0.9806,0.9806,1,0.9612),
                               nrow=2, byrow = TRUE))
colnames(data) <- c("APC","ARID1A","ATRX","BRAF","CASP8","CIC","CTNNB1","KRAS","NFE2L2","PTEN","RNF43","STED2")
rownames(data) <- c('Precision', 'Recall')

# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each variable to show on the plot!
data <- rbind(rep(1,13) , rep(0,13) , data)

# Color vector
colors_border=c( "#ED7D31", "#A6183A" )

# plot with default options:
radarchart2( data  , axistype=1, 
             #custom polygon
             pcol=colors_border , pfcol=NA , plwd=4 , plty=1,
             #custom the grid
             cglcol="darkgrey", cglty=1, axislabcol="black", caxislabels=c(0,0.25,0.5,0.75,1), cglwd=0.8,
             #custom labels
             vlcex=0.8, calcex = 0.65, vlabcol = rep('black',13)
)

# Add a legend
legend(x=1, y=1.3, legend = rownames(data[-c(1,2),]), bty = "n", pch=16 , col=colors_border , cex=1, pt.cex=1)

## exponential confidence decrease

confidence <- read.table('data/confidance_results.tsv', sep = '\t', header = T)
confidence$x <- confidence$key/max(confidence$key)

plot(confidence$key, confidence$auc, col = 'blue', ylim = c(0,1))
points(confidence$key, confidence$f1, col = 'red')
points(confidence$key, confidence$accuracy, col = 'orange')

fit_logistic <- nls(
  accuracy ~ lower + (upper - lower) /
    (1 + exp(k * (x - midpoint))),
  data = confidence,
  start = list(
    lower = 0,
    upper = 1,
    midpoint = confidence$accuracy[1],
    k = 80
  )
)

x_grid <- seq(min(confidence$x), max(confidence$x), length.out = 1000)

plot(confidence$x, confidence$accuracy,
     xlab = "genes removed", ylab = "F1 score")

lines(
  x_grid,
  predict(fit_logistic, newdata = data.frame(x = x_grid)),
  col = "skyblue",
  lwd = 3
)
# saveRDS(fit_logistic, "shiny/data/accuracy_logistic_model.rds")
