# SPOT-Mut
High level information on the bioinformatics pipeline that was used for the analysis as reported in [Predicting clonal somatic mutations from tumour spatial transcriptomics data using a graph convolutional neural network.](https://doi.org/10.64898/2026.07.08.737173)

# Environment
  
Analysis was performed in a Conda/[Mamba](https://mamba.readthedocs.io/en/latest/index.html) environment using both R and Python scripts. See **SPOT-Mut.yml** for details. **scripts/Rpacks** describes R packages that were installed independently from Conda.
```{bash, eval=F}
mamba env create -f SPOT-Mut.yml -p SPOT-Mut          
```

# Data 

## Downloaded data

The following data were downloaded from external sources:
- TCGA mutation data was downloaded as the [MC3 maf file](https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc) from the GDC API.
- TCGA TPM transformed expression counts were obtained from the [GDC data portal](https://portal.gdc.cancer.gov/analysis_page?app=Downloads) as 'Gene Expression Quantification' STAR - counts .tsv files. The 'tmp' column was used to extract TPM expression counts. Only protein coding genes were kept.
- CCLE TPM count data was obtained from the [DepMap portal](https://depmap.org/portal/data_page/?tab=allData), release 22Q2. [download link](https://ndownloader.figshare.com/files/34989919)

## Downstream data

Processed data are available in data/ folder
Spatial transcriptomics data are available from ArrayExpress (Accession number TBA)
  
## Data processing

### Simulation procedure of an in-silico Spatial Transcriptomics expriment

This code was used to create the simulated ST slides from TCGA bulk RNA-seq samples on which SPOT-Mut was trained.
```
scripts/functions/CreateSTGraph.py
```

### TCGA sample selection

Sample selection procedure to select TCGA samples. Next, it shows the train/validation split and generates simulated ST slides.
The input count dataset (`dataset` variable) is the log1p-transformed TCGA expression counts, obtained as described above. Following genes were removed:
- non-coding genes
- Mitochondrial genes
- genes not also profiled in DepMap CCLE (see above for download) or a [representative 10X Genomics Visium experiment](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4565825).

```
scripts/TCGA_selection_and_preprocessing.py
```

### SPOT-Mut model architecture

The SPOT-Mut model can be found in the GCN python class in scripts/functions/GCN_class.py.
Weights for the TP53 model are found in data/TP53_model.pt
Weights for other models for other genes can be downloaded from [the Shiny application](https://ccgg.ugent.be/shiny/spot-mut/).

Load model weights with following code:
```python
import torch
# For the TP53 model
model = torch.load('data/TP53_model.pt', weights_only=False)
# For other driver genes
gene = 'KRAS' # Change the gene to the model of choice
model = GCN(n_features = len(gene_set))
saved_model = torch.load(gene+'_combat_weight_nocnv.pt', weights_only=False)
model.load_state_dict(saved_model["model_state_dict"])
```

### WES Aligment

Code to align the whole exome sequencing data. Both for the bulk sequencing experiment (`Align_WES.sh`), as well as the targeted microdissections (`Align_LCM_WES.sh`).

```{bash, eval=F}
scripts/Align_WES.sh
scripts/Align_LCM_WES.sh
```

### WES Mutation calling and functional annotation

Mutation calling pipeline with [GATK Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2). Functional annotation with [GATK Functotator](https://gatk.broadinstitute.org/hc/en-us/articles/360037224432-Funcotator) (data source version v1.7.20200521s). The tumor-only functionality of Mutect2 was used when running the following script.

```{bash, eval=F}
scripts/Mutect_TN.sh
```

### GLM training and validation

```
scripts/GLM_training.R
```

### Graph Neural Network training and validation

```
scripts/GNN_All_genes_training.py
```
! Beware !  Raw data on which the model trained is too large to be uploaded to GitHub. The data loader is empty. The user can define the `expr` and `expr_test` variables to use your own training and validation data.

# Manuscript Figures

All scripts to generate the figures as reported in the manuscript are available in the `scripts` folder under the appropriate filename. Required data files can be found in the `data` folder and appropriate subfolders. Paths to data files within the scripts should refer to the correct file in the data folder.

