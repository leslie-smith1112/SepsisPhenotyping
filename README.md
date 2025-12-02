Untitled
================

## GitHub Documents

This is an R Markdown format used for publishing markdown documents to
GitHub. When you click the **Knit** button all R code chunks are run and
a markdown file (.md) suitable for publishing to GitHub is generated.

## Including Code

You can include R code in the document as follows:

``` r
summary(cars)
```

    ##      speed           dist       
    ##  Min.   : 4.0   Min.   :  2.00  
    ##  1st Qu.:12.0   1st Qu.: 26.00  
    ##  Median :15.0   Median : 36.00  
    ##  Mean   :15.4   Mean   : 42.98  
    ##  3rd Qu.:19.0   3rd Qu.: 56.00  
    ##  Max.   :25.0   Max.   :120.00

## Including Plots

You can also embed plots, for example:

![](ReadMe_files/figure-gfm/pressure-1.png)<!-- -->

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot. \#### STRART
HERE

This repository is associated with the manuscript: Here we walk through
analysis presented in the paper, as well as provide examples for some
analysis. TODO

<figure>
<img src="images/Method.jpg" alt="Method Overview" />
<figcaption aria-hidden="true">Method Overview</figcaption>
</figure>

## Dataset Processing

We survey refine.bio, GEO, and SRA for transcriptomic datasets with
adult patients meeting sepsis/septic shock criteria. The datasets
included in final analysis are: [GSE100159](), [GSE13015](),
[GSE131761](), [GSE134347](), [GSE137340](), [GSE154918](),
[GSE185263](), [GSE189400](), [GSE196117](), [GSE199816](),
[GSE211210](), [GSE216902](), [GSE222393](), [GSE232404](),
[GSE232753](), [GSE236713](), [GSE236892](), [GSE32707](), [GSE33118](),
[GSE57065](), [GSE63311](), [GSE65682](), [GSE66890](), [GSE69063](),
[GSE74224](), [GSE95233](), [SRP198776](), [Guirgis et al.]().

We preprocess each dataset individually to ensure consistent format and
to maximize comparability. Metadata from all studies is consolidated
into a single file and then cleaned pre- and post- batch correction
using scripts `code/Metadata_Clean_PreBatchCorrection.R` and
`code/Metadata_Clean_PostBatchCorrection.R`.

<figure>
<img src="images/IndividualDatasetProcessing.png"
alt="Dataset Processing Pipeline" />
<figcaption aria-hidden="true">Dataset Processing Pipeline</figcaption>
</figure>

## Batch Correction

All samples from each datasets are used in in batch correction. Datasets
are combined into a single matrix using `code/Create_master_matrix.R`.
Batch correction was done in R using ComBat
(`code/Batch_correction_pca.R`). Processed expression data and metadata
for all 3,713 samples included in batch correction are included as Rdata
in this repo in data/ and provided at TODO add GEO ID labeled
`AllSampleExpressionSubmitted.tsv` and
`AllSampleMetadataSubmitted.xlsx`.

We can look at samples before and after batch correction here: ![Dataset
Batch Correction](images/batchCorrectionHor.png)

Samples values after batch correction are provided in GEO dataset: GSE
TODO.

We can read in all samples here:

``` r
all_metadata <- readRDS(here::here("data","all_samples_metadata.rds"))
all_expr <- readRDS(here::here("data","all_samples_expression.rds"))
```

Lets look at the metadata we have for this dataset. We have 68 columns
in this file, comprising all metadata available for each included study.

``` r
colnames(all_metadata)
```

    ##  [1] "Sample ID"                    "Dataset"                     
    ##  [3] "GEO Accession"                "Platform ID"                 
    ##  [5] "Source"                       "Reanalyzed by another study" 
    ##  [7] "Tissue"                       "Sample Class"                
    ##  [9] "Disease"                      "Disease 2"                   
    ## [11] "Disease Simplified"           "ARDs"                        
    ## [13] "Sepsis/Septic Shock"          "Timepoint"                   
    ## [15] "Is Repeat Sample"             "Low Platelets"               
    ## [17] "RNA-seq/Microarray"           "Treatment"                   
    ## [19] "Responder"                    "Survival"                    
    ## [21] "Age"                          "Gender"                      
    ## [23] "ifng il10 ratio"              "Ethinicity"                  
    ## [25] "Race"                         "Pathogen"                    
    ## [27] "Title"                        "Patient ID"                  
    ## [29] "Covid Status"                 "SOFA (24hr)"                 
    ## [31] "RepSOFA"                      "SAPSII"                      
    ## [33] "Organ Failure Count"          "Cardiovascular Failure"      
    ## [35] "Liver Failure"                "AKI"                         
    ## [37] "Coagulation Failure"          "Respiratory Failure"         
    ## [39] "SIRS Sign"                    "Culture Positive"            
    ## [41] "Pneumonia Diagnosis"          "Thrombocytonpenia"           
    ## [43] "Endotype Cohort"              "PreviouslyDefinedSubtypes"   
    ## [45] "Time to Event"                "ICU Aquired Infection"       
    ## [47] "ICU Aquired Infection Paired" "Diabeties Mellitus"          
    ## [49] "Abdominal Sepsis"             "WBC Cells per UL"            
    ## [51] "Direct Lung"                  "AECC ALI"                    
    ## [53] "Berlin ARDs"                  "Shock"                       
    ## [55] "Severity"                     "AgeSummary"                  
    ## [57] "UnderlyingDisease"            "RenalDysfunction"            
    ## [59] "LiverDysfunction"             "MechanicalVentilation"       
    ## [61] "Group"                        "InfectionOrigin"             
    ## [63] "OnPressors"                   "Immunocompromised"           
    ## [65] "Antibiotics"                  "DaySurvivalEdited"           
    ## [67] "TimetoEventEdited"            "MolecularSubtype"

Importantly, columns used in our analysis include:

| Column Name                 | Description                               |
|-----------------------------|-------------------------------------------|
| `Sample ID`                 | Sample ID                                 |
| `Dataset`                   | Associated Dataset ID                     |
| `Disease`                   | Unedited disease classification           |
| `Disease Simplified`        | Edited/simplified disease classification  |
| `Timepoint`                 | Time sample was taken                     |
| `Is Repeat Sample`          | Whether sample was first from the patient |
| `Sepsis/Septic Shock`       | Sepsis-shock classification               |
| `Shock`                     | Shock classification                      |
| `Survival`                  | Mortality classification                  |
| `Time to Event`             | Time to mortality event                   |
| `Gender`                    | Sample Sex                                |
| `Age`                       | Sample Age                                |
| `Race`                      | Sample Race                               |
| `PreviouslyDefinedSubtypes` | Subtype assigned by previous study        |

Lets look at all diseases present in this file:

``` r
table(all_metadata$`Disease`)
```

    ## 
    ##            Acute Sepsis             Anaphylaxis       Cardiogenic Shock 
    ##                       5                      33                      36 
    ##                 Control         Control/healthy        Control/Recovery 
    ##                      95                       3                       9 
    ## Control/type 2 diabetes            COVID Sepsis                COVID-19 
    ##                       7                      31                       9 
    ##             Head Trauma                 Healthy         Healthy control 
    ##                      30                     189                      33 
    ##               No Sepsis               No-Sepsis        Non-Septic Shock 
    ##                      35                      80                      33 
    ##                    None           Noninfectious           post-surgical 
    ##                       8                      59                      31 
    ##                    SAKI                  Sepsis        Sepsis Follow Up 
    ##                       5                    2030                       4 
    ##      Sepsis-Melioidosis  Sepsis-Other Infection         Sepsis-Recovery 
    ##                      18                      15                       2 
    ##      Sepsis/Melioidosis Sepsis/Other infections               SepsisBSI 
    ##                      24                      24                      56 
    ##            Sepsisperiph             Septic ARDs            Septic Shock 
    ##                      65                      39                     489 
    ##  Septic Shock Follow Up  Septic Shock with ARDs           Severe Sepsis 
    ##                      10                      20                      26 
    ##                    SIRS Uncomplicated Infection               untreated 
    ##                     114                      12                      34

Now lets look at the simplified disease column. We have taken the
diseases in `Disease` and cleaned them in this column for simplicity.

``` r
table(all_metadata$`Disease Simplified`)
```

    ## 
    ##                                              Control 
    ##                                                  344 
    ## Intubated subjects undergoing mechanical ventilation 
    ##                                                   34 
    ##                    Uncomplicated baterical infection 
    ##                                                   12 
    ##                                          Anaphylaxis 
    ##                                                   33 
    ##                          ICU noninfectious condition 
    ##                                                   59 
    ##                                          Head Trauma 
    ##                                                   30 
    ##                                     Critical Illness 
    ##                                                   80 
    ##                                             COVID-19 
    ##                                                    9 
    ##                                    Cardiogenic Shock 
    ##                                                   36 
    ##           Potential for sepsis but was not developed 
    ##                                                   35 
    ##                                                 SIRS 
    ##                                                  178 
    ##                                         Sepsis - AKI 
    ##                                                    5 
    ##                                               Sepsis 
    ##                                                 2263 
    ##                                        Sepsis - ARDs 
    ##                                                   99 
    ##                                       Sepsis - Shock 
    ##                                                  496

## Gene and Sample Selection

# Primary Pipeline Script Summary

    ┌──────────────────────────────────────────────────────────┐
    │                    Create_master_matrix.R                |  
    |                                                          |
    | PURPOSE:                                                 |
    |      ▷ Combine all datasets into a single matrix         |
    | INPUT:                                                   |
    |      ▷ 28 individual preprocessed datasets               |  
    │      ▷ Original manually curated metadata file           |                    
    │ OUTPUT:                                                  |
    |      ▷ data/raw_all_samples_expression.rds               |
    |         ┗ Single matrix with all datasets combined       |
    |      ▷ data/batchIDs_allSamples.rds                      |
    │         ┗ Batch IDs for use in ComBat                    |
    |      ▷ data/all_samples_metadata.rds                     |    
    |         ┗ Cleaned metadata file                          |
    └────────────────────────┬─────────────────────────────────┘
                             │
                             ▼
    ┌──────────────────────────────────────────────────────────┐
    │                  Batch_correction_pca.R                  |  
    |                                                          |
    | PURPOSE:                                                 |
    |      ▷ Batch correct datasets to remove technical        |
    |         variation                                        |
    | INPUT:                                                   |
    |      ▷ data/raw_all_samples_expression.rds               |
    │      ▷ data/batchIDs_allSamples.rds                      |  
    |      ▷ data/all_samples_metadata.rds
    |      ▷ data/disease_metadata.rds
    │ OUTPUT:                                                  |
    |      ▷ data/Init_PCA_values_allSamples.rds               |
    |         ┗ PCA values for all samples prior to batch      |  
    |           correction                                     |
    |      ▷ figures/Initial_PCA_Dataset.png                   |
    │         ┗ Scatterplot of samples prior to batch          |
    |           correction                                     |
    |         ┗ Colored by dataset                             |
    |      ▷ figures/Initial_PCA_Disease.png                   |
    │         ┗ Scatterplot of samples prior to batch          |
    |           correction                                     |
    |         ┗ Colored by disease                             |
    |      ▷ data/all_samples_expression.rds                   |
    │         ┗ Batch corrected expression matrix              |
    |           correction, all samples all genes              |               
    |      ▷ figures/BatchCorrected_PCA_values_allSamples.rds  |
    │         ┗ PCA values for all samples after batch         |
    |           correction                                     |
    |      ▷ figures/BatchCorrected_PCA_Dataset.png            |
    │         ┗ Scatterplot of samples colored by dataset      |  
    |           after batch correction                         |  
    |      ▷ figures/BatchCorrected_PCA_Disease.png            |
    │         ┗ Scatterplot of samples colored by disease      |
    |           after batch correction                         |    
    |      ▷ data/disease_expression_all_genes.rds"                      |
    │         ┗ Expression values for sepsis samples only      |
    |           after batch correction                         |
    |      ▷ figures/BatchCorrected_PCA_DiseaseOnly.png        |
    │         ┗ Scatterplot of sepsis samples colored by       |
    |           disease after batch correction                 |                
    |                                                          |
    └────────────────────────┬─────────────────────────────────┘
                             │
                             ▼
    ┌─────────────────────────────────────────────────────────────────────┐
    │                        CutMixedGenes.R                       |  
    |                                                                     |
    | PURPOSE:                                                            |
    |       ▷ Identifying genes to include for analysis       |
    | INPUT:                                                              |
    |      ▷ data/disease_expression_all_genes.rds
    │      ▷ data/disease_metadata.rds
    │ OUTPUT:                                                             |
    |      ▷ data/disease_expression.rds
                         |
    |         (all diseases/control) and all genes                        |
    |         (data/raw_all_samples_expression.rds)                       |
    |      ▷ Batch IDs for use in ComBat                                  |
    │         (data/batchIDs_allSamples.rds)                              |
    └────────────────────────┬────────────────────────────────────────────┘
                             │
                             ▼
                             
