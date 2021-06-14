# (PART) Analysis steps {-}

# Analysis steps

This part consists of several chapters for steps in a computational analysis pipeline for spatially resolved transcriptomics (ST) data. This includes quality control (QC), normalization, feature selection, dimensionality reduction, clustering, and identifying marker genes.

These steps require that the raw data has been loaded into R. In the previous part, we provide instructions and examples showing how to do this for the 10x Genomics Visium platform.

Throughout these chapters, we follow the Bioconductor principle of modularity -- for steps where multiple alternative analysis tools exist, we show how to substitute these, using the `SpatialExperiment` class to store input and output objects.


## Load data

In the following analysis chapters, we use a pre-prepared dataset where we have previously applied the steps described in [Preprocessing steps], and saved the object in the `SpatialExperiment` format. This is available from the [STexampleData](https://github.com/lmweber/STexampleData) package.

The dataset consists of a single sample of human brain from the dorsolateral prefrontal cortex (DLPFC) region, measured using the 10x Genomics Visium platform, from our publication @Maynard2021. The dataset is described in more detail in [Visium human DLPFC workflow].

Here, we show how to load the data from the [STexampleData](https://github.com/lmweber/STexampleData) package.


```r
library(SpatialExperiment)
library(STexampleData)

# load object
spe <- Visium_humanDLPFC()
```


## SpatialExperiment object structure

Next, we inspect the `SpatialExperiment` object. For more details, see [SpatialExperiment].


```r
# check object
spe
```

```
## class: SpatialExperiment 
## dim: 33538 4992 
## metadata(0):
## assays(1): counts
## rownames(33538): ENSG00000243485 ENSG00000237613 ... ENSG00000277475
##   ENSG00000268674
## rowData names(3): gene_id gene_name feature_type
## colnames(4992): AAACAACGAATAGTTC-1 AAACAAGTATCTCCCA-1 ...
##   TTGTTTGTATTACACG-1 TTGTTTGTGTAAATTC-1
## colData names(3): cell_count ground_truth sample_id
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
## spatialData names(6) : barcode_id in_tissue ... pxl_col_in_fullres
##   pxl_row_in_fullres
## spatialCoords names(2) : x y
## imgData names(4): sample_id image_id data scaleFactor
```

```r
# number of features (rows) and spots (columns)
dim(spe)
```

```
## [1] 33538  4992
```

```r
# names of 'assay' tables
assayNames(spe)
```

```
## [1] "counts"
```

```r
# features metadata
head(rowData(spe))
```

```
## DataFrame with 6 rows and 3 columns
##                         gene_id   gene_name    feature_type
##                     <character> <character>     <character>
## ENSG00000243485 ENSG00000243485 MIR1302-2HG Gene Expression
## ENSG00000237613 ENSG00000237613     FAM138A Gene Expression
## ENSG00000186092 ENSG00000186092       OR4F5 Gene Expression
## ENSG00000238009 ENSG00000238009  AL627309.1 Gene Expression
## ENSG00000239945 ENSG00000239945  AL627309.3 Gene Expression
## ENSG00000239906 ENSG00000239906  AL627309.2 Gene Expression
```

```r
# spots metadata
head(colData(spe))
```

```
## DataFrame with 6 rows and 9 columns
##                    cell_count ground_truth     sample_id         barcode_id
##                     <integer>     <factor>   <character>        <character>
## AAACAACGAATAGTTC-1         NA       NA     sample_151673 AAACAACGAATAGTTC-1
## AAACAAGTATCTCCCA-1          6       Layer3 sample_151673 AAACAAGTATCTCCCA-1
## AAACAATCTACTAGCA-1         16       Layer1 sample_151673 AAACAATCTACTAGCA-1
## AAACACCAATAACTGC-1          5       WM     sample_151673 AAACACCAATAACTGC-1
## AAACAGAGCGACTCCT-1          2       Layer3 sample_151673 AAACAGAGCGACTCCT-1
## AAACAGCTTTCAGAAG-1          4       Layer5 sample_151673 AAACAGCTTTCAGAAG-1
##                    in_tissue array_row array_col pxl_col_in_fullres
##                    <integer> <integer> <integer>          <integer>
## AAACAACGAATAGTTC-1         0         0        16               2435
## AAACAAGTATCTCCCA-1         1        50       102               8468
## AAACAATCTACTAGCA-1         1         3        43               2807
## AAACACCAATAACTGC-1         1        59        19               9505
## AAACAGAGCGACTCCT-1         1        14        94               4151
## AAACAGCTTTCAGAAG-1         1        43         9               7583
##                    pxl_row_in_fullres
##                             <integer>
## AAACAACGAATAGTTC-1               3913
## AAACAAGTATCTCCCA-1               9791
## AAACAATCTACTAGCA-1               5769
## AAACACCAATAACTGC-1               4068
## AAACAGAGCGACTCCT-1               9271
## AAACAGCTTTCAGAAG-1               3393
```

```r
# spatial metadata
head(spatialData(spe))
```

```
## DataFrame with 6 rows and 6 columns
##                            barcode_id in_tissue array_row array_col
##                           <character> <integer> <integer> <integer>
## AAACAACGAATAGTTC-1 AAACAACGAATAGTTC-1         0         0        16
## AAACAAGTATCTCCCA-1 AAACAAGTATCTCCCA-1         1        50       102
## AAACAATCTACTAGCA-1 AAACAATCTACTAGCA-1         1         3        43
## AAACACCAATAACTGC-1 AAACACCAATAACTGC-1         1        59        19
## AAACAGAGCGACTCCT-1 AAACAGAGCGACTCCT-1         1        14        94
## AAACAGCTTTCAGAAG-1 AAACAGCTTTCAGAAG-1         1        43         9
##                    pxl_col_in_fullres pxl_row_in_fullres
##                             <integer>          <integer>
## AAACAACGAATAGTTC-1               2435               3913
## AAACAAGTATCTCCCA-1               8468               9791
## AAACAATCTACTAGCA-1               2807               5769
## AAACACCAATAACTGC-1               9505               4068
## AAACAGAGCGACTCCT-1               4151               9271
## AAACAGCTTTCAGAAG-1               7583               3393
```

```r
# spatial coordinates
head(spatialCoords(spe))
```

```
##                       x    y
## AAACAACGAATAGTTC-1 3913 2435
## AAACAAGTATCTCCCA-1 9791 8468
## AAACAATCTACTAGCA-1 5769 2807
## AAACACCAATAACTGC-1 4068 9505
## AAACAGAGCGACTCCT-1 9271 4151
## AAACAGCTTTCAGAAG-1 3393 7583
```

```r
# image metadata
imgData(spe)
```

```
## DataFrame with 2 rows and 4 columns
##       sample_id    image_id   data scaleFactor
##     <character> <character> <list>   <numeric>
## 1 sample_151673      lowres   ####   0.0450045
## 2 sample_151673       hires   ####   0.1500150
```

