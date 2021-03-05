# (PART) Analysis steps {-}

# Analysis steps

This part consists of several chapters for steps in a computational analysis pipeline for spatially resolved transcriptomics data. This includes quality control (QC), normalization, feature selection, dimensionality reduction, clustering, and identifying marker genes.

These steps require that the raw data has been loaded into R. In the previous part, we provide detailed instructions and examples showing how to do this for the 10x Genomics Visium platform.

Throughout these chapters, we follow the Bioconductor principle of modularity -- if you prefer a different method for one of the steps, this can be substituted, the output stored in the `SpatialExperiment` object, and you can continue with the rest of the analysis pipeline. In several sections we describe several alternative methods within and outside Bioconductor.



## Load data

In the following analysis chapters, we use a pre-prepared dataset where we have previously applied the steps described in [Preprocessing steps], and saved the object in the `SpatialExperiment` format. This is available from the [STexampleData](https://github.com/lmweber/STexampleData) package.

The dataset consists of a single sample of human brain from the dorsolateral prefrontal cortex (DLPFC) region, measured using the 10x Genomics Visium platform, from our publication @Maynard2021. The dataset is described in more detail in [Human DLPFC workflow].

Here, we show how to load the data from [STexampleData](https://github.com/lmweber/STexampleData).


```r
library(SpatialExperiment)
library(STexampleData)

# load object
spe <- load_data("Visium_humanDLPFC")
```



## SpatialExperiment object structure

Inspect the `SpatialExperiment` object to understand the structure. For more details on `SpatialExperiment`, see [SpatialExperiment].


```r
# inspect object
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
## colData names(10): barcode_id imagerow ... pxl_col_in_fullres
##   pxl_row_in_fullres
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
## spatialData names(4) : barcode_id in_tissue x y
## imgData names(4): sample_id image_id data scaleFactor
```

```r
dim(spe)
```

```
## [1] 33538  4992
```

```r
assayNames(spe)
```

```
## [1] "counts"
```

```r
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
head(colData(spe))
```

```
## DataFrame with 6 rows and 10 columns
##                            barcode_id  imagerow  imagecol cell_count
##                           <character> <numeric> <numeric>  <integer>
## AAACAACGAATAGTTC-1 AAACAACGAATAGTTC-1        NA        NA         NA
## AAACAAGTATCTCCCA-1 AAACAAGTATCTCCCA-1   381.098   440.639          6
## AAACAATCTACTAGCA-1 AAACAATCTACTAGCA-1   126.328   259.631         16
## AAACACCAATAACTGC-1 AAACACCAATAACTGC-1   427.768   183.078          5
## AAACAGAGCGACTCCT-1 AAACAGAGCGACTCCT-1   186.814   417.237          2
## AAACAGCTTTCAGAAG-1 AAACAGCTTTCAGAAG-1   341.269   152.700          4
##                    ground_truth     sample_id array_row array_col
##                        <factor>   <character> <integer> <integer>
## AAACAACGAATAGTTC-1       NA     sample_151673         0        16
## AAACAAGTATCTCCCA-1       Layer3 sample_151673        50       102
## AAACAATCTACTAGCA-1       Layer1 sample_151673         3        43
## AAACACCAATAACTGC-1       WM     sample_151673        59        19
## AAACAGAGCGACTCCT-1       Layer3 sample_151673        14        94
## AAACAGCTTTCAGAAG-1       Layer5 sample_151673        43         9
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
head(spatialData(spe))
```

```
##                            barcode_id in_tissue    x     y
## AAACAACGAATAGTTC-1 AAACAACGAATAGTTC-1         0 3913 10880
## AAACAAGTATCTCCCA-1 AAACAAGTATCTCCCA-1         1 9791  4847
## AAACAATCTACTAGCA-1 AAACAATCTACTAGCA-1         1 5769 10508
## AAACACCAATAACTGC-1 AAACACCAATAACTGC-1         1 4068  3810
## AAACAGAGCGACTCCT-1 AAACAGAGCGACTCCT-1         1 9271  9164
## AAACAGCTTTCAGAAG-1 AAACAGCTTTCAGAAG-1         1 3393  5732
```

```r
imgData(spe)
```

```
## DataFrame with 2 rows and 4 columns
##       sample_id    image_id   data scaleFactor
##     <character> <character> <list>   <numeric>
## 1 sample_151673      lowres   ####   0.0450045
## 2 sample_151673       hires   ####   0.1500150
```


