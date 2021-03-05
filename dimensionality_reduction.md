# Dimensionality reduction {#dimensionality_reduction}


## Background

Chapter on dimensionality reduction



## Previous steps

*Code (hidden) to run steps from the previous chapters, to generate the `SpatialExperiment` object required for this chapter.*





## Principal component analysis (PCA)

Apply principal component analysis (PCA) to the set of top highly variable genes (HVGs) to reduce the dimensionality of the dataset, and retain the top 50 principal components (PCs) for further downstream analyses.

This is done for two reasons: (i) to reduce noise due to random variation in expression of biologically uninteresting genes, which are assumed to have expression patterns that are independent of each other, and (ii) to improve computational efficiency during downstream analyses.

We use the computationally efficient implementation of PCA provided in the `scater` package [@McCarthy2017-zd]. This implementation uses randomization, and therefore requires setting a random seed for reproducibility.


```r
# compute PCA
set.seed(123)
spe <- runPCA(spe, subset_row = top_hvgs)

reducedDimNames(spe)
```

```
## [1] "PCA"
```

```r
dim(reducedDim(spe, "PCA"))
```

```
## [1] 3582   50
```



## Uniform Manifold Approximation and Projection (UMAP)

We also run UMAP [@McInnes2018-lx] on the set of top 50 PCs and retain the top 2 UMAP components, which will be used for visualization purposes.


```r
# compute UMAP on top 50 PCs
set.seed(123)
spe <- runUMAP(spe, dimred = "PCA")

reducedDimNames(spe)
```

```
## [1] "PCA"  "UMAP"
```

```r
dim(reducedDim(spe, "UMAP"))
```

```
## [1] 3582    2
```

```r
# update column names for plotting functions
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)
```



## Visualizations

In the next section on clustering, we will use the reduced dimensions (PCA and UMAP) to plot cluster labels in reduced dimension space.


