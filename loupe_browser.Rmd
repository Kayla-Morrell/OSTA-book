# Loupe Browser (Visium)

Draft structure and points to cover (but feel free to change)


## Overview
Loupe Browser is a desktop application from 10x Genomics that allows you to visiual your gene expression data without having to write code. You can utilize this software for most types of single cell transcriptomic data generated from 10x protocols, however we will discuss specifically how to use it for data generated from the Visisun Spatial Gene Expression protocol. In general, you can use the Loupe Browser to align gene expression spots to histology images, look for marker gene expression, annotate populations, and cluster with 3 different clustering methods. 

Below is tutorial from 10x on how to use the Loupe Browser:
https://support.10xgenomics.com/spatial-gene-expression/software/visualization/latest/analysis


## Manual alignment of images
One of the crucial first steps for processing Visium data is to align the gene expression spots to a high resolution (we use 40x) image of the tissue. While Space Ranger can do this automatically, if you want to ensure a high quality alignment with no mistakes it's best to do this manually in the loupe browser. First, you upload the image and enter the serial number for the slide:

(knitr::include_graphics("images/alignment_load_image.png"))

Then, you align the fiducial frame:

(knitr::include_graphics("images/alignment_fiducial.png"))

Then, you manually select the spots that contain tissue:

(knitr::include_graphics("images/alignment_tissue.png"))

## Output files for Space Ranger
The output files from space ranger are detailed in the space ranger section of this book. The .cloupe file is the one you need to import into the Loupe browser in order to explore your processed data. For us they are generally between 1 and 2 GB. 

## Downstream analyses in Loupe Browser
You can look at the dimensionality reduction of the data through either a t-SNE or UMAP as well as apply graph-based or k-means clustering. 

(knitr::include_graphics("images/tsne_kmeans.png"))

Most importantly you can overlay the gene expression data or annotated clusters onto the histology image. 

(knitr::include_graphics("images/spatial_graph_based.png"))

You can make genes lists and plot marker genes spatially.

(knitr::include_graphics("images/gene_expression.png"))

Or make violin plots.

(knitr::include_graphics("images/violin_plots.png"))

Lastly, you can look at differential expression.

(knitr::include_graphics("images/deg.png"))

##add part about annotating cell populations?



