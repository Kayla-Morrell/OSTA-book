# Loupe Browser (Visium)


## Overview

Loupe Browser is a desktop application from 10x Genomics that allows you to visualize your gene expression data without having to write code. You can utilize this software for most types of single cell transcriptomic data generated from 10x protocols, however we will discuss specifically how to use it for data generated from the Visium Spatial Gene Expression protocol. In general, you can use the Loupe Browser to align gene expression spots to histology images, look for marker gene expression, annotate populations, and cluster with three different clustering methods. 

[Here](https://support.10xgenomics.com/spatial-gene-expression/software/visualization/latest/analysis) is a tutorial from 10x Genomics on how to use the Loupe Browser.


## Manual alignment of images

One of the crucial first steps for processing Visium data is to align the gene expression spots to a high resolution (we use 40x) image of the tissue. While Space Ranger can do this automatically, if you want to ensure a high quality alignment with no mistakes it is best to do this manually in the Loupe browser. First, you upload the image and enter the serial number for the slide:

<div class="figure" style="text-align: center">
<img src="images/alignment_load_image.png" alt="Figure caption." width="75%" />
<p class="caption">(\#fig:unnamed-chunk-1)Figure caption.</p>
</div>

Then, you align the fiducial frame:

<div class="figure" style="text-align: center">
<img src="images/alignment_fiducial.png" alt="Figure caption." width="75%" />
<p class="caption">(\#fig:unnamed-chunk-2)Figure caption.</p>
</div>

Then, you manually select the spots that contain tissue:

<div class="figure" style="text-align: center">
<img src="images/alignment_tissue.png" alt="Figure caption." width="75%" />
<p class="caption">(\#fig:unnamed-chunk-3)Figure caption.</p>
</div>


## Output files and Space Ranger

The output files from Space Ranger are detailed in the Space Ranger section of this book. The `.cloupe` file is the one you need to import into the Loupe Browser in order to explore your processed data. For us they are generally between 1 and 2 GB. 


## Downstream analyses in Loupe Browser

You can look at the dimensionality reduction of the data through either a t-SNE or UMAP as well as apply graph-based or k-means clustering. 

<div class="figure" style="text-align: center">
<img src="images/tsne_kmeans.png" alt="Figure caption." width="75%" />
<p class="caption">(\#fig:unnamed-chunk-4)Figure caption.</p>
</div>

Most importantly you can overlay the gene expression data or annotated clusters onto the histology image. 

<div class="figure" style="text-align: center">
<img src="images/spatial_graph_based.png" alt="Figure caption." width="75%" />
<p class="caption">(\#fig:unnamed-chunk-5)Figure caption.</p>
</div>

You can make genes lists and plot marker genes spatially.

<div class="figure" style="text-align: center">
<img src="images/gene_expression.png" alt="Figure caption." width="75%" />
<p class="caption">(\#fig:unnamed-chunk-6)Figure caption.</p>
</div>

Or make violin plots.

<div class="figure" style="text-align: center">
<img src="images/violin_plots.png" alt="Figure caption." width="75%" />
<p class="caption">(\#fig:unnamed-chunk-7)Figure caption.</p>
</div>

Lastly, you can look at differential expression.

<div class="figure" style="text-align: center">
<img src="images/deg.png" alt="Figure caption." width="75%" />
<p class="caption">(\#fig:unnamed-chunk-8)Figure caption.</p>
</div>


<!-- ##add part about annotating cell populations? -->

