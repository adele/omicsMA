## omicsMA R package: Variance-Preserving Estimation and Normalization of M-A Values from Omics Experiments


### Overview

This package provides methods for estimating and normalizing the M (intensity log-ratio) and A (mean log intensity) values from two-channel (or two-color) microarrays. 

Unlike conventional estimation methods which take into account only measures of location (e.g., mean and median) of the pixel intensities of each channel, the provided estimation method takes into account pixel-level variability, which may reflects uncertainties due noise and systematic artifacts. 

To remove array-specific effects, intensity-dependent dye biases, and other systematic trends of the microarray data, the M and A values have to be subjected to a within-slide normalization. The most used within-slide normalization technique is LOWESS. However, the choice of the LOWESS parameters, particularly the smoothing neighborhood parameter (or bandwidth), critically affects the quality of the microarray data normalization. Thus, to preserve relevant variation that may be removed in LOWESS normalization with arbitrarily chosen parameters, it is provided a parameter selection method that is parsimonious and considers intrinsic characteristics of microarray data, such as heteroskedasticity.

### Installation

First, install some dependencies from Bioconductor:

```r
source("https://bioconductor.org/biocLite.R")
biocLite("limma")
```

Then, you can download the latest tar.gz file with the source code of the omicsMA R package, available at https://github.com/adele/omicsMA/releases/latest, and install it with the following command, where `path_to_file` represents the full path and file name of the tar.gz file:
```r
install.packages(path_to_file, repos=NULL, type="source", dependencies=TRUE)
```

Or install the development version directly from GitHub. Make sure you have the devtools R package installed. 
If not, install it with `install.packages("devtools")`.

```r
devtools::install_github("adele/omicsMA", dependencies=TRUE)
```

All releases are available at https://github.com/adele/omicsMA/releases. If you want a specific version of the omicsMA R package, for example, v1.0, you can install it directly from the URL:
```r
install.packages("https://github.com/adele/omicsMA/releases/download/v1.0/omicsMA_1.0.tar.gz", repos=NULL, method="libcurl", dependencies=TRUE)
```



