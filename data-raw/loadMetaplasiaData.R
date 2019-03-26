library(maigesPack)

# to generate data:
# setwd("data-raw/")
# source("loadMetaplasiaData.R")

##########################################################
# Loading metaplasia raw dataset as a Maiges Pack object #
##########################################################

# To load the gastro data, we used the maigesPack R package:
# http://www.bioconductor.org/packages/release/bioc/html/maigesPack.html
# To install maigesPack, uncomment the two following lines:
# source("https://bioconductor.org/biocLite.R")
# biocLite("maigesPack")

# maigesPack object with the metaplasia dataset
# The object has 45015 spots and 10 samples
metaplasiaMP <- loadData(fileConf="loadSelMetaplasiaData.conf")


###########################################
# Processing Maiges Pack object as a list #
###########################################


# Ids of the n slides:
barcodes <- getLabels(metaplasiaMP, "BARCODE")
squares <- getLabels(metaplasiaMP, "QUADRANTE")
array.ids <- paste(barcodes, "_", squares, sep="")

#> array.ids (firts 5 are metaplasia and the last 5 are normal)
# [1] "251485069613_1.1" "251485069611_1.2" "251485069613_1.4" "251485069614_1.4"
# [5] "251485069514_1.2" "251485069484_1.3" "251485069611_1.4" "251485069511_1.1"
# [9] "251485069395_1.1" "251485069395_1.4"


# Ids of the p spots/genes:
gene.ids <- getLabels(metaplasiaMP, "ID", sLabel=FALSE)
gene.names <- getLabels(metaplasiaMP, "GeneName", sLabel=FALSE)

# defining the p x n matrices
# with the pixel-level measurements
R.mean <- metaplasiaMP@Data$rMeanSignal
G.mean <- metaplasiaMP@Data$gMeanSignal
R.bckg <- metaplasiaMP@Data$rBGMedianSignal
G.bckg <- metaplasiaMP@Data$gBGMedianSignal
R.size <- metaplasiaMP@Data$rNumPix
G.size <- metaplasiaMP@Data$gNumPix
R.var <- metaplasiaMP@Data$rPixSDev^2
G.var <- metaplasiaMP@Data$gPixSDev^2
RG.cov <- metaplasiaMP@Data$PixCorrelation * metaplasiaMP@Data$rPixSDev * metaplasiaMP@Data$gPixSDev


metaplasia <- list(array.ids = array.ids, gene.ids = gene.ids, gene.names = gene.names,
                   R.mean = R.mean, G.mean = G.mean,
                   R.var = R.var, G.var = G.var, RG.cov = RG.cov,
                   R.bckg = R.bckg, G.bckg = G.bckg,
                   R.size = R.size, G.size = G.size)

devtools::use_data(metaplasia, compress="xz")
