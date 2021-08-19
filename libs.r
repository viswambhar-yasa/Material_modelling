########
# libs.r
# needed R packages
#
# You can run this R script by running source('libs.r').
########

# Load contributed R package 'spatstat', and, if necessary, install it beforehand:
if(!require(spatstat)) {install.packages("spatstat"); require(spatstat)}  

# Load contributed R package 'EBImage', and, if necessary, install it beforehand:
if(!require(EBImage)) {                                                   
  install.packages("BiocManager")
  BiocManager::install("EBImage") 
  require(EBImage)
} 

# Load contributed R package 'bitops', and, if necessary, install it beforehand:
if(!require(bitops)) {install.packages("bitops"); require(bitops)} 
