# haverö metabarcoding study - additional R packages
# jesse harrison 2020-2021

# set libpath
.libPaths(c("/home/jharriso/seili-singularity/rpackages", .libPaths()))
libpath <- .libPaths()[1]

# load packages
packages <- c("devtools", "BiocManager")
lapply(packages, require, character.only = TRUE)

# install packages
install.packages("ggfortify", lib = libpath)
install("mixOmics", lib = libpath)
install.packages("RVAideMemoire", lib = libpath)
install.packages("colorspace", lib = libpath)
install.packages("ggrepel", lib = libpath)
install_github("gavinsimpson/ggvegan", lib = libpath)
install_github("jfq3/QsRutils", lib = libpath)
install.packages("magick", lib = libpath)
install.packages("rsvg", lib = libpath)
