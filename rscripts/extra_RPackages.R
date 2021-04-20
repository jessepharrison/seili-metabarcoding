# install a few extra R packages

# set libpath
.libPaths(c("/home/jharriso/seili-singularity/rpackages", .libPaths()))
libpath <- .libPaths()[1]

# install packages
install.packages("ggfortify", lib = libpath)
BiocManager::install("mixOmics", lib = libpath, update = FALSE)
install.packages("RVAideMemoire", lib = libpath)

# DESeq2 wants some packages installed separately
install.packages("XML", lib = libpath)
BiocManager::install("annotate", lib = libpath, update = FALSE)
BiocManager::install("genefilter", lib = libpath, update = FALSE)
BiocManager::install("geneplotter", lib = libpath, update = FALSE)
BiocManager::install("DESeq2", lib = libpath, update = FALSE)

# Other packages
install.packages("ggrepel", lib = libpath)
devtools::install_github("gavinsimpson/ggvegan", lib = libpath)
devtools::install_github("jfq3/QsRutils", lib = libpath)
install.packages("magick", lib = libpath)
install.packages("rsvg", lib = libpath)

# QsRutils is from:
# https://github.com/jfq3/QsRutils
