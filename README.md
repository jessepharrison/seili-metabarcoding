# Haverö metabarcoding study

## Publication reference

This repository contains material supporting the following publication:

Harrison JP, Chronopoulou P-M, Salonen IS, Jilbert T, Koho KA. (2021) [16S and 18S rRNA gene metabarcoding provide congruent information on the responses of sediment communities to eutrophication](https://www.frontiersin.org/articles/10.3389/fmars.2021.708716/full). Frontiers in Marine Science 8:708716.

---

## Repository contents

The contents include:

- R scripts used for data analysis and package installations (`/rscripts`)
- phyloseq objects including data processed in micca (`/rdata`)
- environmental metadata (`/envdata`)
- a [Singularity](https://sylabs.io/singularity) definition file and container launch script (`/singularity`)
- output files generated by R scripts (`/figures`, `/stats`)

Descriptions of R scripts:

`PCA.R`
Script for environmental metadata PCA

`16S_singularity_datachecks.R` and `18S_singularity_datachecks.R`:
Data analysis scripts for inspecting data pre-processing steps

`16S_singularity.R` and `18S_singularity.R`:
Data analysis scripts generating most output files 

`taxa_summary.R`:
Phyloseq supporting script, originally available [here](https://github.com/joey711/phyloseq/issues/818)

`imagestitch.R`
Script for combining and editing figure files

`extra_RPackages.R`
Script for additional R package installations

---

## Reproducing the R environment + analyses used

### 1. Clone the repository and move to the `singularity` folder

```
git clone https://github.com/jessepharrison/seili-metabarcoding.git
cd seili-metabarcoding/singularity
```

### 2. Build the Singularity container

This takes some time - grab yourself a coffee! 

```
sudo singularity build seili.sif container.def
```

### 3. Install further R packages

Run the following still in the `singularity` folder:

```
singularity exec seili.sif \
Rscript --no-save ../rscripts/extra_RPackages.R

# note 1:
# the script will install the packages to the current working directory.
# the actual data analysis scripts add /singularity to libpaths.
# hence it's important to make sure you install everything to /singularity!

# note 2:
# the container is configured to use a MRAN snapshot from 2020-11-29.  
```

### 4. Running analysis scripts

You can now move back to the repository root (`path/to/seili-metabarcoding`).

The scripts can be run non-interactively using `Rscript`, for example:

```
singularity exec singularity/seili.sif \
Rscript --no-save rscripts/PCA.R
```

The `imagestitch.R` script should be run last as it relies on figure
files created by other scripts.

c) If wanting to work with RStudio:

```
unset XDG_RUNTIME_DIR
singularity exec singularity/seili.sif \
rstudio
```
