# nclust 

## Installation
### Linux

Install the package (requires C compiler for installing R package from source):

```
library(devtools)
install_github("pwirapati/nclust")
```
If all is okay, load the package

```
library(nclust)
```

### MacOS

Similar to Linux when setup to use `gcc` for compilation. With `clang` it should also work, but see [https://github.com/pwirapat/Rompi] when
controling OpenMP at runtime from inside R is needed (OpenMP shared library is needed).

### Windows

In order to successfuly install nclust from source, Windows users need to have rtools and the BLAS library installed on their system.

Install rtools to build R packages: [rtools40](https://cran.r-project.org/bin/windows/Rtools/)

If you have rtools installed, download the x32 and x64 prebuild BLAS library from the Rtools packages bintray:
- x32: https://dl.bintray.com/rtools/mingw32/mingw-w64-i686-lapack-3.8.0-2-any.pkg.tar.xz
- x64: https://dl.bintray.com/rtools/mingw64/mingw-w64-x86_64-lapack-3.8.0-2-any.pkg.tar.xz

Open a rtools mingw32 shell and run:
```
pacman -Su
pacman -U mingw-w64-i686-lapack-3.8.0-2-any.pkg.tar.xz
```
Open a rtools mingw64 shell and run:
```
pacman -Su
pacman -U mingw-w64-x86_64-lapack-3.8.0-2-any.pkg.tar.xz
```
After successfull BLAS (LAPACK) library installation, install nclust from source:
```
library(devtools)
install_github("pwirapati/nclust")
```

If all is okay, load the package

```
library(nclust)
```

## Demo

We generate a random data with some structure. Just run the code below verbatim.

```
set.seed(666)
design <- data.frame(status=c(rep("disease",50),rep("normal",50)),batch=paste0("batch_",rbinom(100,1,.5)+1),row.names=paste0("patient_",1:100))

m <- 1000   # number of genes
n <- 50 + 50
batcheffect <- cbind("batch_1"=rnorm(m),"batch_2"=rnorm(m))

x <- array(rnorm(m*n),dim=c(m,n)) + 0.5 * batcheffect[,design$batch]
colnames(x) <- rownames(design)
rownames(x) <- paste0("gene_",1:m)
x[1,] <- x[1,] + ifelse(design$status=="disease",2,-2)
```

This creates two main objects:

* `x` is the data matrix (with random normal error, systematic batch effects and the first gene as the disease marker )
* `design` is the sample annotation.

To do the clustering of both rows and columns:

```
myclust <- coldmap(x)
```

This perform the procedure, save the clusters in `myclust` and display the heatmap.

To show the information from the design table:

```R
coldmap(x, clust=myclust, ctag=make_tag( design, varnames=c("status","batch"),cols=c("violet","green3")), ctag.space=3, rmarg=3 )
```

See `?coldmap` for the options descriptions.

To show labels:

```R
coldmap(x, clust=myclust, ctag=make_tag( design, varnames=c("status","batch"),cols=c("violet","green3")),
ctag.space=3, rmarg=3,
       rlab=list("gene_1",200:205),clab=list(c(1,50,n),c("patient_13","patient_33"df)))
```


## Gene Expression Data

`coldmap` can be applied to each stage in the preprocessing pipeline, to visualize the effect on various steps on the pattern and the clustering.

Typically, these are needed to remove the obvious and to highlight biological (or batch/confounding) effects:

* Transform to log scale
* Normalize
* Gene centering

Column and/or row centering can be done by functions provided in `nclust` (see `?center`).
