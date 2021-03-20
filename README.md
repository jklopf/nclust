# nclust 

## Installation
Install the package (requires C compiler for installing R package from source):

```
library(devtools)
install_github("pwirapati/nclust")
```
If all is okay, load the package

```
library(nclust)
```

Note that the core algorithm uses OpenMP for parallelization, and will use the default setting for the maximum number of threads. This is system dependent, and there are many ways to change it:

1. Under virtual machine or containers (such as `docker`) this can be set on the whole virtual environment.
2. Environment variables `OMP_THREAD_LIMIT` or `OMP_NUM_THREADS` that has to be set before running R. In Linux, need to be set outside R (such as via `export OMP_NUM_THREADS=4`  before calling R), and doesn't work by calling `Sys.setenv("OMP_NUM_THREADS=4)`from within R. In Windows, Mac or under RStudio, there might be different ways to do it.
3. Using R packages such as `OpenMPController`or [Rompi](https://github.com/pwirapati/Rompi) which can query/set the number of threads from inside R. They change global OpenMP variables for the process, and can be invoked prior to calling clustering function (`nclust` or `coldmap`). These packages require OpenMP runtime library to install. It is not a problem for `gcc`, but can be tricky for `clang` in make. See the README of [Rompi](https://github.com/pwirapati/Rompi).

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
