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
x[1,] <- x[1,] + ifelse(design$status=="disease",2,-2)
```

This creates two main objects:

* `x` is the data matrix (with random normal error, systematic batch effects and the first gene as the disease marker )
* `design` is the sample annotation.

To do the clustering of both rows and columns:

```
myclust <- coldmap(x)
```

This perform the procedure, save the clusters in `myclust` and display the heatmap. We see two main column clusters.

To show the information in the design:

```
coldmap(x, clust=myclust, ctag=make_tag( design, varnames=c("status","batch"),cols=c("violet","green3")), ctag.space=3, rmarg=3 )
```

See `?coldmap` for the options descriptions.
