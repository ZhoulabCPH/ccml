# CCML: A two-step consensus clustering inputing multiple predictive labels with different sample coverages (missing labels)

## Dependencies
------------
* diceR
* ggplot2
* parallel
* plyr
* tidyr
* SNFtools
* ConsensusClusterPlus



## Installation
------------
Get the latest released version from CRAN:    
``` r
install.packages("ccml")
```
Or from GitHub:
``` r
library(devtools)
devtools::install_github("ZhouSunLab-Workshops/ccml")
```
Or install locally：
``` r
install.packages("ccml_1.4.0.tar.gz",repos = NULL, type="source")
```

## Setup
------------
This package needs to create a folder to store the result. Let’s call the folder
`result_output` and assume the full path to this folder is
`xxxxx/result_output`. Run the function below to set the path for the package.

``` r
library(ccml)
setwd("xxxxx/result_output")
```

Or a folder named "result_output" will be created under the default path when the program is running.



## Introduction
------------
Consensus clustering, also called meta-clustering or cluster ensembles, has been increasingly used in clinicadata. Current consensus clustering methods tend to ensemble a number of different clusters frommathematical replicates with similar sample coverage. As the fact of common variety of sample coverage inthe real-world data, a new consensus clustering strategy dealing with such biological replicates is required.This is a two-step consensus clustering package, which is used to input multiple predictive labels withdifferent sample coverage (missing labels).
Ccml consists of three key steps. The first is calculating normalized consensus weight matrix throughthousands of times of random permutation. Then, compare the difference between original and normalizedconsensus weights. The last step is to use ConsensusClusterPlus to cluster based on normalized consensusweights.



## Functions
------------
* cml     
*A two-step consensus clustering inputing multiple predictive labels with different sample coverages (missing labels).*

* callNCW  
*Calculate normalized consensus weight(NCW) matrix based on permutation.*


* plotCompareCW  
*Plot of original consensus weights vs. normalized consensus weights grouping by the number of co-appeared percent of clustering(non-missing).*


* randConsensusMatrix   
*Calculate consensus weight matrix based on the permuted input label matrix. Internal function used by callNCW.*


* spectralClusteringAffinity  
*Perform spectral clustering algorithms for an affinity matrix, using SNFtool::spectralClustering.*




## Run examples
------------
In the repository, we give some examples to show how to run **ccml**.
*example_data is built into R package and can be loaded with data().*

* load data   
*A matrix or data frame of input labels or a character value of input file name, columns are different clustering results and rows are samples.*
  ``` r
  data(example_data)
  label=example_data
  ```

* output directory      
*This part is mentioned in Setup above.*
  ``` r
  title="result_output"
  ```
  
* run cml
   ```r
    # not estimate stability of permutation numbers.
    res_1<-ccml(title=title,label=label,nperm = 3,ncore=1,stability=FALSE,maxK=5,pItem=0.8)

    # other methods for clustering of distance matrix
    res_2<-ccml(title=title,label=label,nperm = 10,ncore=1,stability=TRUE,
           maxK=3,pItem=0.9,clusterAlg = "hc")

    # not output as "rdata"
    res_3<-ccml(title=title,label=label,output=FALSE,nperm = 5,ncore=1,seedn=150,
           stability=TRUE,maxK=3,pItem=0.9)

    # The output is a list.Some examples are as follws
    # get consensusMatrix
    res$fcluster[[3]]$consensusMatrix

    # get consensusTree
    res$fcluster[[2]]$consensusTree

    # get consensusClass
    res$fcluster[[2]]$consensusClass

   ```

* run callNCW
   ```r
    ncw<-callNCW(title=title,label=label,stability=TRUE,nperm=4,ncore=1)
   ```
   
* run plotCompareCW
   ```r
    plotCompareCW(title, label, ncw, plot = NULL)
   ```
  

## Output Description
------------
* consensus     
  The heatmaps of the consensus matrices(for k= 2, 3, 4),consensus CDF Plot, Delta Area Plot and Tracking Plot are used to determine the best K value of clustering.
  
* icl  
  Item-Consensus Plot show the stability of members.   
  
* plotCompareCW   
  This plot show the original consensus weights vs. normalized consensus weights grouping by the number of co-appeared percent of clustering(non-missing).The number of duplicates sample are indicated by the size of the colored portion, whose color corresponds to percent of co-appeared of clustering.   

* stability.seed100.n10K   
  The name of plot represents 10*1000 permutation using random numbers with seedn=100.The stability of normalized consensus weight is estimated based on permutation numbers.
