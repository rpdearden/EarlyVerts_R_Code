Early vertebrates dispRity pipeline
================
RPDearden
06/11/2018

This is me having a play with the functionalities of dispRity, using the matrix from [Keating *et al.* 2016](http://rspb.royalsocietypublishing.org/content/283/1826/20152917). The tree is a majority rule consensus I obtained by running the dataset in PAUP\*, and the dates are cobbled together from the literature.

### Preliminaries

Read in required packages

``` r
library(dispRity)
library(geomorph)
library(vegan)
library(cluster)
library(Claddis)
```

Set directory and read in nexus file - NOTE ReadMorphNexus reads in both ? and - as NA - instead I have read the matrix in as a csv file (by copying and pasting from Mesquite because I am basic)

``` r
setwd("~/Documents/Side_projects/Masters_writeup/EarlyVerts_R_Code/dispRity_Pipeline")
#Kdata <- ReadMorphNexus("Keating16.nex")
Ktree <- read.nexus("KeatMaj.tre")
Kdates <- read.csv("Keating16Dates.csv", header=TRUE,row.names=1)
Kmatrix <- read.csv("Keating16.csv", header=TRUE, row.names=1, stringsAsFactors = FALSE, na.strings="?") #Reads ? as NA?
```

Change values as detailed in that Deline paper (find)

``` r
#Replace - with -1 (will become 0 later)
Kmatrix[Kmatrix=="-"] <- -1
#Convert all remaining character columns into numeric
indx <- sapply(Kmatrix, is.character)
Kmatrix[indx] <- lapply(Kmatrix[indx], function(Kmatrix) as.numeric(as.character(Kmatrix))) 
 #Add 1 to all values
Kdata <- Kmatrix+1
```

Ordinate - this is the NMDS procedure I used in original project - possible to refine?

``` r
Kdistances=daisy(Kdata, metric=c("gower"))
kdataNMDS<-metaMDS(Kdistances, k=3, zerodist="add")
```

    ## Run 0 stress 0.0604045 
    ## Run 1 stress 0.06528009 
    ## Run 2 stress 0.06040486 
    ## ... Procrustes: rmse 8.933875e-05  max resid 0.0002144964 
    ## ... Similar to previous best
    ## Run 3 stress 0.0604063 
    ## ... Procrustes: rmse 0.000477657  max resid 0.00133977 
    ## ... Similar to previous best
    ## Run 4 stress 0.07103663 
    ## Run 5 stress 0.06538837 
    ## Run 6 stress 0.06040466 
    ## ... Procrustes: rmse 0.0005267043  max resid 0.001238408 
    ## ... Similar to previous best
    ## Run 7 stress 0.06040424 
    ## ... New best solution
    ## ... Procrustes: rmse 0.001450056  max resid 0.00270684 
    ## ... Similar to previous best
    ## Run 8 stress 0.07143277 
    ## Run 9 stress 0.06606869 
    ## Run 10 stress 0.0674164 
    ## Run 11 stress 0.06606065 
    ## Run 12 stress 0.06527616 
    ## Run 13 stress 0.06040474 
    ## ... Procrustes: rmse 0.0003243256  max resid 0.000873928 
    ## ... Similar to previous best
    ## Run 14 stress 0.07135783 
    ## Run 15 stress 0.06040394 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0004252888  max resid 0.001142124 
    ## ... Similar to previous best
    ## Run 16 stress 0.0604056 
    ## ... Procrustes: rmse 0.0008001313  max resid 0.0022685 
    ## ... Similar to previous best
    ## Run 17 stress 0.0604043 
    ## ... Procrustes: rmse 0.0001459228  max resid 0.0002521625 
    ## ... Similar to previous best
    ## Run 18 stress 0.06742327 
    ## Run 19 stress 0.06040459 
    ## ... Procrustes: rmse 0.0005892338  max resid 0.002006138 
    ## ... Similar to previous best
    ## Run 20 stress 0.0713563 
    ## *** Solution reached

``` r
#n.b. Ordination is stored in:
Ordkdata<-kdataNMDS$points
```

Plot first two axes of NMDS - plot looks like the back end of a bus but not much point twiddling currently

``` r
#Normal with text
plot(kdataNMDS)
```

    ## species scores not available

``` r
text(kdataNMDS, cex=0.8)
```

![](dispRity_pipeline_files/figure-markdown_github/unnamed-chunk-5-1.png)

Test stresses

``` r
#Run metaMDS with different numbers of axes
NMDS1=metaMDS(Kdistances, k=1, zerodist="add")
NMDS2=metaMDS(Kdistances, k=2, zerodist="add")
NMDS3=metaMDS(Kdistances, k=3, zerodist="add")
NMDS4=metaMDS(Kdistances, k=4, zerodist="add")
NMDS5=metaMDS(Kdistances, k=5, zerodist="add")
NMDS6=metaMDS(Kdistances, k=6, zerodist="add")
NMDS7=metaMDS(Kdistances, k=7, zerodist="add")
NMDS8=metaMDS(Kdistances, k=8, zerodist="add")
NMDS9=metaMDS(Kdistances, k=9, zerodist="add")
NMDS10=metaMDS(Kdistances, k=10, zerodist="add")

#Make vector containing stresses
Stresses=c(NMDS1$stress, NMDS2$stress, NMDS3$stress, NMDS4$stress, NMDS5$stress, NMDS6$stress, NMDS7$stress, NMDS8$stress, NMDS9$stress, NMDS10$stress)
#Plot it
barplot(Stresses, ylim=c(0, 0.30), names.arg=c(1:10), xlab="No. of axes", ylab="Stress")
```

![](dispRity_pipeline_files/figure-markdown_github/unnamed-chunk-6-1.png)

<!-- Alternatively originally did like below lazily using default Claddis.ordination settings, but have a think about appropriate methods -->
<!-- ```{r warning=FALSE} -->
<!-- OrdKdata <- Claddis.ordination(Kdata) -->
<!-- ``` -->
### Doing disparity through time (disparity metric used the sum of variances as default)

tree requires a root time - for now guessed at start of Cambrian, choose a better one

``` r
Ktree$root.time <- 541
```

Order dates columns to match tip labels in tree (necessary for Claddis{DatePhylo})

``` r
Kdates <- Kdates[match(Ktree$tip.label, rownames(Kdates)),]
```

Date tree using Claddis tree dating - check dating method appropriate

``` r
DatedKtree<-DatePhylo(Ktree, Kdates, method="equal", rlen=1)
```

Make time bins (random 50my)

``` r
time_bins <- rev(seq(from = 250, to = 600, by = 50))
```

Split dataset into subsets

``` r
binned_Kdata <- chrono.subsets(data = Ordkdata, tree = Ktree, method = "discrete", time = time_bins, inc.nodes = FALSE,FADLAD = Kdates)
```

Bootstrapping , first just that and then with bonus rarefaction

``` r
boot_bin_Kdata <- boot.matrix(binned_Kdata)
rare_bin_Kdata <- boot.matrix(binned_Kdata, bootstraps = 100,rarefaction = 6)
```

Calculate disparity for both datasets

``` r
boot_disparity_Kdata <- dispRity(boot_bin_Kdata, metric = c(sum, variances))
rare_disparity_Kdata <- dispRity(rare_bin_Kdata, metric = c(sum, variances))
```

Plot it up

``` r
quartz(width = 15, height = 5) ; par(mfrow = (c(1,2)), bty = "n")
plot(boot_disparity_Kdata, type = "continuous", main = "bootstrapped results")
plot(rare_disparity_Kdata, type = "continuous", main = "rarefied results")
```

![](dispRity_pipeline_files/figure-markdown_github/unnamed-chunk-14-1.png)

<!-- Testing bins difference - won't work with this dataset as not enough data -->
<!-- ```{r warning=FALSE} -->
<!-- #test.dispRity(boot_disparity_Kdata, test = wilcox.test, comparison = "sequential",correction = "bonferroni") -->
<!-- ``` -->
