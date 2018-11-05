**dispRity pipeline**

**Disparity analysis for early vertebrates, example using Keating *et al.* 2016 and the Github release of dispRity**
***
```r
###read in packages

library(dispRity)
library(geomorph)
library(vegan)
library(cluster)
library(Claddis)


###set directory and read in nexus file

setwd("~/Documents/Side_projects/Masters_writeup/EarlyVerts_R_Code/dispRity_Pipeline")
Kdata <- ReadMorphNexus("Keating16.nex")
Ktree <- read.nexus("KeatMaj.tre")
Kdates <- read.csv("Keating16Dates.csv", header=TRUE,row.names=1)



###Ordinate - currently lazily using default Claddis.ordination settings, but have a think

OrdKdata <- Claddis.ordination(Kdata)



###Doing disparity through time (disparity metric used the sum of variances as default)###
##Worked dispRity example##

#tree needs a root time - for now guessed at start of Cambrian, choose a better one
Ktree$root.time <- 541

#Order dates columns to match tip labels in tree
Kdates <- Kdates[match(Ktree$tip.label, rownames(Kdates)),]

#Date tree using Claddis tree dating - check appropriate
DatedKtree<-DatePhylo(Ktree, Kdates, method="equal", rlen=1)

#Make time bins (random 50my)
time_bins <- rev(seq(from = 250, to = 600, by = 50))

#Split dataset into subsets
binned_Kdata <- chrono.subsets(data = OrdKdata, tree = Ktree, method = "discrete", time = time_bins, inc.nodes = FALSE,FADLAD = Kdates)

#Bootstrapping 
(boot_bin_Kdata <- boot.matrix(binned_Kdata))

# Bootstrapping with bonus rarefaction
rare_bin_Kdata <- boot.matrix(binned_Kdata, bootstraps = 100,
    rarefaction = 6)

#Calculate disparity for bootstrapped data
boot_disparity_Kdata <- dispRity(boot_bin_Kdata, metric = c(sum, variances))

#Calculate disparity for rarefied data
rare_disparity_Kdata <- dispRity(rare_bin_Kdata, metric = c(sum, variances))

#Plot it up
quartz(width = 10, height = 5) ; par(mfrow = (c(1,2)), bty = "n")
plot(boot_disparity_Kdata, type = "continuous", main = "bootstrapped results")
plot(rare_disparity_Kdata, type = "continuous", main = "rarefied results")

##Testing bins difference - won't work as not enough data
test.dispRity(boot_disparity_Kdata, test = wilcox.test, comparison = "sequential",
    correction = "bonferroni")

```