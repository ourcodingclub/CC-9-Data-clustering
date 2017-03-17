#Data clustering tutorial
#Clustering data for Bolivia (if Ary allows it)

#Installing packages needed for the analysis
install.packages("recluster")
install.packages("phytools")
install.packages("maps")
library(recluster)
library(phytools)
library(maps)

#Importing our data frames
#species
spp <- read.csv("spp_bol.csv", sep=",", head=TRUE)
head(spp)
dim(spp)

#sites
sites_bol <- read.csv("sites_bolivia.csv", sep=",", head=TRUE) 
dim(sites_bol)
head(sites_bol)

#sppxarea matrix
sppxsites <- read.csv("sppxsites_bol.csv", sep=",", head=TRUE)
dim(sppxsites)
head(sppxsites)
tmp <- unique(sppxsites$AreaID)
length(tmp)
tmp <- unique(sppxsites$SppID)
length(tmp)

#making the species by site matrix (presence and abscence). We'll call it commat
sites_sub <- unique(sppxsites$AreaID)
spp_sub <- unique(sppxsites$SppID)

spp_commat <- matrix(0,length(sites_sub),length(spp_sub))
for (i in 1:nrow(spp_commat)){
  temp_sites <- sppxsites[which(sppxsites$AreaID==sites_sub[i]),]
  spp_commat[i,which(spp_sub%in%temp_sites$SppID)] <- 1
  print(i)
}

rownames(spp_commat) <- as.character(sites$AreaCode[match(sites_sub,sites$AreaID)])
colnames(spp_commat) <- as.character(spp$Species.code[match(spp_sub,spp$SppID)])
dim(spp_commat)
spp_commat[1:6,1:6]

#Removing singletons (species with only one occurrence record on the dataset. They tend to generate noise in the analysis
#without adding more information to the results).

spp_commat_trim <- spp_commat[,which(!colSums(spp_commat) == 1)]
dim(spp_commat_trim)

#Building different clusters with different methods
#recluster uses simpson (best distance metric EVER!) by default
#Why not compare clusters made with all species present in the dataset and the one without singletons

#the write.tree function will let you save your cluster in TRE format so you can open it with figtree and visualize it.

#single-linkage method
#full spp presence x abscence matrix
bol_singlelink <- recluster.cons(spp_commat, tr=100, p=0.5, method = "single")
bol_singlelink_cons <- bol_singlelink$cons
write.tree(bol_singlelink_cons, "bol_singlelink_cons.tre")
plot(bol_singlelink_cons, direction = "downwards")

#trimmed matrix
bol_singlelink_trim <- recluster.cons(spp_commat_trim, tr = 100, p = 0.5, method = "single")
bol_singlelink_trim_cons <- bol_singlelink_trim$cons
write.tree(bol_singlelink_trim_cons, "bol_singlelink_cons_trim.tre")
plot(bol_singlelink_trim_cons, direction = "downwards")

#complete-linkage method
#full spp presence x abscence matrix
bol_completelink <- recluster.cons(spp_commat, tr=100, p=0.5, method = "complete")
bol_completelink_cons <- bol_completelink$cons
write.tree(bol_completelink_cons, "bol_completelink_cons.tre")
plot(bol_completelink_cons, direction = "downwards")

#trimmed matrix
bol_completelink_trim <- recluster.cons(spp_commat_trim, tr=100, p=0.5, method = "complete")
bol_completelink_cons_trim <- bol_completelink_trim$cons
write.tree(bol_completelink_cons_trim, "bol_completelink_cons_trim.tre")
plot(bol_completelink_cons_trim, direction = "downwards")

#clustering using ward's minimum variance
#full spp presence x abscence matrix
bol_ward <- recluster.cons(spp_commat, tr=100, p=0.5, method = "ward.D")
bol_ward_cons <- bol_ward$cons
plot(bol_ward_cons, direction = "downwards")
write.tree(bol_ward_cons, "bol_ward_cons.tre")

#trimmed matrix
bol_ward_trim <- recluster.cons(spp_commat_trim, tr=100, p=0.5, method = "ward.D")
bol_ward_cons_trim <- bol_ward_trim$cons
plot(bol_ward_cons_trim, direction = "downwards")
write.tree(bol_ward_cons_trim, "bol_ward_cons_trim.tre")

#average method (UPGMA)
#full spp presence x abscence matrix
bol_upgma <- recluster.cons(spp_commat, tr=100, p=0.5, method = "average")
bol_upgma_cons <- bol_upgma$cons
write.tree(bol_upgma_cons, "bol_upgma_cons.tre")
plot(bol_upgma_cons, direction = "downwards")

#trimmed spp presence x abscence matrix
bol_upgma_trim <- recluster.cons(spp_commat_trim, tr=100, p=0.5, method = "average")
bol_upgma_trim_cons <- bol_upgma$cons
write.tree(bol_upgma_trim_cons, "bol_upgma_cons_trim.tre")
plot(bol_upgma_trim_cons, direction = "downwards")

#You can calculate support values for your branches
#Bootstrap is one of the main forms of calculating support values, but the values tend to decrease when you add more data
#(our case here). So don't get frustrated by the low values we got as bootstrap is not exactly useful in our case.
#There are other ways of testing cluster consistency and we'll cover that soon.
#the "boot" argument tells the function how many radomizations you want the command to make before estimating the bootstrap
#values. The bigger the number, the longer it will take for your computer to run it. The standard number is 1000, but we'll
#try with only 10 for now.

#Full matrix
bol_upgma_boot <- recluster.boot(bol_upgma_cons, spp_commat, tr = 100, p = 0.5, method = "average", boot = 1000, level = 1)
recluster.plot(bol_upgma_cons, bol_upgma_boot, direction = "downwards")

#Trimmed matrix
bol_upgma_trim_boot <- recluster.boot(bol_upgma_cons, spp_commat_trim, tr = 100, p = 0.5, method = "average", boot = 1000, level = 1)
recluster.plot(bol_upgma_cons, bol_upgma_boot, direction = "downwards")

#Now, for the sake of brevity, let us work with the UPGMA cluster created with the full dataset
#Since these are sites located in geographic space, it would be interesting to plot groups you recognized on your cluster and
#check how they are distributed in multivariate and geographic spaces. However, for you to do that, you need to create a
#vector with the groups you want to explore, sort of a cluster membership. For that, we'll use tools that are usually
#employed by people studying phylogenies. Phylogenies can sound difficult and out of reach, but there's no need to know
#what they are exactly in order to use the tools we'll be using here.

#first of all, we have to make R acknowledge the existence of politomies, which are nodes that lead to more than two nodes/
#tips and/or groups, like the ones we have found on the previous clusters. If you type the name of your consensus tree on 
#your console, you get a quick summary about your cluster (it says that it is a phylogenetic tree, but it's not. Keep that
#in mind) and how many tips (terminal branches) and nodes (ramifications) it has. A cluster with no polytomies will have 
#x tips and x-1 internal nodes.

bol_upgma_cons_nodi <- di2multi(bol_upgma_cons)

#This cluster in particulas has no politomies, but others might have. Always use this command, even though you are pretty
#sure you don't have any politomies. Depending on how big your cluster is, it is really hard to assess this visualy.

#As an example, let's take a look at thow this command will work on a cluster with actual politomies

bol_completelink_cons_nodi <- di2multi(bol_completelink_cons)

#Now we got a reduction in the number of nodes

sort(bol_upgma_cons$tip.label)

#Now that we definetely have a politomy free cluster, we can make that vector of cluster memberships in order to map
#our sites and help us understand how our data is organized in spatial and compositional space.

#Group 1 tips
group1_node <- findMRCA(bol_upgma_cons, c("AmzBO066","AmzBO033"))
group1_tips <- getDescendants(bol_upgma_cons,group1_node)
group1_tips <- na.omit(bol_upgma_cons$tip.label[group1_tips])
length(group1_tips)
#90 tips

#Group 2 tips
group2_node <- findMRCA(bol_upgma_cons, c("CerBO012","AmzBO044"))
group2_tips <- getDescendants(bol_upgma_cons,group2_node)
group2_tips <- na.omit(bol_upgma_cons$tip.label[group2_tips])
length(group2_tips)
#28 tips

#Group 3 tips
group3_node <- findMRCA(bol_upgma_cons, c("AndBO035","CerBO018"))
group3_tips <- getDescendants(bol_upgma_cons,group3_node)
group3_tips <- na.omit(bol_upgma_cons$tip.label[group3_tips])
length(group3_tips)
#52 tips

#Group 4 tips
group4_node <- findMRCA(bol_upgma_cons, c("AndBO046","AndBO018"))
group4_tips <- getDescendants(bol_upgma_cons,group4_node)
group4_tips <- na.omit(bol_upgma_cons$tip.label[group4_tips])
length(group4_tips)
#47

#Now that we have the tips, we need to check if no one was left behind. This is done by siply summing the amount of tips we got for
#each group and seeing if we get all 217 tips that we have on our cluster
90+28+52+47

#There you go. Everyone is here. We good to go.
#Now let's make that vector and bind it to our dataset.

#all_tips
all_tips <- c(group1_tips, group2_tips, group3_tips, group4_tips)
cluster_membership <- vector("character",length(all_tips))
cluster_membership[which(all_tips%in% group1_tips)] <- "Group 1"
cluster_membership[which(all_tips%in% group2_tips)] <- "Group 2"
cluster_membership[which(all_tips%in% group3_tips)] <- "Group 3"
cluster_membership[which(all_tips%in% group4_tips)] <- "Group 4"
length(cluster_membership)
class(cluster_membership)

#matching cluster_membership with areas_lowtrop_1000m_nofrost
cluster_membership <- cluster_membership[match(sites_bol$AreaCode,all_tips)]

#binding cluster membership to areas_lowtrop_1000m_notemperate_norogue
sites_bol_membership <- cbind(sites_bol, cluster_membership)
unique(sites_bol_membership$cluster_membership)
dim(sites_bol_membership)
head(sites_bol_membership)
unique(sites_bol_membership$cluster_membership)

#Let's make the map now

map(xlim=c(-70,-55),ylim=c(-25,-8))
map.axes()
#Atlantic group
points(sites_bol_membership$Long10[which(sites_bol_membership$cluster_membership == "Group 1")]
       ,sites_bol_membership$Lat10[which(sites_bol_membership$cluster_membership == "Group 1")],pch=24,col=rgb
       (t(col2rgb("chartreuse4"))/255,alpha=1), bg=rgb(t(col2rgb("chartreuse4"))/255))
#Cerrado group
points(sites_bol_membership$Long10[which(sites_bol_membership$cluster_membership == "Group 2")]
       ,sites_bol_membership$Lat10[which(sites_bol_membership$cluster_membership == "Group 2")],pch="O",col=rgb
       (t(col2rgb("gray53"))/255,alpha=1), bg=rgb(t(col2rgb("gray53"))/255, alpha=1))
#Amazon group
points(sites_bol_membership$Long10[which(sites_bol_membership$cluster_membership == "Group 3")]
       ,sites_bol_membership$Lat10[which(sites_bol_membership$cluster_membership == "Group 3")],pch=15,col=rgb
       (t(col2rgb("blue"))/255,alpha=1), bg=rgb(t(col2rgb("blue"))/255, alpha=1))
#SDTF group
points(sites_bol_membership$Long10[which(sites_bol_membership$cluster_membership == "Group 4")]
       ,sites_bol_membership$Lat10[which(sites_bol_membership$cluster_membership == "Group 4")],pch=19,col=rgb
       (t(col2rgb("saddlebrown"))/255,alpha=1), bg=rgb(t(col2rgb("saddlebrown"))/255))

#Very interesting results!



























