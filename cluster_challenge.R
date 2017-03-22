# Answer to Challenge to from the Coding Club tutorial on data clustering
# https://ourcodingclub.github.io/2017/03/21/data-clustering.html

# You can get libraries, details on importing data and such from the tutorial
# This script includes just the answer to the second challenge

# Clustering ----
# Group lowlands
lowlands_node <- findMRCA(bol_upgma_cons_nodi, c("AmzBO066", "AmzBO044"))  # Get the node that connects all observations between the interval you set (in this cases, observations between"AmzBO066" and "AmzBO033").
lowlands_tips <- getDescendants(bol_upgma_cons_nodi, lowlands_node)  # Get all nodes and tips linked to the node you have determined above.
lowlands_tips <- na.omit(bol_upgma_cons_nodi$tip.label[lowlands_tips])  # Remove the nodes (NAs) and return just the tips to you.
length(lowlands_tips)  # Count the number of tips you've got for each group. This will be useful to check if the code worked well.
# 118 tips

# Group Andes_Subtropical
Andes_Subtropical_node <- findMRCA(bol_upgma_cons_nodi, c("AndBO034", "AndBO018"))
Andes_Subtropical_tips <- getDescendants(bol_upgma_cons_nodi, Andes_Subtropical_node)
Andes_Subtropical_tips <- na.omit(bol_upgma_cons_nodi$tip.label[Andes_Subtropical_tips])
length(Andes_Subtropical_tips)
# 99 tips

# Summing the values up
118+99

# creating the vector
all_tips2 <- c(lowlands_tips, Andes_Subtropical_tips)  # First we put all tips together
cluster_membership2 <- vector("character",length(all_tips2))  # Then we create a vector as long as the object all_tips
cluster_membership2[which(all_tips2%in% lowlands_tips)] <- "Lowlands"  # Then we assign each set of tips to a subgroup
cluster_membership2[which(all_tips2%in% Andes_Subtropical_tips)] <- "Andes and Subtropical"
length(cluster_membership2)
class(cluster_membership2)

cluster_membership2 <- cluster_membership2[match(sites$AreaCode, all_tips2)]

# Binding "cluster_membership" to "sites"
sites_membership2 <- cbind(sites, cluster_membership2)
unique(sites_membership2$cluster_membership2)  # Checking if all groups are in here
dim(sites_membership2)
head(sites_membership2)

# Making a map ----
map(xlim=c(-70,-55),ylim=c(-25,-8)) # setting the lat,long limits in our map
map.axes()

# Group 1
points(sites_membership2$Long10[which(sites_membership2$cluster_membership2 == "Lowlands")]  # Colour-coding by group.
       ,sites_membership2$Lat10[which(sites_membership2$cluster_membership2 == "Lowlands")], pch = 24, col = rgb
       (t(col2rgb("chartreuse4"))/255, alpha=1), bg = rgb(t(col2rgb("chartreuse4"))/255))

# Group 2
points(sites_membership2$Long10[which(sites_membership2$cluster_membership2 == "Andes and Subtropical")]
       ,sites_membership2$Lat10[which(sites_membership2$cluster_membership2 == "Andes and Subtropical")], pch = "O", col = rgb
       (t(col2rgb("gray53"))/255,alpha=1), bg = rgb(t(col2rgb("gray53"))/255, alpha=1))
