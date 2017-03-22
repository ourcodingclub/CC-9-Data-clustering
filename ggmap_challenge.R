# Use ggmap() to recreate the final cluster map ----
# John Godlee (johngodlee@gmail.com)

## Packages ----
library(ggmap)

## Make a bounding box based on the Lat10 Long10 of sites_membership ----
bbox <- c(min(sites_membership$Long10) - 2,  # Lower left lon
          min(sites_membership$Lat10) - 2,  # Lower left lat
          max(sites_membership$Long10) + 2,  # Upper right lon
          max(sites_membership$Lat10) + 2)  # Upper right lat

## Download map tiles for two different map types ----
google_map <- get_map(location = bbox, source = "google", maptype = "terrain", color = "bw")
stamen_map <- get_map(location = bbox, source = "stamen", maptype = "toner", color = "bw")

## Plot maps using ggplot2 syntax ----
ggmap(google_map) + geom_point(data = sites_membership, 
                               aes(x = Long10, 
                                   y = Lat10, 
                                   colour = cluster_membership, 
                                   shape = cluster_membership),
                               size = 3, alpha = 0.9) + 
  xlab(expression("Longitude ("*degree*")" )) + 
  ylab(expression("Latitude ("*degree*")" )) + 
  theme(legend.title = element_blank())

ggmap(stamen_map) + geom_point(data = sites_membership, 
                               aes(x = Long10, 
                                   y = Lat10, 
                                   colour = cluster_membership, 
                                   shape = cluster_membership), 
                               size = 3, alpha = 0.9) + 
  xlab(expression("Longitude ("*degree*")" )) + 
  ylab(expression("Latitude ("*degree*")" )) + 
  theme(legend.title = element_blank())
