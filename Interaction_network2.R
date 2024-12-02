#This code needs the data frames made from the RNU4 exclusive R script to function
install.packages("dplyr")
library(dplyr)

install.packages("tidyr")
library(tidyr)

install.packages("igraph")
library(igraph)

install.packages("VennDiagram")
library(VennDiagram)

#Selecting for only the interaction columns and renaming them all P1 or P2
LRep1_Filtered_RNU41 <- LRep1_Filtered_RNU41[, c("Part1", "Part2")]

LRep2_Filtered_RNU41 <- LRep2_Filtered_RNU41[, c("Part1", "Part2")]

RNAInter_RNU41 <- RNAInter_RNU41[,c("Interactor1.Symbol", "Interactor2.Symbol")]
colnames(RNAInter_RNU41) <- c("Part1", "Part2")

NPinter_RNU41 <- NPinter_RNU41[, c("Molecule1", "Molecule2")]
colnames(NPinter_RNU41) <- c("Part1", "Part2")

#All RNU4-1 interactions, All filtered RNU4-1 combined
All_RNU41_interactions <- rbind(LRep1_Filtered_RNU41, LRep2_Filtered_RNU41, NPinter_RNU41, RNAInter_RNU41)

All_RNU41_interactions <- All_RNU41_interactions %>% distinct(Part2, .keep_all = TRUE)


#Selecting for only the interaction columns and renaming them all P1 or P2
LRep1_Filtered_RNU42 <- LRep1_Filtered_RNU42[, c("Part1", "Part2")]

LRep2_Filtered_RNU42 <- LRep2_Filtered_RNU42[, c("Part1", "Part2")]

RNAInter_RNU42 <- RNAInter_RNU42[,c("Interactor1.Symbol", "Interactor2.Symbol")]
colnames(RNAInter_RNU42) <- c("Part1", "Part2")

NPinter_RNU42 <- NPinter_RNU42[, c("Molecule1", "Molecule2")]
colnames(NPinter_RNU42) <- c("Part1", "Part2")

#All RNU4-2 interactions, All filtered RNU4-2 combined
All_RNU42_interactions <- rbind(LRep1_Filtered_RNU42, LRep2_Filtered_RNU42, NPinter_RNU42, RNAInter_RNU42)

All_RNU42_interactions <- All_RNU42_interactions %>% distinct(Part2, .keep_all = TRUE)


#combine
All_RNU4 <- rbind(All_RNU41_interactions, All_RNU42_interactions)

#removing redundant values
All_RNU4 <-All_RNU4[!grepl("U4",All_RNU4$Part2),]
All_RNU4 <-All_RNU4[!grepl("U6",All_RNU4$Part2),]
#Removing RN7SL pseudeogenes
Pseudogenes <- c("575P", "5P", "521P", "660P", "735P", "444P")
All_RNU4 <- All_RNU4[!grepl(paste(Pseudogenes, collapse = "|") ,All_RNU4$Part2),]


#Interaction network graph
g <- graph_from_data_frame(All_RNU4, directed = FALSE)

V(g)$color <- ifelse(V(g)$name == "RNU4-1", "red",
                     ifelse(V(g)$name == "RNU4-2", "green","lightgray"))

V(g)$size <- ifelse(V(g)$name == "RNU4-1", 30,
                    ifelse(V(g)$name == "RNU4-2", 30, 20))



plot(g, layout = layout_with_fr(g), vertex.label = V(g)$name, vertex.label.cex = 0.8,
     vertex.size = V(g)$size, vertex.color = V(g)$color,
     edge.arrow.size = 0.5)

tkplot(g)

#Venn Diagram

venn.plot <- venn.diagram(
  x= list(RNU41 = All_RNU41_interactions, RNU42 = All_RNU42_interactions),
  filename= NULL,
  output = TRUE
)

grid.draw(venn.plot)


