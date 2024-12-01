#This code needs the data frames made from the RNU4 exclusive R script to function
install.packages("dplyr")
library(dplyr)

install.packages("tidyr")
library(tidyr)

install.packages("igraph")
library(igraph)


#Interactions that are RNU4-1 Exclusive                                                                       
NW_LRep1_In_RNU41_only <- anti_join(LRep1_Filtered_RNU41, LRep1_Filtered_RNU42, by='Part2')

NW_LRep2_In_RNU41_only <- anti_join(LRep2_Filtered_RNU42, LRep2_Filtered_RNU41, by='Part2')

NW_RNAInter_In_RNU41_only <- anti_join(RNAInter_RNU41, RNAInter_RNU42, by='Interactor2.Symbol')

NW_NPinter_RNU41_only <- anti_join(NPinter_RNU41, NPinter_RNU42, by='Molecule2')


#Selecting for only the interaction columns and renaming them all P1 or P2
NW_LRep1_In_RNU41_only <- NW_LRep1_In_RNU41_only[, c("Part1", "Part2")]

NW_LRep2_In_RNU41_only <- NW_LRep2_In_RNU41_only[, c("Part1", "Part2")]

NW_RNAInter_In_RNU41_only <- NW_RNAInter_In_RNU41_only[,c("Interactor1.Symbol", "Interactor2.Symbol")]
colnames(NW_RNAInter_In_RNU41_only) <- c("Part1", "Part2")

NW_NPinter_RNU41_only <- NW_NPinter_RNU41_only[, c("Molecule1", "Molecule2")]
colnames(NW_NPinter_RNU41_only) <- c("Part1", "Part2")

#combine RNU41 and then remove duplicates
All_RNU41_exclusive <- rbind(NW_LRep1_In_RNU41_only, NW_LRep2_In_RNU41_only, NW_NPinter_RNU41_only, NW_RNAInter_In_RNU41_only)

All_RNU41_exclusive <- All_RNU41_exclusive %>% distinct(Part2, .keep_all = TRUE)



#selecting for only the interaction columns for RNU42
NW_LRep1_In_RNU42_only <- LRep1_In_RNU42_only[, c("Part1", "Part2")]

NW_LRep2_In_RNU42_only <- LRep2_In_RNU42_only[, c("Part1", "Part2")]

NW_RNAInter_RNU42_only <- RNAInter_RNU42_only[,c("Interactor1.Symbol", "Interactor2.Symbol")]
colnames(NW_RNAInter_RNU42_only) <-c("Part1", "Part2")

NW_NPinter_RNU42_only <- NPinter_RNU42_only[, c("Molecule1", "Molecule2")]
colnames(NW_NPinter_RNU42_only) <- c("Part1", "Part2")


#combine RNU42
All_RNU42_exclusive <- rbind(NW_LRep1_In_RNU42_only, NW_LRep2_In_RNU42_only, NW_RNAInter_RNU42_only, NW_NPinter_RNU42_only)

All_RNU42_exclusive <- All_RNU42_exclusive %>% distinct(Part2, .keep_all = TRUE)



#LIGR rep 1, All RNU4 interactions separated
NW_LRep1_Filtered_RNU41 <- LRep1_Filtered_RNU41[, c("Part1", "Part2")]
NW_LRep1_Filtered_RNU42 <- LRep1_Filtered_RNU42[, c("Part1", "Part2")]

##Remove duplicates
NW_LRep1_Filtered_RNU41 <- NW_LRep1_Filtered_RNU41 %>% distinct(Part2, .keep_all = TRUE)
NW_LRep1_Filtered_RNU42 <- NW_LRep1_Filtered_RNU42 %>% distinct(Part2, .keep_all = TRUE)

# Rep 2, All RNu4 interactions separated
NW_LRep2_Filtered_RNU41 <- LRep2_Filtered_RNU41[, c("Part1", "Part2")]
NW_LRep2_Filtered_RNU42 <- LRep2_Filtered_RNU42[, c("Part1", "Part2")]

##Remove duplicates
NW_LRep2_Filtered_RNU41 <- NW_LRep2_Filtered_RNU41 %>% distinct(Part2, .keep_all = TRUE)
NW_LRep2_Filtered_RNU42 <- NW_LRep2_Filtered_RNU42 %>% distinct(Part2, .keep_all = TRUE)

#RNAinter, ALl RNU4 interaction separated
NW_RNAInter_RNU41 <- RNAInter_RNU41[,c("Interactor1.Symbol", "Interactor2.Symbol")]
colnames(NW_RNAInter_RNU41) <-c("Part1", "Part2")
NW_RNAInter_RNU42 <- RNAInter_RNU42[,c("Interactor1.Symbol", "Interactor2.Symbol")]
colnames(NW_RNAInter_RNU42) <-c("Part1", "Part2")

##Remove duplicates
NW_RNAInter_RNU41 <- NW_RNAInter_RNU41 %>% distinct(Part2, .keep_all = TRUE)
NW_RNAInter_RNU42 <- NW_RNAInter_RNU42 %>% distinct(Part2, .keep_all = TRUE)

#NPinter, All RNU4 interactions separated
NW_NPinter_RNU41 <- NPinter_RNU41[, c("Molecule1", "Molecule2")]
colnames(NW_NPinter_RNU41) <- c("Part1", "Part2")
NW_NPinter_RNU42 <- NPinter_RNU42[, c("Molecule1", "Molecule2")]
colnames(NW_NPinter_RNU42) <- c("Part1", "Part2")

##Remove duplicates
NW_NPinter_RNU41 <- NW_NPinter_RNU41 %>% distinct(Part2, .keep_all = TRUE)
NW_NPinter_RNU42 <- NW_NPinter_RNU42 %>% distinct(Part2, .keep_all = TRUE)

#Combine shared

Shared_RNU4<- rbind(NW_LRep1_Filtered_RNU41, NW_LRep1_Filtered_RNU42, NW_LRep2_Filtered_RNU41, 
                    NW_LRep2_Filtered_RNU42, NW_RNAInter_RNU41, NW_RNAInter_RNU42, NW_NPinter_RNU41, NW_NPinter_RNU42)

##Select for duplicates in combined
Shared_RNU4 <- Shared_RNU4[duplicated(Shared_RNU4) | duplicated(Shared_RNU4, fromLast = TRUE), ]

#combine the RNU4-1 and -2 and shared
All_RNU4 <- rbind(All_RNU41_exclusive, All_RNU42_exclusive, Shared_RNU4)

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
