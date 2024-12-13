#Selecting for only the interaction columns and renaming them all P1 or P2
LRep1_Filtered_RNU41 <- LRep1_Filtered_RNU41[, c("Part1", "Part2")]

LRep2_Filtered_RNU41 <- LRep2_Filtered_RNU41[, c("Part1", "Part2")]

RNAInter_RNU41 <- RNAInter_RNU41[, c("Interactor1.Symbol", "Interactor2.Symbol")]
colnames(RNAInter_RNU41) <- c("Part1", "Part2")

NPinter_RNU41 <- NPinter_RNU41[, c("Molecule1", "Molecule2")]
colnames(NPinter_RNU41) <- c("Part1", "Part2")


All_RNU41_interactions <- rbind(LRep1_Filtered_RNU41, LRep2_Filtered_RNU41, NPinter_RNU41, RNAInter_RNU41)

RNU41_interaction_counts <- All_RNU41_interactions %>%
  group_by(Part1, Part2) %>%
  summarize(count = n(), .groups = "drop")


#Selecting for only the interaction columns and renaming them all P1 or P2
LRep1_Filtered_RNU42 <- LRep1_Filtered_RNU42[, c("Part1", "Part2")]

LRep2_Filtered_RNU42 <- LRep2_Filtered_RNU42[, c("Part1", "Part2")]

RNAInter_RNU42 <- RNAInter_RNU42[,c("Interactor1.Symbol", "Interactor2.Symbol")]
colnames(RNAInter_RNU42) <- c("Part1", "Part2")

NPinter_RNU42 <- NPinter_RNU42[, c("Molecule1", "Molecule2")]
colnames(NPinter_RNU42) <- c("Part1", "Part2")

#All RNU4-2 interactions, All filtered RNU4-2 combined
All_RNU42_interactions <- rbind(LRep1_Filtered_RNU42, LRep2_Filtered_RNU42, NPinter_RNU42, RNAInter_RNU42)


RNU42_interaction_counts <- All_RNU42_interactions %>%
  group_by(Part1, Part2) %>%
  summarize(count = n(), .groups = "drop")

Counts <- rbind(RNU41_interaction_counts, RNU42_interaction_counts)


#Fix the encode values
Counts[85, 2] <- "lnc-AC007952.1.1-1"
Counts[82, 2] <- "lnc-LRR1-1"
Counts[84, 2] <- "lnc-AC007952.2-2"
Counts[83, 2] <- "lnc-GRAP-1"

#Remove the duplicate values
Counts[68, 3] <- 2
Counts <- Counts[-73, ]

Counts[92, 3] <- 2
Counts <- Counts[-73, ]

Counts[84, 3] <- 2
Counts <- Counts[-76, ]

#Remove RNU4-1 interactions with RNU$-2 interactions
Counts <- Counts[-106, ]


row.names(Counts) <- NULL

#Rename FAM166C to CIMIP2C
Counts[83, 2] <- "CIMIP2C"


# Create the RNA vector with no extra commas
RNA_vector <- c(
  "rRNA", "lincRNA", "mRNA", "lincRNA", "mRNA", "lincRNA", "lincRNA", "srpRNA", "srpRNA", "srpRNA", "rRNA",
  "rRNA", "rRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA",
  "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA",
  "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA",
  "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA",
  "snRNA", "snRNA", "snRNA", "mRNA", "snoRNA", "snoRNA", "snoRNA", "mRNA", "mRNA", "snRNA", "mRNA", "lincRNA",
  "mRNA", "mRNA", "mRNA", "mRNA", "mRNA", "lincRNA", "lincRNA", "lincRNA", "lincRNA", "mRNA", "lincRNA", "mRNA",
  "mRNA", "mRNA", "mRNA", "mRNA", "mRNA", "lincRNA", "mRNA", "srpRNA", "srpRNA", "srpRNA", "srpRNA", "srpRNA",
  "srpRNA", "srpRNA", "rRNA", "rRNA", "rRNA", "rRNA", "rRNA", "rRNA", "snRNA", "snRNA", "snRNA", "snRNA",
  "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA",
  "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA",
  "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA",
  "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA",
  "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA",
  "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", 
  "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNa", "snRNA", "snRNA", "snRNA",
  "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "rRNA", "mRNA", "snoRNA", "snoRNA", "snoRNA", "snoRNA", "snoRNA",
  "snoRNA", "snoRNA", "snoRNA", "mRNA", "mRNA", "mRNA", "mRNA", "snRNA", "mRNA", "mRNA" 
)
  

# Add the RNA_vector as a new column to the All_RNU4 dataframe
Counts$RNA <- RNA_vector


write.csv(Counts, "/home/jessicaaucott/Desktop/Counts.csv", row.names = FALSE)

  