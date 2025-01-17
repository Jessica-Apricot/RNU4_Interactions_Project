install.packages("dplyr")
library(dplyr)

install.packages("tidyr")
library(tidyr)

# Specify the file path
LIGR_seq_Rep1_file_path <- "large_data/GSM2113739_S1_AMT_Ligase_Rep1_v19-48.rmsk.uniq.unfilt.lig.txt"

# Read the tab-delimited file
LIGR_seq_Rep1 <- read.delim(LIGR_seq_Rep1_file_path)

#Set column names
colnames(LIGR_seq_Rep1) <- c("Interaction_status", 
                             "Chimeric_read_structure", 
                             "Local_alignments_of_chimera", 
                             "Gene_Symbols", "Gene_IDs", 
                             "Transcript_IDs", 
                             "RNA_catagory", 
                             "RepeatMasker_annotation", 
                             "RepeatMasker_annotation_catagory", 
                             "Read_name", "Read_Sequence", 
                             "Alignment_Score", "Deprecated_Chimera_uniqueness",
                             "Alignment_positions_in_transcripts", 
                             "RepeatMask_offset_positions", "Alignment_lengths", 
                             "Genome_start_position_of_alignment_segment", "1",
                             "2", "3", "4", "5", "6", "7", "8")


#Filter out only rows with RNU4
# Filter only rows with exact RNU4-1 in Gene_Symbols
LRep1_Filtered_RNU41 <- LIGR_seq_Rep1[grepl("\\bRNU4-1\\b", LIGR_seq_Rep1$Gene_Symbols),]

# Filter only rows with exact RNU4-2 in Gene_Symbols
LRep1_Filtered_RNU42 <- LIGR_seq_Rep1[grepl("\\bRNU4-2\\b", LIGR_seq_Rep1$Gene_Symbols),]



#only RNU4 that is interacting with molecules other than self
LRep1_Filtered_RNU41 <- LRep1_Filtered_RNU41[!grepl("S", LRep1_Filtered_RNU41$Interaction_status),]

LRep1_Filtered_RNU42 <- LRep1_Filtered_RNU42[!grepl("S", LRep1_Filtered_RNU42$Interaction_status),]

#Separate the interactions into two columns
LRep1_Filtered_RNU41 <- LRep1_Filtered_RNU41 %>%
  separate(Gene_Symbols, into = c("Part1", "Part2"), sep = ":")

LRep1_Filtered_RNU42 <- LRep1_Filtered_RNU42 %>%
  separate(Gene_Symbols, into = c("Part1", "Part2"), sep = ":")

#Put RNU4 interactions in the same column
LRep1_Filtered_RNU41 <- LRep1_Filtered_RNU41 %>%
  mutate(
    Part2 = ifelse(Part2 == "RNU4-1", Part1, Part2),
    Part1 = ifelse(Part1 != "RNU4-1", "RNU4-1", Part1)
  )

LRep1_Filtered_RNU42 <- LRep1_Filtered_RNU42 %>%
  mutate(
    Part2 = ifelse(Part2 == "RNU4-2", Part1, Part2),
    Part1 = ifelse(Part1 != "RNU4-2", "RNU4-2", Part1)
  )


#Whats exlusive to RNU4-2?
LRep1_In_RNU42_only <- anti_join(LRep1_Filtered_RNU42, LRep1_Filtered_RNU41, by='Part2')




LIGR_seq_Rep2_file_path <- "large_data/GSM2113743_S5_AMT_Ligase_Rep2_v19-48.rmsk.uniq.unfilt.lig.txt"

# Read the tab-delimited file
LIGR_seq_Rep2 <- read.delim(LIGR_seq_Rep2_file_path)

#Set column names
colnames(LIGR_seq_Rep2) <- c("Interaction_status", "Chimeric_read_structure",
                             "Local_alignments_of_chimera", "Gene_Symbols",
                             "Gene_IDs", "Transcript_IDs", "RNA_catagory",
                             "RepeatMasker_annotation", 
                             "RepeatMasker_annotation catagory", "Read_name",
                             "Read_Sequence", "Alignment_Score", 
                             "Deprecated_Chimera_uniqueness", 
                             "Alignment_positions_in_transcripts", 
                             "RepeatMask_offset_positions", "Alignment_lengths",
                             "Genome_start_position_of_alignment_segment", "1",
                             "2", "3", "4", "5", "6", "7","8")

#Only RNU4-1 and -2 interactions
LRep2_Filtered_RNU41 <- LIGR_seq_Rep2[grepl("\\bRNU4-1\\b",LIGR_seq_Rep2$Gene_Symbols),]

LRep2_Filtered_RNU42 <- LIGR_seq_Rep2[grepl("\\bRNU4-2\\b",LIGR_seq_Rep2$Gene_Symbols),]

#only RNU4 that is interacting with molecules other than self
LRep2_Filtered_RNU41 <- LRep2_Filtered_RNU41[!grepl("S", LRep2_Filtered_RNU41$Interaction_status),]

LRep2_Filtered_RNU42 <- LRep2_Filtered_RNU42[!grepl("S", LRep2_Filtered_RNU42$Interaction_status),]

#Separate the interactions into two columns
LRep2_Filtered_RNU41 <- LRep2_Filtered_RNU41 %>%
  separate(Gene_Symbols, into = c("Part1", "Part2"), sep = ":")

LRep2_Filtered_RNU42 <- LRep2_Filtered_RNU42 %>%
  separate(Gene_Symbols, into = c("Part1", "Part2"), sep = ":")

#Put RNU4 interactions in the same column
LRep2_Filtered_RNU41 <- LRep2_Filtered_RNU41 %>%
  mutate(
    Part2 = ifelse(Part2 == "RNU4-1", Part1, Part2),
    Part1 = ifelse(Part1 != "RNU4-1", "RNU4-1", Part1)
  )

LRep2_Filtered_RNU42 <- LRep2_Filtered_RNU42 %>%
  mutate(
    Part2 = ifelse(Part2 == "RNU4-2", Part1, Part2),
    Part1 = ifelse(Part1 != "RNU4-2", "RNU4-2", Part1)
  )

#Remove any left over U4 interactions
LRep2_Filtered_RNU41 <- LRep2_Filtered_RNU41[!grepl("U4",LRep2_Filtered_RNU41$Part2),]
LRep2_Filtered_RNU42 <- LRep2_Filtered_RNU42[!grepl("U4",LRep2_Filtered_RNU42$Part2),]

#Whats exclusive to RNU4-2?
LRep2_In_RNU42_only <- anti_join(LRep2_Filtered_RNU42, LRep2_Filtered_RNU41, by='Part2') 




   
RNAInter_file_path <- "large_data/RNAInter_full_interactions_data.txt"

# Read the tab-delimited file
RNAInter <- read.delim(RNAInter_file_path)

#Only human
RNAInter <- RNAInter[grepl("Homo sapiens",RNAInter$Species1),]

#RNU4-1 in any column
RNAInter_RNU41 <- RNAInter %>%
  filter(apply(., 1, function(row) any(grepl("\\bRNU4-1\\b", row))))

RNAInter_RNU42 <- RNAInter %>%
  filter(apply(., 1, function(row) any(grepl("\\bRNU4-2\\b", row))))


#RNU4-2 in the same row
RNAInter_RNU42 <- RNAInter_RNU42 %>%
  mutate(
    Interactor2.Symbol = ifelse(Interactor2.Symbol == "RNU4-2", Interactor1.Symbol, Interactor2.Symbol),
    Interactor1.Symbol = ifelse(Interactor1.Symbol != "RNU4-2", "RNU4-2", Interactor1.Symbol)
  )

#RNU4-1 in the same row
RNAInter_RNU41 <- RNAInter_RNU41 %>%
  mutate(
    Interactor2.Symbol = ifelse(Interactor2.Symbol == "RNU4-1", Interactor1.Symbol, Interactor2.Symbol),
    Interactor1.Symbol = ifelse(Interactor1.Symbol != "RNU4-1", "RNU4-1", Interactor1.Symbol) 
  )

#Whats exlusive to RNU4-2?
RNAInter_RNU42_only <- anti_join(RNAInter_RNU42, RNAInter_RNU41, by='Interactor2.Symbol') 


#NPinter interactions

NPinter_file_path <- "large_data/snRNA_interaction.txt"

# Read the tab-delimited file
NPinter <- read.delim(NPinter_file_path)

#Column Names
colnames(NPinter) <- c("M1_Interaction_ID", "Molecule1", 
                       "Gene_identifier", "RNA_catagory",
                       "Molecule2", "M2_Interaction_ID", 
                       "Molecule_type", "Data_source", 
                       "Assay_type", "database_ID", "Species",
                       "Cell_type", "Interaction_catagory", 
                       "Interaction_class", "Interaction_catagory_2", "Source")

#Human only
NPinter <- NPinter[grepl("Homo sapiens",NPinter$Species),]

#RNU4-1 in any column
NPinter_RNU41 <- NPinter %>%
  filter(apply(., 1, function(row) any(grepl("\\bRNU4-1\\b",row))))


#RNU4-2 in any column

NPinter_RNU42 <- NPinter %>%
  filter(apply(., 1, function(row) any(grepl("\\bRNU4-1\\b",row))))

#RNU4-2 in the same row
NPinter_RNU42 <- NPinter_RNU42 %>%
  mutate(
    Molecule2 = ifelse(Molecule2 == "RNU4-2", Molecule1, Molecule2),
    Molecule1 = ifelse(Molecule1 != "RNU4-2", "RNU4-2", Molecule1)
  )

#RNU4-1 in the same row
NPinter_RNU41 <- NPinter_RNU41 %>%
  mutate(
    Molecule2 = ifelse(Molecule2 == "RNU4-1", Molecule1, Molecule2),
    Molecule1 = ifelse(Molecule1 != "RNU4-1", "RNU4-1", Molecule1)
  )

#Whats exclusive to RNU4-2?
NPinter_RNU42_only <- anti_join(NPinter_RNU42, NPinter_RNU41, by='Molecule2') 





#Fixing up and separating RNU4-1 interactions

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
U4_RNA_vector <- c(
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
Counts$RNA <- U4_RNA_vector


