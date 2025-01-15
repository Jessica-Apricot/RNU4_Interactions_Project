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


#Filter out only rows with RNU6
# Filter only rows with exact RNU4-1 in Gene_Symbols
LRep1_Filtered_RNU6_1 <- LIGR_seq_Rep1[grepl("\\bRNU6-1\\b", LIGR_seq_Rep1$Gene_Symbols),]

# Filter only rows with exact RNU6-2 in Gene_Symbols
LRep1_Filtered_RNU6_2 <- LIGR_seq_Rep1[grepl("\\bRNU6-2\\b", LIGR_seq_Rep1$Gene_Symbols),]

# Filter only rows with exact RNU6-7 in Gene_Symbols
LRep1_Filtered_RNU6_7 <- LIGR_seq_Rep1[grepl("\\bRNU6-7\\b", LIGR_seq_Rep1$Gene_Symbols),]

# Filter only rows with exact RNU6-8 in Gene_Symbols
LRep1_Filtered_RNU6_8 <- LIGR_seq_Rep1[grepl("\\bRNU6-8\\b", LIGR_seq_Rep1$Gene_Symbols),]

# Filter only rows with exact RNU6-9 in Gene_Symbols
LRep1_Filtered_RNU6_9 <- LIGR_seq_Rep1[grepl("\\bRNU6-9\\b", LIGR_seq_Rep1$Gene_Symbols),]



#only RNU6 that is interacting with molecules other than self
LRep1_Filtered_RNU6_1 <- LRep1_Filtered_RNU6_1[!grepl("S", LRep1_Filtered_RNU6_1$Interaction_status),]
LRep1_Filtered_RNU6_2 <- LRep1_Filtered_RNU6_2[!grepl("S", LRep1_Filtered_RNU6_2$Interaction_status),]
LRep1_Filtered_RNU6_7 <- LRep1_Filtered_RNU6_7[!grepl("S", LRep1_Filtered_RNU6_7$Interaction_status),]
LRep1_Filtered_RNU6_8 <- LRep1_Filtered_RNU6_8[!grepl("S", LRep1_Filtered_RNU6_8$Interaction_status),]
LRep1_Filtered_RNU6_9 <- LRep1_Filtered_RNU6_9[!grepl("S", LRep1_Filtered_RNU6_9$Interaction_status),]

#Separate the interactions into two columns
LRep1_Filtered_RNU6_1 <- LRep1_Filtered_RNU6_1 %>%
  separate(Gene_Symbols, into = c("Part1", "Part2"), sep = ":")

LRep1_Filtered_RNU6_2 <- LRep1_Filtered_RNU6_2 %>%
  separate(Gene_Symbols, into = c("Part1", "Part2"), sep = ":")

LRep1_Filtered_RNU6_7 <- LRep1_Filtered_RNU6_7 %>%
  separate(Gene_Symbols, into = c("Part1", "Part2"), sep = ":")

LRep1_Filtered_RNU6_8 <- LRep1_Filtered_RNU6_8 %>%
  separate(Gene_Symbols, into = c("Part1", "Part2"), sep = ":")

LRep1_Filtered_RNU6_9 <- LRep1_Filtered_RNU6_9 %>%
  separate(Gene_Symbols, into = c("Part1", "Part2"), sep = ":")


#Put RNU4 interactions in the same column
LRep1_Filtered_RNU6_1 <- LRep1_Filtered_RNU6_1 %>%
  mutate(
    Part2 = ifelse(Part2 == "RNU6-1", Part1, Part2),
    Part1 = ifelse(Part1 != "RNU6-1", "RNU6-1", Part1)
  )

LRep1_Filtered_RNU6_2 <- LRep1_Filtered_RNU6_2 %>%
  mutate(
    Part2 = ifelse(Part2 == "RNU6-2", Part1, Part2),
    Part1 = ifelse(Part1 != "RNU6-2", "RNU6-2", Part1)
  )

LRep1_Filtered_RNU6_7 <- LRep1_Filtered_RNU6_7 %>%
  mutate(
    Part2 = ifelse(Part2 == "RNU6-7", Part1, Part2),
    Part1 = ifelse(Part1 != "RNU6-7", "RNU6-7", Part1)
  )

LRep1_Filtered_RNU6_8 <- LRep1_Filtered_RNU6_8 %>%
  mutate(
    Part2 = ifelse(Part2 == "RNU6-8", Part1, Part2),
    Part1 = ifelse(Part1 != "RNU6-8", "RNU6-8", Part1)
  )

LRep1_Filtered_RNU6_9 <- LRep1_Filtered_RNU6_9 %>%
  mutate(
    Part2 = ifelse(Part2 == "RNU6-9", Part1, Part2),
    Part1 = ifelse(Part1 != "RNU6-9", "RNU6-9", Part1)
  )

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

#Filter out only rows with RNU6
# Filter only rows with exact RNU4-1 in Gene_Symbols
LRep2_Filtered_RNU6_1 <- LIGR_seq_Rep2[grepl("\\bRNU6-1\\b", LIGR_seq_Rep2$Gene_Symbols),]

# Filter only rows with exact RNU6-2 in Gene_Symbols
LRep2_Filtered_RNU6_2 <- LIGR_seq_Rep2[grepl("\\bRNU6-2\\b", LIGR_seq_Rep2$Gene_Symbols),]

# Filter only rows with exact RNU6-7 in Gene_Symbols
LRep2_Filtered_RNU6_7 <- LIGR_seq_Rep2[grepl("\\bRNU6-7\\b", LIGR_seq_Rep2$Gene_Symbols),]

# Filter only rows with exact RNU6-8 in Gene_Symbols
LRep2_Filtered_RNU6_8 <- LIGR_seq_Rep2[grepl("\\bRNU6-8\\b", LIGR_seq_Rep2$Gene_Symbols),]

# Filter only rows with exact RNU6-9 in Gene_Symbols
LRep2_Filtered_RNU6_9 <- LIGR_seq_Rep2[grepl("\\bRNU6-9\\b", LIGR_seq_Rep2$Gene_Symbols),]



#only RNU4 that is interacting with molecules other than self
LRep2_Filtered_RNU6_1 <- LRep2_Filtered_RNU6_1[!grepl("S", LRep2_Filtered_RNU6_1$Interaction_status),]
LRep2_Filtered_RNU6_2 <- LRep2_Filtered_RNU6_2[!grepl("S", LRep2_Filtered_RNU6_2$Interaction_status),]
LRep2_Filtered_RNU6_7 <- LRep2_Filtered_RNU6_7[!grepl("S", LRep2_Filtered_RNU6_7$Interaction_status),]
LRep2_Filtered_RNU6_8 <- LRep2_Filtered_RNU6_8[!grepl("S", LRep2_Filtered_RNU6_8$Interaction_status),]
LRep2_Filtered_RNU6_9 <- LRep2_Filtered_RNU6_9[!grepl("S", LRep2_Filtered_RNU6_9$Interaction_status),]


#Separate the interactions into two columns
LRep2_Filtered_RNU6_1 <- LRep2_Filtered_RNU6_1 %>%
  separate(Gene_Symbols, into = c("Part1", "Part2"), sep = ":")
LRep2_Filtered_RNU6_2 <- LRep2_Filtered_RNU6_2 %>%
  separate(Gene_Symbols, into = c("Part1", "Part2"), sep = ":")
LRep2_Filtered_RNU6_7 <- LRep2_Filtered_RNU6_7 %>%
  separate(Gene_Symbols, into = c("Part1", "Part2"), sep = ":")
LRep2_Filtered_RNU6_8 <- LRep2_Filtered_RNU6_8 %>%
  separate(Gene_Symbols, into = c("Part1", "Part2"), sep = ":")
LRep2_Filtered_RNU6_9 <- LRep2_Filtered_RNU6_9 %>%
  separate(Gene_Symbols, into = c("Part1", "Part2"), sep = ":")


#Put RNU4 interactions in the same column
LRep2_Filtered_RNU6_1 <- LRep2_Filtered_RNU6_1 %>%
  mutate(
    Part2 = ifelse(Part2 == "RNU6-1", Part1, Part2),
    Part1 = ifelse(Part1 != "RNU6-1", "RNU6-1", Part1)
  )
LRep2_Filtered_RNU6_2 <- LRep2_Filtered_RNU6_2 %>%
  mutate(
    Part2 = ifelse(Part2 == "RNU6-2", Part1, Part2),
    Part1 = ifelse(Part1 != "RNU6-2", "RNU6-2", Part1)
  )
LRep2_Filtered_RNU6_7 <- LRep2_Filtered_RNU6_7 %>%
  mutate(
    Part2 = ifelse(Part2 == "RNU6-7", Part1, Part2),
    Part1 = ifelse(Part1 != "RNU6-7", "RNU6-7", Part1)
  )
LRep2_Filtered_RNU6_8 <- LRep2_Filtered_RNU6_8 %>%
  mutate(
    Part2 = ifelse(Part2 == "RNU6-8", Part1, Part2),
    Part1 = ifelse(Part1 != "RNU6-8", "RNU6-8", Part1)
  )
LRep2_Filtered_RNU6_9 <- LRep2_Filtered_RNU6_9 %>%
  mutate(
    Part2 = ifelse(Part2 == "RNU6-9", Part1, Part2),
    Part1 = ifelse(Part1 != "RNU6-9", "RNU6-9", Part1)
  )



RNAInter_file_path <- "large_data/RNAInter_full_interactions_data.txt"

# Read the tab-delimited file
RNAInter <- read.delim(RNAInter_file_path)

#Only human
RNAInter <- RNAInter[grepl("Homo sapiens",RNAInter$Species1),]

#Filter out only rows with RNU6
RNAInter_RNU6_1 <- RNAInter %>%
  filter(apply(., 1, function(row) any(grepl("\\bRNU6-1\\b", row))))

RNAInter_RNU6_2 <- RNAInter %>%
  filter(apply(., 1, function(row) any(grepl("\\bRNU6-2\\b", row))))

RNAInter_RNU6_7 <- RNAInter %>%
  filter(apply(., 1, function(row) any(grepl("\\bRNU6-7\\b", row))))

RNAInter_RNU6_8 <- RNAInter %>%
  filter(apply(., 1, function(row) any(grepl("\\bRNU6-8\\b", row))))

RNAInter_RNU6_9 <- RNAInter %>%
  filter(apply(., 1, function(row) any(grepl("\\bRNU6-9\\b", row))))



#RNU6 in the same row
RNAInter_RNU6_1 <- RNAInter_RNU6_1 %>%
  mutate(
    Interactor2.Symbol = ifelse(Interactor2.Symbol == "RNU6-1", Interactor1.Symbol, Interactor2.Symbol),
    Interactor1.Symbol = ifelse(Interactor1.Symbol != "RNU6-1", "RNU6-1", Interactor1.Symbol)
  )

RNAInter_RNU6_2 <- RNAInter_RNU6_2 %>%
  mutate(
    Interactor2.Symbol = ifelse(Interactor2.Symbol == "RNU6-2", Interactor1.Symbol, Interactor2.Symbol),
    Interactor1.Symbol = ifelse(Interactor1.Symbol != "RNU6-2", "RNU6-2", Interactor1.Symbol)
  )

RNAInter_RNU6_7 <- RNAInter_RNU6_7 %>%
  mutate(
    Interactor2.Symbol = ifelse(Interactor2.Symbol == "RNU6-7", Interactor1.Symbol, Interactor2.Symbol),
    Interactor1.Symbol = ifelse(Interactor1.Symbol != "RNU6-7", "RNU6-7", Interactor1.Symbol)
  )

RNAInter_RNU6_8 <- RNAInter_RNU6_8 %>%
  mutate(
    Interactor2.Symbol = ifelse(Interactor2.Symbol == "RNU6-8", Interactor1.Symbol, Interactor2.Symbol),
    Interactor1.Symbol = ifelse(Interactor1.Symbol != "RNU6-8", "RNU6-8", Interactor1.Symbol)
  )

RNAInter_RNU6_9 <- RNAInter_RNU6_9 %>%
  mutate(
    Interactor2.Symbol = ifelse(Interactor2.Symbol == "RNU6-9", Interactor1.Symbol, Interactor2.Symbol),
    Interactor1.Symbol = ifelse(Interactor1.Symbol != "RNU6-9", "RNU6-9", Interactor1.Symbol)
  )



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

#RNU6 in any column
NPinter_RNU6_1 <- NPinter %>%
  filter(apply(., 1, function(row) any(grepl("\\bRNU6-1\\b",row))))

NPinter_RNU6_2 <- NPinter %>%
  filter(apply(., 1, function(row) any(grepl("\\bRNU6-2\\b",row))))

NPinter_RNU6_7 <- NPinter %>%
  filter(apply(., 1, function(row) any(grepl("\\bRNU6-7\\b",row))))

NPinter_RNU6_8 <- NPinter %>%
  filter(apply(., 1, function(row) any(grepl("\\bRNU6-8\\b",row))))

NPinter_RNU6_9 <- NPinter %>%
  filter(apply(., 1, function(row) any(grepl("\\bRNU6-9\\b",row))))

#RNU6 in the same row
NPinter_RNU6_1 <- NPinter_RNU6_1 %>%
  mutate(
    Molecule2 = ifelse(Molecule2 == "RNU6-1", Molecule1, Molecule2),
    Molecule1 = ifelse(Molecule1 != "RNU6-1", "RNU6-1", Molecule1)
  )
NPinter_RNU6_2 <- NPinter_RNU6_2 %>%
  mutate(
    Molecule2 = ifelse(Molecule2 == "RNU6-2", Molecule1, Molecule2),
    Molecule1 = ifelse(Molecule1 != "RNU6-2", "RNU6-2", Molecule1)
  )
NPinter_RNU6_7 <- NPinter_RNU6_7 %>%
  mutate(
    Molecule2 = ifelse(Molecule2 == "RNU6-7", Molecule1, Molecule2),
    Molecule1 = ifelse(Molecule1 != "RNU6-7", "RNU6-7", Molecule1)
  )
NPinter_RNU6_8 <- NPinter_RNU6_8 %>%
  mutate(
    Molecule2 = ifelse(Molecule2 == "RNU6-8", Molecule1, Molecule2),
    Molecule1 = ifelse(Molecule1 != "RNU6-8", "RNU6-8", Molecule1)
  )
NPinter_RNU6_9 <- NPinter_RNU6_9 %>%
  mutate(
    Molecule2 = ifelse(Molecule2 == "RNU6-9", Molecule1, Molecule2),
    Molecule1 = ifelse(Molecule1 != "RNU6-9", "RNU6-9", Molecule1)
  )

#Fixing up and separating RNU4-1 interactions

#Selecting for only the interaction columns and renaming them all P1 or P2
LRep1_Filtered_RNU6_1 <- LRep1_Filtered_RNU6_1[, c("Part1", "Part2")]
LRep1_Filtered_RNU6_2 <- LRep1_Filtered_RNU6_2[, c("Part1", "Part2")]
LRep1_Filtered_RNU6_7 <- LRep1_Filtered_RNU6_7[, c("Part1", "Part2")]
LRep1_Filtered_RNU6_8 <- LRep1_Filtered_RNU6_8[, c("Part1", "Part2")]
LRep1_Filtered_RNU6_9 <- LRep1_Filtered_RNU6_9[, c("Part1", "Part2")]

LRep2_Filtered_RNU6_1 <- LRep2_Filtered_RNU6_1[, c("Part1", "Part2")]
LRep2_Filtered_RNU6_2 <- LRep2_Filtered_RNU6_2[, c("Part1", "Part2")]
LRep2_Filtered_RNU6_7 <- LRep2_Filtered_RNU6_7[, c("Part1", "Part2")]
LRep2_Filtered_RNU6_8 <- LRep2_Filtered_RNU6_8[, c("Part1", "Part2")]
LRep2_Filtered_RNU6_9 <- LRep2_Filtered_RNU6_9[, c("Part1", "Part2")]

RNAInter_RNU6_1 <- RNAInter_RNU6_1[, c("Interactor1.Symbol", "Interactor2.Symbol")]
colnames(RNAInter_RNU6_1) <- c("Part1", "Part2")
RNAInter_RNU6_2 <- RNAInter_RNU6_2[, c("Interactor1.Symbol", "Interactor2.Symbol")]
colnames(RNAInter_RNU6_2) <- c("Part1", "Part2")
RNAInter_RNU6_7 <- RNAInter_RNU6_7[, c("Interactor1.Symbol", "Interactor2.Symbol")]
colnames(RNAInter_RNU6_7) <- c("Part1", "Part2")
RNAInter_RNU6_8 <- RNAInter_RNU6_8[, c("Interactor1.Symbol", "Interactor2.Symbol")]
colnames(RNAInter_RNU6_8) <- c("Part1", "Part2")
RNAInter_RNU6_9 <- RNAInter_RNU6_9[, c("Interactor1.Symbol", "Interactor2.Symbol")]
colnames(RNAInter_RNU6_9) <- c("Part1", "Part2")

NPinter_RNU6_1 <- NPinter_RNU6_1[, c("Molecule1", "Molecule2")]
colnames(NPinter_RNU6_1) <- c("Part1", "Part2")
NPinter_RNU6_2 <- NPinter_RNU6_2[, c("Molecule1", "Molecule2")]
colnames(NPinter_RNU6_2) <- c("Part1", "Part2")
NPinter_RNU6_7 <- NPinter_RNU6_7[, c("Molecule1", "Molecule2")]
colnames(NPinter_RNU6_7) <- c("Part1", "Part2")
NPinter_RNU6_8 <- NPinter_RNU6_8[, c("Molecule1", "Molecule2")]
colnames(NPinter_RNU6_8) <- c("Part1", "Part2")
NPinter_RNU6_9 <- NPinter_RNU6_9[, c("Molecule1", "Molecule2")]
colnames(NPinter_RNU6_9) <- c("Part1", "Part2")


All_RNU6_1_interactions <- rbind(LRep1_Filtered_RNU6_1, LRep2_Filtered_RNU6_1, NPinter_RNU6_1, RNAInter_RNU6_1)
All_RNU6_2_interactions <- rbind(LRep1_Filtered_RNU6_2, LRep2_Filtered_RNU6_2, NPinter_RNU6_2, RNAInter_RNU6_2)
All_RNU6_7_interactions <- rbind(LRep1_Filtered_RNU6_7, LRep2_Filtered_RNU6_7, NPinter_RNU6_7, RNAInter_RNU6_7)
All_RNU6_8_interactions <- rbind(LRep1_Filtered_RNU6_8, LRep2_Filtered_RNU6_8, NPinter_RNU6_8, RNAInter_RNU6_8)
All_RNU6_9_interactions <- rbind(LRep1_Filtered_RNU6_9, LRep2_Filtered_RNU6_9, NPinter_RNU6_9, RNAInter_RNU6_9)

U6_Counts <- rbind(All_RNU6_1_interactions, All_RNU6_2_interactions, All_RNU6_7_interactions, All_RNU6_8_interactions, All_RNU6_9_interactions)

U6_Counts <- U6_Counts %>%
  group_by(Part1, Part2) %>%
  summarize(count = n(), .groups = "drop")

U6_Counts[9,2] <- "lnc-AC087650.1-1"

# Create the RNA vector with no extra commas
U6_RNA_vector <- c(
  "lincRNA", "rRNA", "lincRNA", "lincRNA", "mRNA", "mRNA", "mRNA", "mRNA", "lincRNA", "mRNA", "mRNA", "mRNA",
  "lincRNA", "lincRNA", "ncRNA", "lincRNA", "srpRNA", "rRNA", "rRNA", "rRNA", "snRNA", "snRNA", "snRNA",
  "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", 
  "snRNA", "snRNA", "rRNA", "snoRNA", "snoRNA", "snoRNA", "snoRNA", "snoRNA", "snoRNA", "snoRNA", "snRNA", 
  "mRNA", "mRNA", "ncRNA", "lincRNA")

U6_Counts$RNA <- U6_RNA_vector

write.csv(U6_Counts, "U6_Counts.csv", row.names = TRUE)
