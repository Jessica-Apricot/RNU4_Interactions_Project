install.packages("dplyr")
library(dplyr)

install.packages("tidyr")
library(tidyr)

# Specify the file path
LIGR_seq_Rep1_file_path <- "large_data/GSM2113739_S1_AMT_Ligase_Rep1_v19-48.rmsk.uniq.unfilt.lig.txt"

# Read the tab-delimited file
LIGR_seq_Rep1 <- read.delim(LIGR_seq_file_path)

#Set column names
colnames(LIGR_seq_Rep1) <- c("Interaction status", 
                             "Chimeric read structure", 
                             "Local alignments of chimera", 
                             "Gene Symbols", "Gene IDs", 
                             "Transcript IDs", 
                             "RNA catagory", 
                             "RepeatMasker annotation", 
                             "RepeatMasker annotation catagory", 
                             "Read name", "Read Sequence", 
                             "Alignment Score", "Deprecated Chimera uniqueness",
                             "Alignment positions in transcripts", 
                             "RepeatMask offset positions", "Alignment lengths", 
                             "Genome start position of alignment segment", "1",
                             "2", "3", "4", "5", "6", "7", "8")


#Filter out only rows with RNU4
LRep1_Filtered_RNU41 <- LIGR_seq_Rep1[grepl("RNU4-1",LIGR_seq_Rep1$`Gene Symbols`),]


LRep1_Filtered_RNU42 <- LIGR_seq_Rep1[grepl("RNU4-2",LIGR_seq_Rep1$`Gene Symbols`),]

#only RNU4 that is interacting with molecules other than self
LRep1_Only_Intermolecular_RNU41 <- LRep1_Filtered_RNU41[!grepl("S", LRep1_Filtered_RNU41$`Interaction status`),]

LRep1_Only_Intermolecular_RNU42 <- LRep1_Filtered_RNU42[!grepl("S", LRep1_Filtered_RNU42$`Interaction status`),]

#Separate the interactions into two columns
LRep1_Only_Intermolecular_RNU41 <- LRep1_Only_Intermolecular_RNU41 %>%
  separate(`Gene Symbols`, into = c("Part1", "Part2"), sep = ":")

LRep1_Only_Intermolecular_RNU42 <- LRep1_Only_Intermolecular_RNU42 %>%
  separate(`Gene Symbols`, into = c("Part1", "Part2"), sep = ":")

#Put RNU4 interactions in the same column
LRep1_Only_Intermolecular_RNU41 <- LRep1_Only_Intermolecular_RNU41 %>%
  mutate(
    Part2 = ifelse(Part2 == "RNU4-1", Part1, Part2),
    Part1 = ifelse(Part1 != "RNU4-1", "RNU4-1", Part1)
  )

LRep1_Only_Intermolecular_RNU42 <- LRep1_Only_Intermolecular_RNU42 %>%
  mutate(
    Part2 = ifelse(Part2 == "RNU4-2", Part1, Part2),
    Part1 = ifelse(Part1 != "RNU4-2", "RNU4-2", Part1)
  )

#Whats exlusive to RNU4-2?
LRep1_In_RNU42_only <- anti_join(LRep1_Only_Intermolecular_RNU42, LRep1_Only_Intermolecular_RNU41, by='Part2')

#Fixing lncRNA values = ENSG00000258486.2
LRep1_Only_Intermolecular_RNU42[15, "Part2"] <- "lnc-LRR1-1"

#Remove U6 interactions
LRep1_In_RNU42_only <- LRep1_In_RNU42_only[!grepl("U6",LRep1_In_RNU42_only$Part2),]



LIGR_seq_Rep2_file_path <- "large_data/GSM2113743_S5_AMT_Ligase_Rep2_v19-48.rmsk.uniq.unfilt.lig.txt"

# Read the tab-delimited file
LIGR_seq_Rep2 <- read.delim(LIGR_seq_Rep2_file_path)

#Set column names
colnames(LIGR_seq_Rep2) <- c("Interaction status", "Chimeric read structure",
                             "Local alignments of chimera", "Gene Symbols",
                             "Gene IDs", "Transcript IDs", "RNA catagory",
                             "RepeatMasker annotation", 
                             "RepeatMasker annotation catagory", "Read name",
                             "Read Sequence", "Alignment Score", 
                             "Deprecated Chimera uniqueness", 
                             "Alignment positions in transcripts", 
                             "RepeatMask offset positions", "Alignment lengths",
                             "Genome start position of alignment segment", "1",
                             "2", "3", "4", "5", "6", "7","8")

#Only RNU4-1 and -2 interactions
LRep2_Filtered_RNU41 <- LIGR_seq_Rep2[grepl("RNU4-1",LIGR_seq_Rep2$`Gene Symbols`),]

LRep2_Filtered_RNU42 <- LIGR_seq_Rep2[grepl("RNU4-2",LIGR_seq_Rep2$`Gene Symbols`),]

#only RNU4 that is interacting with molecules other than self
LRep2_Only_Intermolecular_RNU41 <- LRep2_Filtered_RNU41[!grepl("S", LRep2_Filtered_RNU41$`Interaction status`),]

LRep2_Only_Intermolecular_RNU42 <- LRep2_Filtered_RNU42[!grepl("S", LRep2_Filtered_RNU42$`Interaction status`),]

#Separate the interactions into two columns
LRep2_Only_Intermolecular_RNU41 <- LRep2_Only_Intermolecular_RNU41 %>%
  separate(`Gene Symbols`, into = c("Part1", "Part2"), sep = ":")

LRep2_Only_Intermolecular_RNU42 <- LRep2_Only_Intermolecular_RNU42 %>%
  separate(`Gene Symbols`, into = c("Part1", "Part2"), sep = ":")

#Put RNU4 interactions in the same column
LRep2_Only_Intermolecular_RNU41 <- LRep2_Only_Intermolecular_RNU41 %>%
  mutate(
    Part2 = ifelse(Part2 == "RNU4-1", Part1, Part2),
    Part1 = ifelse(Part1 != "RNU4-1", "RNU4-1", Part1)
  )

LRep2_Only_Intermolecular_RNU42 <- LRep2_Only_Intermolecular_RNU42 %>%
  mutate(
    Part2 = ifelse(Part2 == "RNU4-2", Part1, Part2),
    Part1 = ifelse(Part1 != "RNU4-2", "RNU4-2", Part1)
  )

#Whats exlusive to RNU4-2?
LRep2_In_RNU42_only <- anti_join(LRep2_Only_Intermolecular_RNU42, LRep2_Only_Intermolecular_RNU41, by='Part2') 

#Remove U6 interactions
LRep2_In_RNU42_only <- LRep2_In_RNU42_only[!grepl("U6",LRep2_In_RNU42_only$Part2),]

#sorting lncRNA - ENSG00000262202.2
LRep2_In_RNU42_only[3, "Part2"] <- "lnc-AC007952.2-2"

#ENSG00000262074.3
LRep2_In_RNU42_only[13, "Part2"] <- "lnc-GRAP-1"

#ENSG00000265185.1
LRep2_In_RNU42_only[17, "Part2"] <- "lnc-AC007952.1.1-1"

   
RNAInter_file_path <- "large_data/RNAInter_full_interactions_data.txt"

# Read the tab-delimited file
RNAInter <- read.delim(RNAInter_file_path)

#Only human
RNAInter <- RNAInter[grepl("Homo sapiens",RNAInter$Species1),]

#RNU4-1 in any column
RNAInter_RNU41 <- RNAInter %>%
  filter(apply(., 1, function(row) any(grepl("RNU4-1",row))))

#RNU4-2 in any column

RNAInter_RNU42 <- RNAInter %>%
  filter(apply(., 1, function(row) any(grepl("RNU4-2",row))))


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
colnames(NPinter) <- c("M1 Interaction ID", "Molecule 1", 
                       "Gene identifier", "RNA catagory",
                       "Molecule 2", "M2 Interaction ID", 
                       "Molecule type", "Data source", 
                       "Assay type", "database ID", "Species",
                       "Cell type", "Interaction catagory", 
                       "Interaction class", "Interaction catagory 2", "Source")

#Human only
NPinter <- NPinter[grepl("Homo sapiens",NPinter$Species),]

#RNU4-1 in any column
NPinter_RNU41 <- NPinter %>%
  filter(apply(., 1, function(row) any(grepl("RNU4-1",row))))


#RNU4-2 in any column

NPinter_RNU42 <- NPinter %>%
  filter(apply(., 1, function(row) any(grepl("RNU4-2",row))))

#RNU4-2 in the same row
NPinter_RNU42 <- NPinter_RNU42 %>%
  mutate(
    `Molecule 2` = ifelse(`Molecule 2` == "RNU4-2", `Molecule 1`, `Molecule 2`),
    `Molecule 1` = ifelse(`Molecule 1` != "RNU4-2", "RNU4-2", `Molecule 1`)
  )

#RNU4-1 in the same row
NPinter_RNU41 <- NPinter_RNU41 %>%
  mutate(
    `Molecule 2` = ifelse(`Molecule 2` == "RNU4-1", `Molecule 1`, `Molecule 2`),
    `Molecule 1` = ifelse(`Molecule 1` != "RNU4-1", "RNU4-1", `Molecule 1`)
  )

#Whats exlusive to RNU4-2?
NPinter_RNU42_only <- anti_join(NPinter_RNU42, NPinter_RNU41, by='Molecule 2') 


