install.packages("dplyr")
library(dplyr)

install.packages("tidyr")
library(tidyr)

# Specify the file path
LIGR_seq_Rep1_file_path <- "large_data/GSM2113739_S1_AMT_Ligase_Rep1_v19-48.rmsk.uniq.unfilt.lig.txt"

# Read the tab-delimited file
LIGR_seq_Rep1 <- read.delim(LIGR_seq_file_path)

# View the first few rows of the data
head(LIGR_seq_Rep1)

#Filter out only rows with RNU4
LRep1_Filtered_RNU41 <- LIGR_seq_Rep1[grepl("RNU4-1",LIGR_seq_Rep1$TET1.TET1),]

LRep1_Filtered_RNU42 <- LIGR_seq_Rep1[grepl("RNU4-2",LIGR_seq_Rep1$TET1.TET1),]

#only RNU4 that is interacting with molecules other than self
LRep1_Only_Intermolecular_RNU41 <- LRep1_Filtered_RNU41[!grepl("S", LRep1_Filtered_RNU41$S),]

LRep1_Only_Intermolecular_RNU42 <- LRep1_Filtered_RNU42[!grepl("S", LRep1_Filtered_RNU42$S),]

#Separate the interactions into two columns
LRep1_Only_Intermolecular_RNU41 <- LRep1_Only_Intermolecular_RNU41 %>%
  separate(TET1.TET1, into = c("Part1", "Part2"), sep = ":")

LRep1_Only_Intermolecular_RNU42 <- LRep1_Only_Intermolecular_RNU42 %>%
  separate(TET1.TET1, into = c("Part1", "Part2"), sep = ":")

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

#Fixing lncRNA values = ENSG00000258486.2
LRep1_Only_Intermolecular_RNU42[15, "Part2"] <- "lnc-LRR1-1"

#Values only present in RNU4-2
LRep1_In_RNU42_only <- anti_join(LRep1_Only_Intermolecular_RNU42, LRep1_Only_Intermolecular_RNU41, by = "Part2") 

#Remove U6 interactions
LRep1_In_RNU42_only <- LRep1_In_RNU42_only[!grepl("U6",LRep1_In_RNU42_only$Part2),]



LIGR_seq_Rep2_file_path <- "large_data/GSM2113743_S5_AMT_Ligase_Rep2_v19-48.rmsk.uniq.unfilt.lig.txt"

# Read the tab-delimited file
LIGR_seq_Rep2 <- read.delim(LIGR_seq_Rep2_file_path)

#Only RNU4-1 and -2 interactions
LRep2_Filtered_RNU41 <- LIGR_seq_Rep2[grepl("RNU4-1",LIGR_seq_Rep2$XLOC.010820.XLOC.010820),]

LRep2_Filtered_RNU42 <- LIGR_seq_Rep2[grepl("RNU4-2",LIGR_seq_Rep2$XLOC.010820.XLOC.010820),]

#only RNU4 that is interacting with molecules other than self
LRep2_Only_Intermolecular_RNU41 <- LRep2_Filtered_RNU41[!grepl("S", LRep2_Filtered_RNU41$S),]

LRep2_Only_Intermolecular_RNU42 <- LRep2_Filtered_RNU42[!grepl("S", LRep2_Filtered_RNU42$S),]

#Separate the interactions into two columns
LRep2_Only_Intermolecular_RNU41 <- LRep2_Only_Intermolecular_RNU41 %>%
  separate(XLOC.010820.XLOC.010820, into = c("Part1", "Part2"), sep = ":")

LRep2_Only_Intermolecular_RNU42 <- LRep2_Only_Intermolecular_RNU42 %>%
  separate(XLOC.010820.XLOC.010820, into = c("Part1", "Part2"), sep = ":")

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
LRep2_In_RNU42_only[15, "Part2"] <- "lnc-AC007952.2-2"

#ENSG00000262074.3
LRep2_In_RNU42_only[57, "Part2"] <- "lnc-GRAP-1"

#ENSG00000265185.1
LRep2_In_RNU42_only[74, "Part2"] <- "lnc-AC007952.1.1-1"



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
RNAInter_RNU42_only <- anti_join(RNAInter_RNU42_change, RNAInter_RNU41_change, by='Interactor2.Symbol') 


#NPinter interactions

NPinter_file_path <- "large_data/snRNA_interaction.txt"

# Read the tab-delimited file
NPinter <- read.delim(NPinter_file_path)

#Human only
NPinter <- NPinter[grepl("Homo sapiens",NPinter$Homo.sapiens),]

#RNU4-1 in any column
NPinter_RNU41 <- NPinter %>%
  filter(apply(., 1, function(row) any(grepl("RNU4-1",row))))


#RNU4-2 in any column

NPinter_RNU42 <- NPinter %>%
  filter(apply(., 1, function(row) any(grepl("RNU4-2",row))))

#RNU4-2 in the same row
NPinter_RNU42_change <- NPinter_RNU42 %>%
  mutate(
    C21orf33 = ifelse(C21orf33 == "RNU4-2", U11, C21orf33),
    U11 = ifelse(U11 != "RNU4-2", "RNU4-2", U11)
  )

#RNU4-1 in the same row
NPinter_RNU41_change <- NPinter_RNU41 %>%
  mutate(
    C21orf33 = ifelse(C21orf33 == "RNU4-1", U11, C21orf33),
    U11 = ifelse(U11 != "RNU4-1", "RNU4-1", U11)
  )

#Whats exlusive to RNU4-2?
NPinter_RNU42_only <- anti_join(NPinter_RNU42_change, NPinter_RNU41_change, by='C21orf33') 


