# Specify the file path
file_path <- "/home/jessicaaucott/Desktop/LIGR-seq/GSM2113739_S1_AMT_Ligase_Rep1_v19-48.rmsk.uniq.unfilt.lig.txt"

# Read the tab-delimited file
data <- read.delim(file_path)

# View the first few rows of the data
head(data)

Filtered_RNU41 <- data[grepl("RNU4-1",data$TET1.TET1),]

Filtered_RNU42 <- data[grepl("RNU4-2",data$TET1.TET1),]


Only_Intermolecular_RNU41 <- Filtered_RNU41[!grepl("S", Filtered_RNU41$S),]

Only_Intermolecular_RNU42 <- Filtered_RNU42[!grepl("S", Filtered_RNU42$S),]


install.packages("compare")
library(compare)

install.packages("dplyr")
library(dplyr)

Only_Intermolecular_RNU41_TEST

In_RNU42_only <- anti_join(Only_Intermolecular_RNU42_Sep_Change, Only_Intermolecular_RNU41_Sep_Change, by='Part2') 

Only_Intermolecular_RNU42_Sep_TEST <- Only_Intermolecular_RNU42_Sep



Only_Intermolecular_RNU41_Sep_Change <- Only_Intermolecular_RNU41_Sep %>%
  mutate(
    Part2 = ifelse(Part2 == "RNU4-1", Part1, Part2),
    Part1 = ifelse(Part1 != "RNU4-1", "RNU4-1", Part1)
  )

Only_Intermolecular_RNU42_Sep_Change <- Only_Intermolecular_RNU42_Sep %>%
  mutate(
    Part2 = ifelse(Part2 == "RNU4-2", Part1, Part2),
    Part1 = ifelse(Part1 != "RNU4-2", "RNU4-2", Part1)
  )

#Missing Value = ENSG00000258486.2
Only_Intermolecular_RNU42_Sep_Change[15, "Part2"] <- "lnc-LRR1-1"




print(Only_Intermolecular_RNU42_Sep_Change)

rm(Only_Intermolecular_RNU42_Sep_Change)

install.packages("tidyr")
library(tidyr)

Only_Intermolecular_RNU42_Sep <- Only_Intermolecular_RNU42 %>%
  separate(TET1.TET1, into = c("Part1", "Part2"), sep = ":")



install.packages(c("cli", "glue", "lifecycle", "magrittr", "pillar", "rlang", "tibble", "tidyselect", "vctrs"))

rm(df)



Rep2_file_path <- "/home/jessicaaucott/Desktop/LIGR-seq/GSM2113743_S5_AMT_Ligase_Rep2_v19-48.rmsk.uniq.unfilt.lig.txt"

# Read the tab-delimited file
Rep2 <- read.delim(Rep2_file_path)

#Only RNU4-1 and -2 interactions
Rep2_Filtered_RNU41 <- Rep2[grepl("RNU4-1",Rep2$XLOC.010820.XLOC.010820),]

Rep2_Filtered_RNU42 <- Rep2[grepl("RNU4-2",Rep2$XLOC.010820.XLOC.010820),]

#Only intermolecular interactions
Rep2_Only_Intermolecular_RNU41 <- Rep2_Filtered_RNU41[!grepl("S", Rep2_Filtered_RNU41$S),]

Rep2_Only_Intermolecular_RNU42 <- Rep2_Filtered_RNU42[!grepl("S", Rep2_Filtered_RNU42$S),]

#Seperate the interaction columns
Rep2_Only_Intermolecular_RNU41_Sep <- Rep2_Only_Intermolecular_RNU41 %>%
  separate(XLOC.010820.XLOC.010820, into = c("Part1", "Part2"), sep = ":")

Rep2_Only_Intermolecular_RNU42_Sep <- Rep2_Only_Intermolecular_RNU42 %>%
  separate(XLOC.010820.XLOC.010820, into = c("Part1", "Part2"), sep = ":")

#RNU4 in the same row
Rep2_Only_Intermolecular_RNU41_Sep_Change <- Rep2_Only_Intermolecular_RNU41_Sep %>%
  mutate(
    Part2 = ifelse(Part2 == "RNU4-1", Part1, Part2),
    Part1 = ifelse(Part1 != "RNU4-1", "RNU4-1", Part1)
  )

Rep2_Only_Intermolecular_RNU42_Sep_Change <- Rep2_Only_Intermolecular_RNU42_Sep %>%
  mutate(
    Part2 = ifelse(Part2 == "RNU4-2", Part1, Part2),
    Part1 = ifelse(Part1 != "RNU4-2", "RNU4-2", Part1)
  )

#Whats exlusive to RNU4-2?
Rep2_In_RNU42_only <- anti_join(Rep2_Only_Intermolecular_RNU42_Sep_Change, Rep2_Only_Intermolecular_RNU41_Sep_Change, by='Part2') 

#sorting lncRNA - ENSG00000262202.2
Rep2_In_RNU42_only[15, "Part2"] <- "lnc-AC007952.2-2"

#ENSG00000262074.3
Rep2_In_RNU42_only[57, "Part2"] <- "lnc-GRAP-1"

#ENSG00000265185.1
Rep2_In_RNU42_only[74, "Part2"] <- "lnc-AC007952.1.1-1"

#Whats in Rep 2 but not 1
Rep2_Exclusive<- anti_join(Rep2_In_RNU42_only, In_RNU42_only, by='Part2') 

#what values are in 1 and 2
common_values <- intersect(Rep2_In_RNU42_only$Part2, In_RNU42_only$Part2)
print(common_values)

common_rows <-merge(Rep2_In_RNU42_only, In_RNU42_only, by = "Part2")

rm


Rep1_background_file_path <- "/home/jessicaaucott/Desktop/LIGR-seq/GSM2113741_S3_no-AMT_Ligase_Rep1_v19-48.rmsk.uniq.unfilt.lig.txt"

# Read the tab-delimited file
NO_AMT_REP1 <- read.delim(Rep1_background_file_path)

#Only intermolecular interactions

Rep1_Background_RNU41_intermolecular <- Rep1_Background_filtered_RNU41[!grepl("S", Rep1_Background_filtered_RNU41$S),]

Rep1_Background_RNU42_intermolecular <- Rep1_Background_filtered_RNU42[!grepl("S", Rep1_Background_filtered_RNU41$S),]

#Seperate the interaction columns
Rep1_Background_RNU41_intermolecular_Sep <- Rep1_Background_RNU41_intermolecular %>%
  separate(Y.RNA.Y.RNA, into = c("Part1", "Part2"), sep = ":")

Rep1_Background_RNU42_intermolecular_Sep <- Rep1_Background_RNU42_intermolecular %>%
  separate(Y.RNA.Y.RNA, into = c("Part1", "Part2"), sep = ":")

#RNU4 in the same row
Rep1_Background_RNU41_intermolecular_Sep_Change <- Rep1_Background_RNU41_intermolecular_Sep %>%
  mutate(
    Part2 = ifelse(Part2 == "RNU4-1", Part1, Part2),
    Part1 = ifelse(Part1 != "RNU4-1", "RNU4-1", Part1)
  )

Rep1_Background_RNU42_intermolecular_Sep_Change <- Rep1_Background_RNU42_intermolecular_Sep %>%
  mutate(
    Part2 = ifelse(Part2 == "RNU4-2", Part1, Part2),
    Part1 = ifelse(Part1 != "RNU4-2", "RNU4-2", Part1)
  )

Remove_background_rep1_RNU41 <- anti_join(Only_Intermolecular_RNU42_Sep_Change, Rep1_Background_RNU42_intermolecular_Sep_Change, by='Part2')


#Remove RNU6

In_RNU42_only <- In_RNU42_only[!grepl("RNU6", In_RNU42_only$Part2),]

Rep2_Exclusive <- Rep2_Exclusive[!grepl("RNU6", Rep2_Exclusive$Part2),]

RNAInter_file_path <- "/home/jessicaaucott/Desktop/RNAinter database/Download_data_RR.txt"

# Read the tab-delimited file
RNAInter <- read.delim(RNAInter_file_path)

#Only human
RNAInter_Human <- RNAInter[grepl("Homo sapiens",RNAInter$Species1),]

#RNU4-1 in any column

RNAInter_RNU41 <- RNAInter_Human %>%
  filter(apply(., 1, function(row) any(grepl("RNU4-1",row))))


#RNU4-2 in any column

RNAInter_RNU42 <- RNAInter_Human %>%
  filter(apply(., 1, function(row) any(grepl("RNU4-2",row))))

#Remove pseudogenes
value_to_remove <- "RNU4-23P"
# Remove rows where any column contains the value 
RNAInter_RNU42 <- RNAInter_RNU42 %>% 
  filter(!apply(., 1, function(row) any(row == value_to_remove)))


#RNU4-2 in the same row
RNAInter_RNU42_change <- RNAInter_RNU42 %>%
  mutate(
    Interactor2.Symbol = ifelse(Interactor2.Symbol == "RNU4-2", Interactor1.Symbol, Interactor2.Symbol),
    Interactor1.Symbol = ifelse(Interactor1.Symbol != "RNU4-2", "RNU4-2", Interactor1.Symbol)
  )

#RNU4-1 in the same row
RNAInter_RNU41_change <- RNAInter_RNU41 %>%
  mutate(
    Interactor2.Symbol = ifelse(Interactor2.Symbol == "RNU4-1", Interactor1.Symbol, Interactor2.Symbol),
    Interactor1.Symbol = ifelse(Interactor1.Symbol != "RNU4-1", "RNU4-1", Interactor1.Symbol) 
  )

#Whats exlusive to RNU4-2?
RNAInter_RNU42_only <- anti_join(RNAInter_RNU42_change, RNAInter_RNU41_change, by='Interactor2.Symbol') 

rm(RNAInter)

write.table(RNAInter_Human, file = "/home/jessicaaucott/Desktop/RNAinter_database.txt", sep = "\t", row.names = FALSE)\

#NPinter interactions

NPinter_file_path <- "/home/jessicaaucott/Desktop/RNAinter database/snRNA_interaction.txt"

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


#SAving RNU4-2 exclusive data frames
write.table(In_RNU42_only, file = "/home/jessicaaucott/Desktop/RNU4_Interactions_Project/LIGR-seq_Rep1_Processed_RNU42exl.txt", sep = "\t", row.names = FALSE)

write.table(Rep2_In_RNU42_only, file = "/home/jessicaaucott/Desktop/RNU4_Interactions_Project/LIGR-seq_Rep2_Processed_RNU42exc.txt", sep = "\t", row.names = FALSE)

write.table(NPinter_RNU42_only, file = "/home/jessicaaucott/Desktop/RNU4_Interactions_Project/NPinterdata_Processed_RNU42exc.txt", sep = "\t", row.names = FALSE)

write.table(RNAInter_RNU42_only, file = "/home/jessicaaucott/Desktop/RNU4_Interactions_Project/RNAinterdata_Processed_RNU42exc.txt", sep = "\t", row.names = FALSE)
