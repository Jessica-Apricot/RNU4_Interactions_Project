row.names(All_RNU4) <- NULL

#Fix the encode values
All_RNU4[176, 2] <- "lnc-AC007952.1.1-1"
All_RNU4[86, 2] <- "lnc-LRR1-1"
All_RNU4[141, 2] <- "lnc-AC007952.2-2"
All_RNU4[165, 2] <- "lnc-GRAP-1"

#Remove the duplicate values
All_RNU4 <- All_RNU4[-207, ]
All_RNU4 <- All_RNU4[-61, ]


row.names(All_RNU4) <- NULL

# Create the RNA vector with no extra commas
RNA_vector <- c(
  "snRNA", "snRNA", "srpRNA", "snRNA", "snRNA", "srpRNA", "snRNA", "lincRNA", "snRNA", "snRNA", "snRNA",
  "snRNA", "rRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA",
  "snRNA", "snRNA", "snRNA", "snRNA", "snoRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA",
  "snRNA", "snRNA", "snRNA", "rRNA", "snRNA", "snRNA", "snRNA", "srpRNA", "snRNA", "snRNA", "snRNA", "snRNA", 
  "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "mRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snoRNA", 
  "snRNA", "mRNA", "snRNA", "mRNA", "mRNA", "mRNA", "rRNA", "ncRNA", "lincRNA", "lincRNA", "lincRNA", "rRNA",
  "snoRNA", "srpRNA", "snRNA", "snRNA", "srpRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA",
  "rRNA", "snRNA", "lincRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "srpRNA", "snRNA", "snRNA",
  "rRNA", "srpRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA",
  "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snoRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA",
  "srpRNA" , "snRNA", "snRNA", "rRNA", "snRNA", "snRNA", "snRNA", "snRNA", "srpRNA", "snRNA", "snRNA", "snRNA",
  "snRNA", "snRNA", "snoRNA", "snRNA", "snRNA", "snRNA", "srpRNA", "snRNA", "snRNA", "lincRNA", "snRNA", "snoRNA",
  "snRNA", "snRNA", "snoRNA", "rRNA", "snRNA", "snRNA", "snoRNA", "rRNA", "snRNA", "snRNA", "snRNA", "snoRNA",
  "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "lincRNA", "snRNA", "snRNA", 
  "snRNA", "snRNA", "snRNA", "snRNA", "snoRNA", "snRNA", "snRNA", "snRNA", "lincRNA", "snRNA", "snRNA", "snRNA",
  "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "snRNA", "rRNA", "mRNA", "mRNA", "mRNA",
  "mRNA", "mRNA", "mRNA", "mRNA", "mRNA", "mRNA", "ncRNA", "mRNA", "mRNA", "mRNA", "mRNA", "mRNA", "mRNA",
  "lincRNA", "mRNA",  "lincRNA", "mRNA", "rRNA", "mRNA", "snoRNA", "mRNA", "mRNA", "snRNA", "mRNA"  
)

# Add the RNA_vector as a new column to the All_RNU4 dataframe
All_RNU4$RNA <- RNA_vector


write.csv(All_RNU4, "/home/jessicaaucott/Desktop/RNU4_Interactions.csv", row.names = FALSE)




rRNA

  