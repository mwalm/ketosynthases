###KS domain extraction script###
##Function of this script is to extract KS domain sequences from mastersheet output of 20230302_Looped_MiBIG_scrape.R##
##Second loop combines these with SMILES data and information on KS type##

##Dependencies:##
library(dplyr)

##Inputs:##
mastersheet <- read.csv("yymmdd_mastersheet.csv")
SMILES_sheet <- read.csv("YYYYMMDD_Smiles_and_AT_reviewed.csv")

#KS extraction loop
#filter for KS columns#
x3 <- c(2)
for (i in 0:9) {
  x1 <- paste0("KS.", i)
  x2 <- grepl(pattern = x1, x = colnames(mastersheet))
  
  for(ii in 1:length(x2)){
    if(x2[ii] == TRUE){x3 <- c(x3, ii)}
  }
}
KS_df <- mastersheet[,x3]

#remove AT-less rows#
empty_rows <- c()
for (i in 1:nrow(KS_df)){
  if(sum(is.na(KS_df[i,2:10])) == 9){empty_rows <- c(empty_rows, i)}
}
KS_df <- KS_df[-empty_rows,]
print.noquote(paste0("Total KS domains: ", sum(!is.na(KS_df))))

#Break KSs to modules to line up with the smiles data. 
KS_SMILE_df <- data.frame(NA, NA)
colnames(KS_SMILE_df) <- c("Gene", "Sequence")
for (i in 1:nrow(KS_df)) {
  for (ii in 2:ncol(KS_df)){
    if(is.na(KS_df[i,ii]) == FALSE){
      KS_SMILE_df <- rbind(KS_SMILE_df, c(paste0(KS_df[i,1],"_Mod.", ii - 2), KS_df[i,ii]))
    }
  }  
}
KS_SMILE_df <- KS_SMILE_df[-1,]
KS_SMILE_df$SMILES <- NA


#remove duplicates
KS_SMILE_df <- KS_SMILE_df[!duplicated(KS_SMILE_df$Gene),]

#identify genes with a module 0 in the MiBIG group#
modzeros <- c()
for (i in 1:nrow(KS_SMILE_df)) {
  if(sum(grepl(pattern = paste0(sub("\\..*", "", x = KS_SMILE_df[i, 1]), ".0"), KS_SMILE_df$Gene)) >= 1)
  {modzeros <- c(modzeros, KS_SMILE_df[i, 1])}
  
}

#for genes not in modzeros, add 1 to the smiles_sheet[,6]
j <-  sub("\\..*", "", SMILES_sheet[grep("\\.0", SMILES_sheet$renamed_modules),6])
for(i in j){
  Addto<- grep(pattern = i, x = SMILES_sheet$renamed_modules)
  SMILES_sheet[Addto, 6] <- paste0(sub("\\..*", "", x = SMILES_sheet[Addto,6]), ".", (as.integer(sub(".*\\.", "", SMILES_sheet[Addto,6])) + 1))
}

#Replace the missing module 0s. Note: script is not perfect.
SMILES_sheet2 <- read.csv("YYYYMMDD_Smiles_and_AT_reviewed.csv")
for(i in modzeros){
  ii <- grep(pattern = i, SMILES_sheet2$renamed_modules)
  SMILES_sheet[ii, 6] <- SMILES_sheet2[ii, 6]
}

#combine elements from KS sequence and SMILES information
for (i in 1:nrow(KS_SMILE_df)) {
  ii <- grep(KS_SMILE_df[i, 1] , x = SMILES_sheet$renamed_modules)
  iii <- 0
  KS_SMILE_df[i,3] <- SMILES_sheet[(ii[1] - iii), 3]
}

KS_SMILE_df$KS_substrate <- NA
for (i in 1:nrow(KS_SMILE_df)) {
  ii <- grep(KS_SMILE_df[i, 1] , x = SMILES_sheet$renamed_modules)
  iii <- 0
  KS_SMILE_df[i,4] <- SMILES_sheet[(ii[1] - iii), 5]
}

KS_SMILE_df$BGC_compound <- NA
for (i in 1:nrow(KS_SMILE_df)) {
  ii <- grep(KS_SMILE_df[i, 1] , x = SMILES_sheet$renamed_modules)
  iii <- 0
  KS_SMILE_df[i,5] <- SMILES_sheet[(ii[1] - iii), 4]
}

#write to file:
write.csv(KS_SMILE_df, file = paste0(date, "_KS_SMILES_unified.csv"))

          