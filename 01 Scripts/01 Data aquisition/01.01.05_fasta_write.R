###Fasta file script###

library(seqinr)

#csv <- "20230322_AT_SMILES_unified.csv"
csv <- "YYYYMMDD_KS_SMILES_unified.csv"
#csv <- "C:\\Users\\q31032mw\\Dropbox (The University of Manchester)\\Max\\17_ML_Project\\Webscraping\\ClusterCAD Webscraping\\20231107_Smiles_and_AT_reviewed_and_unreviewed_extras.csv"
sequences_file <- read.csv(csv)
sequences_file <- sequences_file[!is.na(sequences_file$SMILES),]

# As monomers #

for (i in 1:nrow(sequences_file)) {try(
  
  write.fasta(sequences  = sequences_file$Sequence[i], file.out = paste0("C:\\Users\\q31032mw\\Dropbox (The University of Manchester)\\Max\\17_ML_Project\\Webscraping\\ClusterCAD Webscraping\\KS_extras\\", sequences_file$Gene[i], "_SNAKE.fasta"), open = "w", names = sequences_file$Gene[i])
)}


# As dimers # 

for (i in 1:nrow(sequences_file)) {
  ii <- paste0(sequences_file$Sequence[i], ":", sequences_file$Sequence[i])
  print(ii)
  write.fasta(sequences  = ii, file.out = paste0("C:\\Users\\q31032mw\\Dropbox (The University of Manchester)\\Max\\17_ML_Project\\Webscraping\\ClusterCAD Webscraping\\KS_extras_dimers\\", sequences_file$Gene[i], "_dimer.fasta"), open = "w", names = sequences_file$Gene[i])
}
