#### Script for Ketosynthases ####

### B-score graph ###
# Install and load the msa package
library(msa)
library(Biostrings)
library(bio3d)

## Define ML output file location and names ##
date <- "20240411"
experiment_type <- "_raw_KRvsDH_binary_GCN_"
output_folder <- "/home/q31032mw/Dropbox (The University of Manchester)/Max/17_ML_Project/Figures scripts/20240411_true_KRvsDH_GCN_binary/"
#output_folder <- "C:\\Users\\q31032mw\\Dropbox (The University of Manchester)\\Max\\17_ML_Project\\Figures scripts\\"

## amino acid sequences import ## 
route <-  "/home/q31032mw/Dropbox (The University of Manchester)/Max/17_ML_Project/Rank_1/20231211_KS_dimers_aligned/"
#route <- "C:\\Users\\q31032mw\\Dropbox (The University of Manchester)\\Max\\17_ML_Project\\Rank_1\\20231211_KS_dimers_aligned\\"

# Extract file names in route # 
files_list <- list.files(path = route)

# NetworkX coordinates ## 
netX_coords <- read.csv("/home/q31032mw/Dropbox (The University of Manchester)/Max/17_ML_Project/Figures scripts/20240306_binary_classification_KRvsDH_graphviz.csv")
#netX_coords <- read.csv("C:\\Users\\q31032mw\\Dropbox (The University of Manchester)\\Max\\17_ML_Project\\Figures scripts\\20240208_KR_DH_graphviz.csv")

# SeEryAIII_Mod.1
ery2 <- paste0(route, "nysC_Mod.2.pdb_aligned.pdb")

# Define your reference sequence and other sequences # 
# You can change this if you need to do per-type alignments. See the 
reference_sequence_NR2 <- read.pdb(paste0(route, "NcmAI_Mod.3.pdb_aligned.pdb")) 
reference_sequence_NR <- pdbseq(reference_sequence_NR2)

#reference_sequence_KR2 <- read.pdb(paste0(route, "eryAI_Mod.2.pdb_aligned.pdb")) 
#reference_sequence_KR2 <- read.pdb(paste0(route, "scnS2_Mod.2.pdb_aligned.pdb")) 

reference_sequence_KR2 <- read.pdb(ery2) 
reference_sequence_KR <- pdbseq(reference_sequence_KR2)

#reference_sequence_DH2 <- read.pdb(paste0(route, "nysC_Mod.2.pdb_aligned.pdb")) 
reference_sequence_DH2 <- read.pdb(ery2) 
reference_sequence_DH <- pdbseq(reference_sequence_DH2)

reference_sequence_ER2 <- read.pdb(paste0(route, "eryAIII_Mod.1.pdb_aligned.pdb")) 
reference_sequence_ER <- pdbseq(reference_sequence_ER2)

# Create data structures for B-score and residue #
total_data <- list()
mat1_data <- list()
for (i in files_list){
  print(i)
  name_in_list <- paste0(route, i)
  pdb <- read.pdb(name_in_list)
  sequence <- pdbseq(pdb, aa1 = TRUE )
  mat1 <- matrix()
  pdb2 <- cbind(pdb[1], as.data.frame(pdb[3]))
  colnames(pdb2) <- gsub(x = colnames(pdb2),pattern = "atom.", replacement = "") 
  total_data[[i]] <- pdb2
  b_values <- pdb2[["b"]][which(pdb2[["elety"]] == "CA")]
  mat1 <- cbind(sequence, b_values)
  mat1_data[[i]] <- mat1
}

# apply netX as column #
mat1_data_netX <- c()
for (i in netX_coords$Gene..AT.) {
  index <- grep(pattern = i, x = names(mat1_data))
  mat1_data_netX <- c(mat1_data_netX, index)
}
countr <- 0L
mat3_data <- mat1_data[mat1_data_netX]
for(i in names(mat3_data)){
  countr <- countr + 1L
  strong <- as.numeric(unlist(strsplit(netX_coords$Residues.strong[countr], ",")))
  weak <-   as.numeric(unlist(strsplit(netX_coords$Residues.weak[countr], ",")))
  
  # Create vectors filled with zeros
  strong_nodes <- rep(0, nrow(mat3_data[[i]]))
  weak_nodes <- rep(0, nrow(mat3_data[[i]]))
  combined_nodes <- rep(0, nrow(mat3_data[[i]]))
  # Set 1s in rows designated by strong values
  strong_nodes[strong] <- 100L
  # Set 0.5s in rows designated by weak values
  weak_nodes[weak] <- 50L
  # Create a combined column # 
  combined_nodes[strong] <- 100L
  combined_nodes[weak] <- 50L
  
  # Append new columns to the data frame
  mat3_data[[i]] <- cbind(mat3_data[[i]], strong_nodes, weak_nodes, combined_nodes)
  
  # Rename the new columns
  colnames(mat3_data[[i]])[ncol(mat3_data[[i]]) - 2] <- "strong_nodes"
  colnames(mat3_data[[i]])[ncol(mat3_data[[i]]) - 1] <- "weak_nodes"
  colnames(mat3_data[[i]])[ncol(mat3_data[[i]])] <- "combined_nodes"
  print(i)
}


# Adjust netX coordinates to account for dimer #
countr <- 0
mat4_data <- mat3_data
for (i in mat3_data) {
  countr <- countr + 1L
  n_res_mono <- (nrow(i)/2)
  row_coords <- which(i[,5] != 0)
  
  positive <- c()
  negative <- c()
  for(ii in row_coords){
    positive <- c(positive, (ii + n_res_mono))
    negative <- c(negative, (ii - n_res_mono))
  }
  positive <- positive[which(positive <= nrow(i))]
  negative <- negative[which(negative >= 0L)]
  combined_posneg <- c(positive, negative)
  i[combined_posneg, 5] <- 100L
  print(mat3_data[countr])
  mat4_data[[countr]] <- i
}
mat3_data <- mat4_data



Substrates <-  read.csv("/home/q31032mw/Dropbox (The University of Manchester)/Max/17_ML_Project/Rank_1/20231108_KS_expandedset_reduction_state_labels_new_borders")
#Substrates <-  read.csv("C:\\Users\\q31032mw\\Dropbox (The University of Manchester)\\Max\\17_ML_Project\\Rank_1\\20231108_KS_expandedset_reduction_state_labels_new_borders")

NR_index <- which(Substrates$NR == 1L)
KR_index <- which(Substrates$KR == 1L)
DH_index <- which(Substrates$DH == 1L)
ER_index <- which(Substrates$ER == 1L)

NR_names <- Substrates$X[NR_index]
KR_names <- Substrates$X[KR_index]
DH_names <- Substrates$X[DH_index]
ER_names <- Substrates$X[ER_index]




## B_score indexing ##
NR_names_index <- c()
for (i in NR_names) {
  print(i)
  index <- grep(x = names(mat1_data), pattern = i)
  if (length(index) != 0) {
    NR_names_index <- c(NR_names_index, index)
  }
}

KR_names_index <- c()
for (i in KR_names) {
  print(i)
  index <- grep(x = names(mat1_data), pattern = i)
  if (length(index) != 0) {
    KR_names_index <- c(KR_names_index, index)
  }
}

DH_names_index <- c()
for (i in DH_names) {
  print(i)
  index <- grep(x = names(mat1_data), pattern = i)
  if (length(index) != 0) {
    DH_names_index <- c(DH_names_index, index)
  }
}

ER_names_index <- c()
for (i in ER_names) {
  print(i)
  index <- grep(x = names(mat1_data), pattern = i)
  if (length(index) != 0) {
    ER_names_index <- c(ER_names_index, index)
  }
}




# Function for row insertion as part of alignment #
insert_rows <- function(df, rows_to_insert, ref_seq) {
  #rows_to_insert <- which(rows_to_insert <= ref_seq)
  for (row_num in rows_to_insert) {
    # Add a new row with NA values
    new_row <- c(sequence = NA, b_values = NA)
    if(row_num >= nrow(df)){
      row_num <- (as.numeric(row_num) - 1)
      print("BLAM")
    }
    df <- rbind(df[1:(row_num - 1), ], new_row, df[row_num:nrow(df), ])
    
    # Update row numbers
    rownames(df) <- NULL
  }
  return(df)
}

## Perform alignment, add and remove rows from query sequences ##
alignment_function <- function(mat1_data, names_index, reference_sequence){
  countr <- 0
  namez <- names(mat1_data[names_index])
  mat2_data <- list()
  for (i in mat1_data[names_index]) {
    alignment <- pairwiseAlignment(
      subject = paste(reference_sequence, collapse = ""),
      pattern = paste(i[,1], collapse = ''),
      type = "local-global"
    )
    alignment2 <- pairwiseAlignment(
      pattern = paste(reference_sequence, collapse = ""),
      subject = paste(i[,1], collapse = ''),
      type = "global-local"
    )
    countr <- countr + 1
    name <- namez[countr]
    print(name)
    
    input_mat <- insert_rows(
      i, 
      which(as.matrix(alignment) == "-"), 
      length(reference_sequence)
    )
    mat2_data[[name]] <- input_mat
    
    if(length(grep(x = as.matrix(alignment2), pattern = "-")) == 0){next}
    
    else{mat2_data[[name]] <- mat2_data[[name]][-grep(x = as.matrix(alignment2), pattern = "-"),]}
    
    if(nrow(mat2_data[[name]]) >= length(reference_sequence)){
      mat2_data[[name]] <- mat2_data[[name]][1:length(reference_sequence),]
    
    }
  }
  return(mat2_data)
}
## Alignments for B-score ##
# NR # 
NR_out<- alignment_function(mat1_data = mat1_data, names_index = NR_names_index, reference_sequence = reference_sequence_NR)
# KR # 
KR_out<- alignment_function(mat1_data = mat1_data, names_index = KR_names_index, reference_sequence = reference_sequence_KR)
# DH # 
DH_out<- alignment_function(mat1_data = mat1_data, names_index = DH_names_index, reference_sequence = reference_sequence_DH)
# ER # 
ER_out<- alignment_function(mat1_data = mat1_data, names_index = ER_names_index, reference_sequence = reference_sequence_ER)




## B-score averages ## 
b_score_avg_function <- function(alignment_out, reference_sequence_2 ){
  b_scores <- matrix(rep(1:nrow(alignment_out[[1]])))
  for (i in alignment_out) {
    b_scores <- cbind(b_scores, i[1:nrow(alignment_out[[1]]),2])
  }
  b_scores <- b_scores[,-1]
  # Coerce each element to numeric
  b_scores_numeric <- apply(b_scores, 2, as.numeric)
  
  # Calculate row means while ignoring NA values
  b_average <- rowMeans(b_scores_numeric, na.rm = TRUE)
  
  atom_map <- reference_sequence_2[["atom"]][["resno"]]
  atom_matrix <- matrix(nrow = length(atom_map), ncol = 2)
  atom_matrix[,1] <- atom_map
  atom_matrix[,2] <- b_average[atom_matrix[,1]]
  output_b_score <- reference_sequence_2
  output_b_score[["atom"]][["b"]] <- atom_matrix[,2]
  return(output_b_score)
}



# NR b_avg write # 
NR_b_score_avg <- b_score_avg_function(NR_out, reference_sequence_NR2)
write.pdb(pdb = NR_b_score_avg, file = paste0(output_folder, paste0(date, "_NR_b_score_NcmAI_Mod.3.pdb")) )
          
# KR b_avg write # 
KR_b_score_avg <- b_score_avg_function(KR_out, reference_sequence_KR2)
write.pdb(pdb = KR_b_score_avg, file = paste0(output_folder, paste0(date, "_KR_b_score_eryAI_Mod.2.pdb")) )

# DH b_avg write #
DH_b_score_avg <- b_score_avg_function(DH_out, reference_sequence_DH2)
write.pdb(pdb = DH_b_score_avg, file = paste0(output_folder, paste0(date, "_DH_b_score_nysC_Mod.2.pdb")) )

# ER b_avg write # 
ER_b_score_avg <- b_score_avg_function(KR_out, reference_sequence_ER2)
write.pdb(pdb = ER_b_score_avg, file = paste0(output_folder, paste0(date, "_ER_b_score_eryAIII_Mod.1.pdb")) )



## NetX averages ##
## Perform alignment, add and remove rows from query sequences ##
# NR in test set #
NR_test_set <- c()
for(i in NR_names){
  NR_test_set <- c(NR_test_set, grep(pattern = i, names(mat3_data)))
}

# KR in test set #
KR_test_set <- c()
for(i in KR_names){
  KR_test_set <- c(KR_test_set, grep(pattern = i, names(mat3_data)))
}

# DH in test set #
DH_test_set <- c()
for(i in DH_names){
  DH_test_set <- c(DH_test_set, grep(pattern = i, names(mat3_data)))
}

# KR in test set #
ER_test_set <- c()
for(i in ER_names){
  ER_test_set <- c(ER_test_set, grep(pattern = i, names(mat3_data)))
}

## Alignments for NetX ##
# NR # 
NR_out<- alignment_function(mat1_data = mat3_data, names_index = NR_test_set, reference_sequence = reference_sequence_NR)
# KR # 
KR_out<- alignment_function(mat1_data = mat3_data, names_index = KR_test_set, reference_sequence = reference_sequence_KR)
# DH # 
DH_out<- alignment_function(mat1_data = mat3_data, names_index = DH_test_set, reference_sequence = reference_sequence_DH)
# ER # 
ER_out<- alignment_function(mat1_data = mat3_data, names_index = ER_test_set, reference_sequence = reference_sequence_ER)


## NetX averages ## 
netx_avg_function <- function(alignment_out, reference_sequence_2 ){
  netx_values <- matrix(rep(1:nrow(alignment_out[[1]])))
  for (i in alignment_out) {
    netx_values <- cbind(netx_values, i[1:nrow(alignment_out[[1]]),5])
  }
  netx_values <- netx_values[,-1]
  # Coerce each element to numeric
  netx_values_numeric <- apply(netx_values, 2, as.numeric)
  
  # Calculate row means while ignoring NA values
  netx_average <- rowMeans(netx_values_numeric, na.rm = TRUE)
  
  atom_map <- reference_sequence_2[["atom"]][["resno"]]
  atom_matrix <- matrix(nrow = length(atom_map), ncol = 2)
  atom_matrix[,1] <- atom_map
  atom_matrix[,2] <- netx_average[atom_matrix[,1]]
  output_netx_score <- reference_sequence_2
  output_netx_score[["atom"]][["b"]] <- atom_matrix[,2]
  return(output_netx_score)
}

netx_avg_function(KR_out, reference_sequence_KR2)

# NR netx_avg write # 
NR_netx_score_avg <- netx_avg_function(NR_out, reference_sequence_NR2)
write.pdb(pdb = NR_netx_score_avg, file = paste0(output_folder, paste0(date, experiment_type ,"_NR_netX_score_NcmAI_Mod.3.pdb")) )

# KR netx_avg write # 
KR_netx_score_avg <- netx_avg_function(KR_out, reference_sequence_KR2)
write.pdb(pdb = KR_netx_score_avg, file = paste0(output_folder, paste0(date, experiment_type, "_KR_netx_score_SeEryAIII_Mod.2.pdb")) )


# DH netx_avg write # 
DH_netx_score_avg <- netx_avg_function(DH_out, reference_sequence_DH2)
write.pdb(pdb = DH_netx_score_avg, file = paste0(output_folder, paste0(date, experiment_type, "_DH_netx_score_SeEryAIII_Mod.2.pdb")) )

# ER netx_avg write # 
ER_netx_score_avg <- netx_avg_function(ER_out, reference_sequence_ER2)
write.pdb(pdb = ER_netx_score_avg, file = paste0(output_folder, paste0(date, experiment_type, "_ER_netx_score_eryAIII_Mod.1.pdb")) )


#####################
### Logo Diagrams ###
#####################

library(ggplot2)
library(ggseqlogo)

# list of residues from cross words # 
KRr1 <- 160:170 
DHr1 <- 166:176

KRr2 <- 289:302  
DHr2 <- 295:308

KRr3 <- 105:120 
DHr3 <- 113:127 

KRsvm <- 190:195 
DHsvm <- 196:201

KRr4 <- 256:266  
DHr4 <- 261:272

r5 <- 1:18
r6 <- 40:50

r1adjKR <- 230:238
r1adjDH <- 236:244
  
r1.2adjDH <- 245:266
r1.2adjKR <- 239:260

r4adKR <- 399:418
r4adDH <- 404:423

r4.2adKR <- 261:271
r4.2adDH <- 267:277

r4.3KR <- 179:184
r4.3DH <- 185:190

r4.3KR <- 104:108
r4.3DH <- 112:116

helix1_KR <- 80:96
helix1_DH <- 88:104
  
helix2_KR <- 341:358
helix2_DH <- 347:364

sageAKR <- 111:116
sageBKR <- 185:191
sageCKR <- 287:293

sageADH <- 111:116
sageBDH <- 185:191
sageCDH <- 287:293

total1 <- 1:20

target <- DH_out
region <- total1




out <- c()
for (i in 1:length(target)) {
  i2 <- target[[i]][, 1]
  i2[is.na(i2)] <- "-"
  out <- c(out, paste(i2, collapse = ""))
}

out2 <- c()
start_position <-min(region) 
end_position <- max(region)
for (i in 1:length(out)) {
  out2 <- c(out2, substr(out[i], start_position, end_position))
}

logo <- ggseqlogo(out2, seq_type = 'aa', method = 'prob')
logo

filename <- paste0("/home/q31032mw/Dropbox (The University of Manchester)/Max/17_ML_Project/Figures scripts/", "20240327_","KRvsDH_DH_NysC_Mod.2","_r4c_",".jpeg")
ggsave(filename, logo, device = "jpeg", width = 10, height = 4)


# loop for all in pdb
KS_length <- 424
ranges_list <- c()
for (i in seq(1, KS_length, by = 20)) {
  ranges_list <- c(ranges_list, i)
}

for(iii in ranges_list){
  target <- KR_out
  print(iii)
  out <- c()
  for (i in 1:length(target)) {
    i2 <- target[[i]][, 1]
    i2[is.na(i2)] <- "-"
    out <- c(out, paste(i2, collapse = ""))
  }
  
  out2 <- c()
  start_position <-iii 
  end_position <- iii + 19
  print(end_position)
  for (ii in 1:length(out)) {
    out2 <- c(out2, substr(out[ii], start_position, end_position))
  }
  
  logo <- ggseqlogo(out2, seq_type = 'aa', method = 'prob')
  filename <- paste0("/home/q31032mw/Dropbox (The University of Manchester)/Max/17_ML_Project/Figures scripts/", "20240410","_NysC_Mod.2_KRalign_", iii ,".jpeg")
  ggsave(filename, logo, device = "jpeg", width = 10, height = 4)
}





##################################
# True .pdb b-factor colouration #
##################################
for (i in 1:length(mat3_data[DH_test_set])) {
  protname <- names(mat3_data[i])
  name_in_list <- paste0(route, protname)
  pdb <- read.pdb(name_in_list)
  
  m3 <- mat3_data[[i]][,5]
  m3 <- as.character(m3)  # Convert to character string
  m3_numbers <- as.numeric(unlist(strsplit(m3, split = " ")))
  m3_numbers <- m3_numbers[1:(length(m3_numbers)/2)]
  
  
  frequency_table <- table(pdb$atom$resno)  # Assuming pdb is a data frame
  b_values <- c()
  for(ii in 1:length(frequency_table)){
    f <- frequency_table[ii]
    b_values <- c(b_values, rep(m3_numbers[ii], times = (f/2)))
  }
  pdb$atom$b <- b_values
  write.pdb(pdb = pdb, file = paste0(output_folder, paste0("DH",date, experiment_type , protname)) )
}


##################################


