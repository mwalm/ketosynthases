#### Script for AT domains ####

### B-score graph ###
# Install and load the msa package
library(msa)
library(Biostrings)
library(bio3d)

## amino acid sequences import ## 
route <-  "/home/q31032mw/Dropbox (The University of Manchester)/Max/17_ML_Project/Rank_1/Aligned_big_set3_SNAKELESS/"

# Extract file names in route # 
files_list <- list.files(path = route)


# NetworkX coordinates ## 
netX_coords <- read.csv("/home/q31032mw/Dropbox (The University of Manchester)/Max/17_ML_Project/Figures scripts/20240201_AT_graphViz.csv")



# Create data structures for B-score and residue #
total_data <- list()
mat1_data <- list()
for (i in files_list){
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
  colnames(mat3_data[[i]])[ncol(mat3_data[[i]]) - 1] <- "combined_nodes"
}



Substrates <-  read.csv("/home/q31032mw/Dropbox (The University of Manchester)/Max/17_ML_Project/Webscraping/ClusterCAD Webscraping/20230707_AT_SMILES_unified.csv")
mmal_index <- which(Substrates$AT_substrate == "mmal")
mal_index <- which(Substrates$AT_substrate == "mal")
mmal_names <- Substrates$Gene[mmal_index]
mal_names <- Substrates$Gene[mal_index]

# Define your reference sequence and other sequences
reference_sequence_mmal2 <- read.pdb(paste0(route, "/eryAI_Mod.2.pdb_aligned.pdb")) 
reference_sequence_mmal <- pdbseq(reference_sequence_mmal2)

reference_sequence_mal2 <- read.pdb(paste0(route, "/pikAI_Mod.3.pdb_aligned.pdb")) 
reference_sequence_mal <- pdbseq(reference_sequence_mal2)


## B_score indexing ##
mmal_names_index <- c()
for (i in mmal_names) {
  print(i)
  index <- grep(x = names(mat1_data), pattern = i)
  if (length(index) != 0) {
    mmal_names_index <- c(mmal_names_index, index)
  }
}

mal_names_index <- c()
for (i in mal_names) {
  print(i)
  index <- grep(x = names(mat1_data), pattern = i)
  if (length(index) != 0) {
    mal_names_index <- c(mal_names_index, index)
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
# Mal #
{
countr <- 0
namez <- names(mat1_data[mal_names_index])
mat2_data_mal <- list()
for (i in mat1_data[mal_names_index]) {
  alignment <- pairwiseAlignment(
    subject = paste(reference_sequence_mal, collapse = ""),
    pattern = paste(i[,1], collapse = ''),
    type = "local-global"
  )
  alignment2 <- pairwiseAlignment(
    pattern = paste(reference_sequence_mal, collapse = ""),
    subject = paste(i[,1], collapse = ''),
    type = "global-local"
  )
  countr <- countr + 1
  name <- namez[countr]
  print(name)
  
  input_mat <- insert_rows(
    i, 
    which(as.matrix(alignment) == "-"), 
    length(reference_sequence_mal)
  )
  mat2_data_mal[[name]] <- input_mat
  
  if(length(grep(x = as.matrix(alignment2), pattern = "-")) == 0){next}
  
  else{mat2_data_mal[[name]] <- mat2_data_mal[[name]][-grep(x = as.matrix(alignment2), pattern = "-"),]}
  
  if(nrow(mat2_data_mal[[name]]) >= length(reference_sequence_mal)){
    mat2_data_mal[[name]] <- mat2_data_mal[[name]][1:length(reference_sequence_mal),]
    }
}
}



# Mmal #
{
countr <- 0
namez <- names(mat1_data[mmal_names_index])
mat2_data_mmal <- list()
for (i in mat1_data[mmal_names_index]) {
  alignment <- pairwiseAlignment(
    subject = paste(reference_sequence_mmal, collapse = ""),
    pattern = paste(i[,1], collapse = ''),
    type = "local-global"
  )
  alignment2 <- pairwiseAlignment(
    pattern = paste(reference_sequence_mmal, collapse = ""),
    subject = paste(i[,1], collapse = ''),
    type = "global-local"
  )
  countr <- countr + 1
  name <- namez[countr]
  print(name)
  
  input_mat <- insert_rows(
    i, 
    which(as.matrix(alignment) == "-"), 
    length(reference_sequence_mmal)
  )
  mat2_data_mmal[[name]] <- input_mat
  
  if(length(grep(x = as.matrix(alignment2), pattern = "-")) == 0){next}
  
  else{mat2_data_mmal[[name]] <- mat2_data_mmal[[name]][-grep(x = as.matrix(alignment2), pattern = "-"),]}
  
  if(nrow(mat2_data_mmal[[name]]) >= length(reference_sequence_mmal)){
    mat2_data_mmal[[name]] <- mat2_data_mmal[[name]][1:length(reference_sequence_mmal),]
  }
}
}

## B-score averages ## 
{
# Mmal #
mmal_b_scores <- matrix(rep(1:nrow(mat2_data_mmal[[1]])))
for (i in mat2_data_mmal) {
  mmal_b_scores <- cbind(mmal_b_scores, i[,2])
}
mmal_b_scores <- mmal_b_scores[,-1]
# Coerce each element to numeric
mmal_b_scores_numeric <- apply(mmal_b_scores, 2, as.numeric)

# Calculate row means while ignoring NA values
mmal_b_average <- rowMeans(mmal_b_scores_numeric, na.rm = TRUE)


# Mal #
mal_b_scores <- matrix(rep(1:nrow(mat2_data_mal[[1]])))
for (i in mat2_data_mal) {
  mal_b_scores <- cbind(mal_b_scores, i[1:nrow(mat2_data_mal[[1]]),2])
}
mal_b_scores <- mal_b_scores[,-1]
# Coerce each element to numeric
mal_b_scores_numeric <- apply(mal_b_scores, 2, as.numeric)

# Calculate row means while ignoring NA values
mal_b_average <- rowMeans(mal_b_scores_numeric, na.rm = TRUE)


# Apply score to atomic map #
atom_map_mal <- reference_sequence_mal2[["atom"]][["resno"]]
atom_matrix_mal <- matrix(nrow = length(atom_map_mal), ncol = 2)
atom_matrix_mal[,1] <- atom_map_mal
atom_matrix_mal[,2] <- mal_b_average[atom_matrix_mal[,1]]
output_b_score_mal <- reference_sequence_mal2
output_b_score_mal[["atom"]][["b"]] <- atom_matrix_mal[,2]
write.pdb(pdb = output_b_score_mal, file = paste0(route, "/mal_average_b_score_pikAI_mod.3.pdb") )

atom_map_mmal <- reference_sequence_mmal2[["atom"]][["resno"]]
atom_matrix_mmal <- matrix(nrow = length(atom_map_mmal), ncol = 2)
atom_matrix_mmal[,1] <- atom_map_mmal
atom_matrix_mmal[,2] <- mmal_b_average[atom_matrix_mmal[,1]]
output_b_score_mmal <- reference_sequence_mmal2
output_b_score_mmal[["atom"]][["b"]] <- atom_matrix_mmal[,2]
write.pdb(pdb = output_b_score_mmal, file = paste0(route, "/mmal_average_b_score_eryAI_mod.2.pdb") )
}


## NetX averages ##
## Perform alignment, add and remove rows from query sequences ##
mmal_names_index <- c()
for (i in mmal_names) {
  print(i)
  index <- grep(x = names(mat3_data), pattern = i)
  if (length(index) != 0) {
    mmal_names_index <- c(mmal_names_index, index)
  }
}

mal_names_index <- c()
for (i in mal_names) {
  print(i)
  index <- grep(x = names(mat3_data), pattern = i)
  if (length(index) != 0) {
    mal_names_index <- c(mal_names_index, index)
  }
}

# Mal #
{
  countr <- 0
  namez <- names(mat3_data[mal_names_index])
  mat2_data_mal <- list()
  for (i in mat3_data[mal_names_index]) {
    alignment <- pairwiseAlignment(
      subject = paste(reference_sequence_mal, collapse = ""),
      pattern = paste(i[,1], collapse = ''),
      type = "local-global"
    )
    alignment2 <- pairwiseAlignment(
      pattern = paste(reference_sequence_mal, collapse = ""),
      subject = paste(i[,1], collapse = ''),
      type = "global-local"
    )
    countr <- countr + 1
    name <- namez[countr]
    print(name)
    
    input_mat <- insert_rows(
      i, 
      which(as.matrix(alignment) == "-"), 
      length(reference_sequence_mal)
    )
    mat2_data_mal[[name]] <- input_mat
    
    if(length(grep(x = as.matrix(alignment2), pattern = "-")) == 0){next}
    
    else{mat2_data_mal[[name]] <- mat2_data_mal[[name]][-grep(x = as.matrix(alignment2), pattern = "-"),]}
    
    if(nrow(mat2_data_mal[[name]]) >= length(reference_sequence_mal)){
      mat2_data_mal[[name]] <- mat2_data_mal[[name]][1:length(reference_sequence_mal),]
    }
  }
}

# Mmal #
{
  countr <- 0
  namez <- names(mat3_data[mmal_names_index])
  mat2_data_mmal <- list()
  for (i in mat3_data[mmal_names_index]) {
    alignment <- pairwiseAlignment(
      subject = paste(reference_sequence_mmal, collapse = ""),
      pattern = paste(i[,1], collapse = ''),
      type = "local-global"
    )
    alignment2 <- pairwiseAlignment(
      pattern = paste(reference_sequence_mmal, collapse = ""),
      subject = paste(i[,1], collapse = ''),
      type = "global-local"
    )
    countr <- countr + 1
    name <- namez[countr]
    print(name)
    
    input_mat <- insert_rows(
      i, 
      which(as.matrix(alignment) == "-"), 
      length(reference_sequence_mmal)
    )
    mat2_data_mmal[[name]] <- input_mat
    
    if(length(grep(x = as.matrix(alignment2), pattern = "-")) == 0){next}
    
    else{mat2_data_mmal[[name]] <- mat2_data_mmal[[name]][-grep(x = as.matrix(alignment2), pattern = "-"),]}
    
    if(nrow(mat2_data_mmal[[name]]) >= length(reference_sequence_mmal)){
      mat2_data_mmal[[name]] <- mat2_data_mmal[[name]][1:length(reference_sequence_mmal),]
    }
  }
}

## NetX averages ## 
{
  # Mmal #
  mmal_netx_scores <- matrix(rep(1:nrow(mat2_data_mmal[[1]])))
  for (i in mat2_data_mmal) {
    mmal_netx_scores <- cbind(mmal_netx_scores, i[,5])
  }
  mmal_netx_scores <- mmal_netx_scores[,-1]
  # Coerce each element to numeric
  mmal_netx_scores_numeric <- apply(mmal_netx_scores, 2, as.numeric)
  
  # Calculate row means while ignoring NA values
  mmal_netx_average <- rowMeans(mmal_netx_scores_numeric, na.rm = TRUE)
  
  
  # Mal #
  mal_netx_scores <- matrix(rep(1:nrow(mat2_data_mal[[1]])))
  for (i in mat2_data_mal) {
    mal_netx_scores <- cbind(mal_netx_scores, i[,5])
  }
  mal_netx_scores <- mal_netx_scores[,-1]
  # Coerce each element to numeric
  mal_netx_scores_numeric <- apply(mal_netx_scores, 2, as.numeric)
  
  # Calculate row means while ignoring NA values
  mal_netx_average <- rowMeans(mal_netx_scores_numeric, na.rm = TRUE)
  
  
  # Apply score to atomic map #
  atom_map_mal <- reference_sequence_mal2[["atom"]][["resno"]]
  atom_matrix_mal <- matrix(nrow = length(atom_map_mal), ncol = 2)
  atom_matrix_mal[,1] <- atom_map_mal
  atom_matrix_mal[,2] <- mal_netx_average[atom_matrix_mal[,1]]
  output_netx_score_mal <- reference_sequence_mal2
  output_netx_score_mal[["atom"]][["b"]] <- atom_matrix_mal[,2]
  write.pdb(pdb = output_netx_score_mal, file = paste0(route, "/mal_average_netx_score_pikAI_mod.3.pdb") )
  
  atom_map_mmal <- reference_sequence_mmal2[["atom"]][["resno"]]
  atom_matrix_mmal <- matrix(nrow = length(atom_map_mmal), ncol = 2)
  atom_matrix_mmal[,1] <- atom_map_mmal
  atom_matrix_mmal[,2] <- mmal_netx_average[atom_matrix_mmal[,1]]
  output_netx_score_mmal <- reference_sequence_mmal2
  output_netx_score_mmal[["atom"]][["b"]] <- atom_matrix_mmal[,2]
  write.pdb(pdb = output_netx_score_mmal, file = paste0(route, "/mmal_average_netx_score_eryAI_mod.2.pdb") )
}















### Logo Diagrams ###
library(ggplot2)
library(ggseqlogo)

# list of residues from cross words # 
r1 <- 1:13
r2 <- 85:100
r3 <- 190:210

## Mal ## 
out <- c()
for (i in 1:length(mat2_data_mal)) {
  i2 <- mat2_data_mal[[i]][, 1]
  i2[is.na(i2)] <- "-"
  out <- c(out, paste(i2, collapse = ""))
}

out2 <- c()
start_position <-min(r3) 
end_position <- max(r3)
for (i in 1:length(out)) {
  out2 <- c(out2, substr(out[i], start_position, end_position))
}

logo <- ggseqlogo(out2, seq_type = 'aa', method = 'prob')
logo

filename <- paste0("/home/q31032mw/Dropbox (The University of Manchester)/Max/17_ML_Project/Figures scripts/","20240320_","mal_","190-210",".jpeg")
ggsave(filename, logo, device = "jpeg", width = 8, height = 6)


## Mmal ## 
r1 <- 1:13
r2 <- 100:116
r3 <- 185:197

out <- c()
for (i in 1:length(mat2_data_mmal)) {
  i2 <- mat2_data_mmal[[i]][, 1]
  i2[is.na(i2)] <- "-"
  out <- c(out, paste(i2, collapse = ""))
}

out2 <- c()
start_position <-min(r3) 
end_position <- max(r3)
for (i in 1:length(out)) {
  out2 <- c(out2, substr(out[i], start_position, end_position))
}

logo <- ggseqlogo(out2, seq_type = 'aa', method = 'prob')
logo

filename <- paste0("/home/q31032mw/Dropbox (The University of Manchester)/Max/17_ML_Project/Figures scripts/","20240320_","mmal_","185-197",".jpeg")
ggsave(filename, logo, device = "jpeg", width = 8, height = 6)
