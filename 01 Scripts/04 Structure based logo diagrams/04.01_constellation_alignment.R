### Protein Constellation viewer script ###

## Dependencies ##
library(bio3d)
library(ggseqlogo)
library(ggplot2)

## Inputs ##
# define alignment tolerance range on alpha carbons #
angstrom_range <- 2.5 # units = Angstroms

# set reference proteins .pdb #
ref_prot_NR <- "C:\\Users\\q31032mw\\The University of Manchester Dropbox\\Max Walmsley\\Max\\17_ML_Project\\Rank_1\\20240911_KS_dimers_extras_aligned_reviewed\\AceP4_Mod.6.pdb_aligned.pdb"
ref_prot_KR <- "C:\\Users\\q31032mw\\The University of Manchester Dropbox\\Max Walmsley\\Max\\17_ML_Project\\Rank_1\\20240911_KS_dimers_extras_aligned_reviewed\\SeEryAI_Mod.2.pdb_aligned.pdb"
ref_prot_DH <- "C:\\Users\\q31032mw\\The University of Manchester Dropbox\\Max Walmsley\\Max\\17_ML_Project\\Rank_1\\20240911_KS_dimers_extras_aligned_reviewed\\nysC_Mod.2.pdb_aligned.pdb"
ref_prot_ER <- "C:\\Users\\q31032mw\\The University of Manchester Dropbox\\Max Walmsley\\Max\\17_ML_Project\\Rank_1\\20240911_KS_dimers_extras_aligned_reviewed\\nysC_Mod.4.pdb_aligned.pdb"

ref_prot_NR <- read.pdb(ref_prot_NR)
ref_prot_KR <- read.pdb(ref_prot_KR)
ref_prot_DH <- read.pdb(ref_prot_DH)
ref_prot_ER <- read.pdb(ref_prot_ER)


# Set a directory of pre-aligned structures #
folder <- "C:\\Users\\q31032mw\\The University of Manchester Dropbox\\Max Walmsley\\Max\\17_ML_Project\\Rank_1\\20240911_KS_dimers_extras_aligned_reviewed\\"
directory <- list.files(folder, full.names = TRUE)

# Reduction type #
type_matrix <- read.csv("C:\\Users\\q31032mw\\The University of Manchester Dropbox\\Max Walmsley\\Max\\17_ML_Project\\Rank_1\\20231108_KS_expandedset_reduction_state_labels_new_borders.csv")


### Map directory to type ###
# protein type names #
NR_prots_names <- type_matrix[which(type_matrix$NR == 1), 1]
KR_prots_names <- type_matrix[which(type_matrix$KR == 1), 1]
DH_prots_names <- type_matrix[which(type_matrix$DH == 1), 1]
ER_prots_names <- type_matrix[which(type_matrix$ER == 1), 1]

## file paths for proteins by type ##
{
# NR # 
NR_dir <- c()
for (i in NR_prots_names) {
  NR_dir <- c(NR_dir, grep(pattern = paste0(i,".pdb"), directory))
}
NR_files <- directory[NR_dir] 

NR_pdb <- list()
countr <- 0L
for (i in NR_files) {
  countr <- countr + 1L
  print(paste0(countr, " of ", length(NR_files)))
  pdb <- read.pdb(i)
  pdb_CA <- which(pdb$atom$elety == "CA")
  prot_name <- gsub(pattern = paste0(folder), replacement = "", i, fixed = TRUE)
  prot_name <- gsub(pattern = ".pdb_aligned.pdb", replacement = "", prot_name, fixed = TRUE)
  NR_pdb[[prot_name]] <- pdb$atom[pdb_CA,]
}

# KR # 
KR_dir <- c()
for (i in KR_prots_names) {
  KR_dir <- c(KR_dir, grep(pattern = paste0(i,".pdb"), directory))
}
KR_files <- directory[KR_dir] 

KR_pdb <- list()
countr <- 0L
for (i in KR_files) {
  countr <- countr + 1L
  print(paste0(countr, " of ", length(KR_files)))
  pdb <- read.pdb(i)
  pdb_CA <- which(pdb$atom$elety == "CA")
  prot_name <- gsub(pattern = paste0(folder), replacement = "", i, fixed = TRUE)
  prot_name <- gsub(pattern = ".pdb_aligned.pdb", replacement = "", prot_name, fixed = TRUE)
  KR_pdb[[prot_name]] <- pdb$atom[pdb_CA,]
}

# DH # 
DH_dir <- c()
for (i in DH_prots_names) {
  DH_dir <- c(DH_dir, grep(pattern = paste0(i,".pdb"), directory))
}
DH_files <- directory[DH_dir] 

DH_pdb <- list()
countr <- 0L
for (i in DH_files) {
  countr <- countr + 1L
  print(paste0(countr, " of ", length(DH_files)))
  pdb <- read.pdb(i)
  pdb_CA <- which(pdb$atom$elety == "CA")
  prot_name <- gsub(pattern = paste0(folder), replacement = "", i, fixed = TRUE)
  prot_name <- gsub(pattern = ".pdb_aligned.pdb", replacement = "", prot_name, fixed = TRUE)
  DH_pdb[[prot_name]] <- pdb$atom[pdb_CA,]
}

# ER # 
ER_dir <- c()
for (i in ER_prots_names) {
  ER_dir <- c(ER_dir, grep(pattern = paste0(i,".pdb"), directory))
}
ER_files <- directory[ER_dir] 

ER_pdb <- list()
countr <- 0L
for (i in ER_files) {
  countr <- countr + 1L
  print(paste0(countr, " of ", length(ER_files)))
  pdb <- read.pdb(i)
  pdb_CA <- which(pdb$atom$elety == "CA")
  prot_name <- gsub(pattern = paste0(folder), replacement = "", i, fixed = TRUE)
  prot_name <- gsub(pattern = ".pdb_aligned.pdb", replacement = "", prot_name, fixed = TRUE)
  ER_pdb[[prot_name]] <- pdb$atom[pdb_CA,]
}

}

## reference tables ##
CA_NR_table <- ref_prot_NR$atom[which(ref_prot_NR$atom$elety == "CA"), c(5,7,9,10,11)]
CA_KR_table <- ref_prot_KR$atom[which(ref_prot_KR$atom$elety == "CA"), c(5,7,9,10,11)]
CA_DH_table <- ref_prot_DH$atom[which(ref_prot_DH$atom$elety == "CA"), c(5,7,9,10,11)]
CA_ER_table <- ref_prot_ER$atom[which(ref_prot_ER$atom$elety == "CA"), c(5,7,9,10,11)]

# Amino acid trip code #
acids <- c("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "None")

# Add columns # 
add_amino_acid_columns <- function(table) {
  for (acid in acids) {
    table[[acid]] <- 0  
  }
  return(table)
}

CA_NR_table <- add_amino_acid_columns(CA_NR_table)
CA_KR_table <- add_amino_acid_columns(CA_KR_table)
CA_DH_table <- add_amino_acid_columns(CA_DH_table)
CA_ER_table <- add_amino_acid_columns(CA_ER_table)


## Identify nearest residues ## 
INR <- function(CA_table, type_pdb){
  for(i in 1:nrow(CA_table)){
    print(paste0(i, " of ", nrow(CA_table)))
    residue <- CA_table[i,]
    
    residue_x <- c((residue$x - angstrom_range), (residue$x + angstrom_range))
    residue_y <- c((residue$y - angstrom_range), (residue$y + angstrom_range))
    residue_z <- c((residue$z - angstrom_range), (residue$z + angstrom_range))
    
    for(ii in 1:length(type_pdb)){
      # x coords # 
      x_range <- which(type_pdb[[ii]]$x > residue_x[1] & type_pdb[[ii]]$x < residue_x[2])
      
      # y coords #
      y_range <- which(type_pdb[[ii]]$y > residue_y[1] & type_pdb[[ii]]$y < residue_y[2])
      
      # z coords # 
      z_range <- which(type_pdb[[ii]]$z > residue_z[1] & type_pdb[[ii]]$z < residue_z[2])
      
      range_xyz <- intersect(intersect(x_range, y_range), z_range)

      if(length(range_xyz) == 0){
        CA_table[i, 26] <- CA_table[i, 26] + 1 
        next
        }
      
      if(length(range_xyz) > 1){
        proximal_res <- c()
        
        for(iiii in range_xyz){
          diff_x <- residue$x - type_pdb[[ii]][iiii, 9]
          diff_y <- residue$y - type_pdb[[ii]][iiii, 10]
          diff_z <- residue$z - type_pdb[[ii]][iiii, 11]
          
          proximal_res <- c(proximal_res, sum( (((diff_x)^2)^0.5), (((diff_y)^2)^0.5), (((diff_z)^2)^0.5)) )
        }
        
      min_proximal_res <- min(proximal_res)
      range_xyz <- range_xyz[which(proximal_res == min_proximal_res)]
      if(length(range_xyz) > 1){print("fuck!")
        print(proximal_res)
        print(range_xyz)}
      }
      aa <- type_pdb[[ii]][range_xyz,5]
      aa_col <- grep(pattern = aa, x = acids) + 5
      CA_table[i, aa_col] <- CA_table[i, aa_col] + 1
    }
  }
  return(CA_table)
}

CA_NR_table <- INR(CA_NR_table, NR_pdb)
CA_KR_table <- INR(CA_KR_table, KR_pdb)
CA_DH_table <- INR(CA_DH_table, DH_pdb)
CA_ER_table <- INR(CA_ER_table, ER_pdb)


## CA_table to % residue ## 
CA_percentage_table <- function(CA_table){
  for(i in 1:nrow(CA_table)){
    total_res <- sum(CA_table[i, 6:26])
    for(ii in 1:length(acids)){
      if(CA_table[i,(ii + 5)] == 0){next}
      CA_table[i,(ii + 5)] <- (CA_table[i,(ii + 5)]/total_res)*100
      
    }
  }
  colnames(CA_table) <- c("resid", "resno", "xc", "y", "z", "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L","K", "M", "F", "P", "S", "T", "W", "Y", "V", "X")
  return(CA_table)
}

CA_NR_table_per <- CA_percentage_table(CA_NR_table)
CA_KR_table_per <- CA_percentage_table(CA_KR_table)
CA_DH_table_per <- CA_percentage_table(CA_DH_table)
CA_ER_table_per <- CA_percentage_table(CA_ER_table)


## create arbitrary sequences for logo diagrams ##
resolution <- 0 # Decimal places on percentage. e.g 83.33%  ==> 83 if resolution = 0, 83.3 if resolution = 1 # 

seq_table <- function(CA_table_per, resolution){
  # Blank sheet #
  seq_length <- nrow(CA_table_per)
  na_table <- matrix(NA, nrow = 100*(10^resolution), ncol = seq_length)
  
  # populate with imaginary sequences based on probability score #
  for(i in 1:ncol(na_table)){
    
    residues <- c()
    for(ii in 6:26){
      residues <- c(residues, rep(colnames(CA_table_per)[ii], (round(CA_table_per[i,ii], digits = resolution))))

    }

    if(length(residues) < nrow(na_table)){residues <- c(residues, rep("X", (nrow(na_table)-length(residues))))}
    for(iii in 1:nrow(na_table)){
     na_table[iii,i] <- residues[iii]
    }
  }
  
  return(na_table)
}

NR_seqs_table <- seq_table(CA_NR_table_per, resolution)
KR_seqs_table <- seq_table(CA_KR_table_per, resolution)
DH_seqs_table <- seq_table(CA_DH_table_per, resolution)
ER_seqs_table <- seq_table(CA_ER_table_per, resolution)

## create string of sequences ##
seqvec <- function(seqs_table){
  seqvec <- c()
  for(i in 1:nrow(seqs_table)){
    seqvec <- c(seqvec, paste0(seqs_table[1,], collapse = ""))
  }
  return(seqvec)
}

NR_seqs <- seqvec(NR_seqs_table)
KR_seqs <- seqvec(KR_seqs_table)
DH_seqs <- seqvec(DH_seqs_table)
ER_seqs <- seqvec(ER_seqs_table)


out_dir <- "C:\\Users\\q31032mw\\The University of Manchester Dropbox\\Max Walmsley\\Max\\17_ML_Project\\Figures scripts\\20240926_constellation_maps\\"
date <- "20241001"


logos <- function(type,type_seqs_table, CA_type_table_per){
  ref_length <- nrow(CA_type_table_per)
  
  ranges_list <- c()
  
  for (i in seq(1, ref_length, by = 20)) {
    ranges_list <- c(ranges_list, i)
  }
  
  for(ii in ranges_list){
    start_position <- ii
    end_position <- ii + 19
    if(end_position > nrow(CA_type_table_per)){end_position <- nrow(CA_type_table_per)}
    seq_vector <- c()
    for(iii in 1:nrow(type_seqs_table)){
      seq_vector <- c(seq_vector, paste0(type_seqs_table[iii , (start_position:end_position)], collapse = "")) 
    }
    logo <- ggseqlogo(seq_vector, seq_type = 'aa', method = 'bits')
    filename <- paste0(out_dir, date, "_", type, "_", ii, "_bits" ,".jpeg")
    ggsave(filename, logo, device = "jpeg", width = 10, height = 4)
  }
}
logos("NR", NR_seqs_table, CA_NR_table_per)
logos("KR", KR_seqs_table, CA_KR_table_per)
logos("DH", DH_seqs_table, CA_DH_table_per)
logos("ER", ER_seqs_table, CA_ER_table_per)

  

