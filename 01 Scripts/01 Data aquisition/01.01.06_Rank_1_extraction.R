### Top structure extraction ###
## 

# Inputs # 
#Directory <- 'C:\\Users\\q31032mw\\Dropbox (The University of Manchester)\\Max\\17_ML_Project\\Test_fa_output'
Destination <- 'C:\\Users\\q31032mw\\Dropbox (The University of Manchester)\\Max\\17_ML_Project\\Rank_1\\20231211_KS_dimers'
Directory <- 'D:\\03 HPC downloads\\20231206_KS_extras_dimers_output'

# Move files #
filenames <- list.files(Directory)
filenames <- filenames[grepl("unrelaxed_rank_001_", x = filenames)]

for (i in filenames){
  file.copy(from = paste0(Directory,"\\", i), to = paste0(Destination))
  
}


# Rename files # 
filenames <- list.files(Destination)
for (i in filenames) {
  file.rename(from = paste0(Destination ,"\\" ,i), 
              #to = paste0(Destination, "\\", gsub(pattern = "_unrelaxed_rank_001_alphafold2_ptm_model_._seed_000", replacement = "", x = i) )
              to = paste0(Destination, "\\", gsub(pattern = "_dimer_unrelaxed_rank_001_alphafold2_multimer_v3_model_._seed_000", replacement = "", x = i) )
  )

              }

