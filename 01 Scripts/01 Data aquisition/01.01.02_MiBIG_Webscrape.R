###MiBIG scraping tool and module information extraction tool###
##Function of this script is to download .gbK (Genbank) files from MiBIG##
##PKS domains within .gbks are ordered into modules.##

##Dependencies:##
library(rvest)
library(polite)
library(genbankr)
library(stringr)
library(Biostrings)

##inputs:##
#date is used to label outputs#
date <- "yymmdd"

#df3 is a list of BGCs to be worked on.#
df3 <- read.csv("C:\\Users\\q31032mw\\Dropbox (The University of Manchester)\\Max\\17_ML_Project\\GitHub\\01 Scripts\\01 Data aquisition\\01 Webscraping\\df3.csv", row.names = 1)

mibig <- "https://mibig.secondarymetabolites.org/repository/"
mibig1 <- bow("https://mibig.secondarymetabolites.org/repository/")

#script for downloading .gbk files from MiBIG, then extracting module information.#
framesize <- c()
for(bgc_name in df3[,1])
  try(
    {
  #get antiSMASH documentation .gbk files#
  dir.create(paste0(date, "_genbank file"))
  download.file(url = paste("https://mibig.secondarymetabolites.org/repository/", 
                            bgc_name, "/", 
                            bgc_name, 
                            ".gbk", 
                            sep = ""), 
                destfile = paste(date, "_genbank file/", bgc_name, ".gbk", sep = ""))
  
  
  #parse .gbk file to a dataframe#
  x <- parseGenBank(paste(date, "_genbank file/", bgc_name, ".gbk", sep = ""))
  
  #Identify features containing NRPS/PKS features#
  x2 <- grep("NRPS_PKS", x[["FEATURES"]])
  
  #Retrieve gene names for PKS genes. This will also pick up Trans-thioesterases and NRPS/PKS associated genes#
  #This loop is a recurring break point. If dfy2 has "features ##" for gene names, the below loop needs to be modified to include 
  x4 <- c()
  for(x3 in x2){
    genename <- x[["FEATURES"]][[x3]][["gene"]]
    if(is.null(x[["FEATURES"]][[x3]][["gene"]]) == TRUE){genename <- x[["FEATURES"]][[x3]][["locus_tag"]]}
    if(is.null(x[["FEATURES"]][[x3]][["gene"]]) == TRUE){genename <- x[["FEATURES"]][[x3]][["product"]]}
    if(is.null(genename) == TRUE){genename <- paste("Feature ID_", x3, sep = "")}
    x4<- c(x4, genename)
    
  }
  
  # Function to add unique numbers
  add_unique_numbers <- function(x) {
    counts <- table(x)
    for (item in unique(x)) {
      if (counts[item] > 1) {
        indices <- which(x == item)
        for (i in 1:length(indices)) {
          x[indices[i]] <- paste0(item, i)
        }
      }
    }
    return(x)
  }
  
  # Apply the function to your vector
  x4 <- add_unique_numbers(x4)

  
  
  #Create dfx, a dataframe for holding domain information#
  dfx <- data.frame(row.names = c(1:length(x4)))
  
  dfx$Gene <- x4
  dfx$FramePosition <- x2
  
  #Expand dfx to include parsed .gbk coordinates for PKS domains#
  y2 <- c()
  for(y1 in x2){
    y2<- c(y2, grep("Domain", x[["FEATURES"]][[y1]]))
    countr <- 0
    y2 <- length(y2)
    for (y3 in 1:y2) {
      countr = countr + 1
      y4 <- paste("Domain", ".", countr, sep="")
      dfx[,y4] <- NA
    }  
  }
  dfx
  y2 <- c()
  y3 <- c()
  countr = 0
  for (z in x2) {
    y2<- grep("Domain", x[["FEATURES"]][[z]])
    print(y2)
    y3 <- c(y3,y2)
    if(length(y2) <= 3){y2 = c(y2, NA, NA, NA) } else {y2==y2}
    print(y2)
    countr = countr + 1
    leng <- length(y2) + 2
    dfx[countr, 3:leng] <- y2
  }
  
  #Create dfy, a dataframe for holding amino acid sequences of domains#
  Domain_types <-  c('KS', 'AT', 'KR', 'DH', 'ER', 'MT', 'ACP', 'TE') 
  Docks <- c('N.term.dock', 'C.term.dock' )
  module_limit <- 0:9
  df_headers <- c()
  for (i in module_limit) {
    df_headers <- c(df_headers, paste(Domain_types,".", i, sep = "" ))
  }
  dfy <- data.frame(row.names = dfx$Gene)
  dfy[,Docks] <- NA
  for (i in df_headers) {
    dfy[,i] <- NA
  }
  
  #populate dfy
  #create dfz, a dataframe containing the order and names of domains in genes#
  dfz <- dfx
  for (i in 1:nrow(dfx)) 
  {
    for (ii in 3:ncol(dfx)) 
    {
      if (is.na(dfx[i,ii])==TRUE){next}
      print(dfx[i,ii])
      p1 <- scan(text = x[["FEATURES"]][[dfx[i,2]]][[dfx[i,ii]]], what = "")
      p2 <- sub(".*nrpspksdomains_","", p1[length(p1)])
      dfz[i,ii] <- p2
    }
  }
  
  #create dfz2, a dataframe with the same structure as dfz, where domain names are replaced with amino acid coordinate ranges.#
  
  dfz2 <- dfz
  for (i in 1:nrow(dfx)) {
    for (ii in 3:ncol(dfx)) 
    {
      if (is.na(dfx[i,ii])==TRUE){next}
      print(dfx[i,ii])
      p1 <- scan(text = x[["FEATURES"]][[dfx[i,2]]][[dfx[i,ii]]], what = "")
      p2 <- sub(".*nrpspksdomains_","", p1[3])
      p2 <- gsub("\\(|\\).", "", p2)
      dfz2[i,ii] <- p2
    }
  }
  
  #create domain_order dataframe. This is used as a template for renaming domains in the gbk file, based on their module#
  domain_order <- data.frame(dfx$Gene)
  test2 <- c()
  for (i in 1:nrow(dfz)){
    test <- !is.na(dfz[i,])
    test <- c(dfz[i,test])
    test <- sub(paste(dfz[i,1], "_", sep = ""),"" ,  test)
    test <- sub("PKS_","", test)
    #test <- sub("PP", "ACPx", test)
    test <- test[3:length(test)]
    test <- paste( unlist(test), collapse=' ')
    test2 <- c(test2, test)
  }
  domain_order$domainstring <- test2
  modulenumbers <- c()
  for(i in 1:nrow(dfz)){
    test <- as.list(strsplit(domain_order[i,2], split = "KS..."))
    for(ii in 1:length(test)){
      modulenumbers <- c(modulenumbers, paste("module.", (i-1), sep = ""))
    }
    for(iii in 1:length(modulenumbers)){
      domain_order[ ,iii + 2] <- NA
      colnames(domain_order)[iii+2] <- modulenumbers[iii]
    }
  }
  for (i in 1:nrow(dfz) ) {
    test <- as.list(strsplit(domain_order[i,2], split = "KS..."))
    countr <- 2
    for (ii in test) {
      for (iii in ii) {
        countr <- countr + 1 
        domain_order[i, countr] <- iii
        if(countr == length(test)){break}
      }
    }
  }
  
  #clean up domain names
  for(i in 1:nrow(domain_order)){
    for (ii in 3:ncol(domain_order)) {
      replacement_string <- c()
      
      if(is.na(domain_order[i,ii]) == TRUE){next}
      else{NULL}
      
      #order of operations is based on domain structure in the PKS. If the domains are not ordered in the gbk at follows, domains will be assigned the wrong amino acids.
      if(is.na(str_extract(string = domain_order[i,ii], pattern = "Docking_Nterm..")) == FALSE)
      {replacement_string <- c(replacement_string, str_extract(string = domain_order[i,ii], pattern = "Docking_Nterm..") )}
      if(is.na(str_extract(string = domain_order[i,ii], pattern = "AT..")) == FALSE)
      {replacement_string <- c(replacement_string, str_extract(string = domain_order[i,ii], pattern = "AT..") )}
      if(is.na(str_extract(string = domain_order[i,ii], pattern = "DH..")) == FALSE)
      {replacement_string <- c(replacement_string, str_extract(string = domain_order[i,ii], pattern = "DH..") )}
      if(is.na(str_extract(string = domain_order[i,ii], pattern = "ER..")) == FALSE)
      {replacement_string <- c(replacement_string, str_extract(string = domain_order[i,ii], pattern = "ER..") )}
      if(is.na(str_extract(string = domain_order[i,ii], pattern = "KR..")) == FALSE)
      {replacement_string <- c(replacement_string, str_extract(string = domain_order[i,ii], pattern = "KR..") )}
      if(is.na(str_extract(string = domain_order[i,ii], pattern = "PP..")) == FALSE)
      {replacement_string <- c(replacement_string, str_extract(string = domain_order[i,ii], pattern = "PP..") )}
      if(is.na(str_extract(string = domain_order[i,ii], pattern = "ACP..")) == FALSE)
      {replacement_string <- c(replacement_string, str_extract(string = domain_order[i,ii], pattern = "ACP..") )}
      if(is.na(str_extract(string = domain_order[i,ii], pattern = "Thioesterase..")) == FALSE)
      {replacement_string <- c(replacement_string, str_extract(string = domain_order[i,ii], pattern = "Thioesterase..") )}
      if(is.na(str_extract(string = domain_order[i,ii], pattern = "Docking_Cterm..")) == FALSE)
      {replacement_string <- noquote(c(replacement_string, str_extract(string = domain_order[i,ii], pattern = "Docking_Cterm..") ) )}
      
      if(is.null(replacement_string) == TRUE){replacement_string <- domain_order[i,ii]}
      
      domain_order[i,ii] <- paste(replacement_string, collapse = " ")
    }
  }
  
  
  #split the modules to domains, replace the domain numbers with the module numbers
  for(i in 1:nrow(domain_order)){
    for (ii in 3:ncol(domain_order)) {
      modnumb <- sub("module\\.", "", (colnames(domain_order[ii])))
      test <- as.list(domain_order[i,ii])
      test <- gsub("[.].", paste(".", replacement = modnumb, sep= ""), x = test )
      domain_order[i,ii] <- test
    }}
  print("WARNING: If the PKS module has more than 9 modules you are going to have a problem here.")
  print("This includes an N dock as a domain")
  print("This includes an N dock as a domain")
  #adds KS back to domain list
  for (i in 1:nrow(domain_order)) {
    for (ii in 4:ncol(domain_order)) 
    {
      if(grepl(pattern = "NA", domain_order[i,ii]) == TRUE){domain_order[i,ii] <- NA}
      if(is.na(domain_order[i,ii]) == FALSE){
        domain_order[i,ii] <- paste0("KS.", ii-3, " ", domain_order[i,ii])}
      else{NULL}
    }
  }
  
  #Rename NRPS_PKS domains using domain_order frame.#
    #Create DO_List, a character sting of modules in sequence order of genes in BGC#
    #Note: this is not in order of enzyme activity in the PKS#
  DO_list <- c()
  for (i in 1:nrow(domain_order)){
    for (ii in 3:ncol(domain_order)) 
    {
      if(grepl(pattern= "^$", domain_order[i,ii]) == TRUE){domain_order[i,ii] <- NA}
      for(iii in 1:length(strsplit(domain_order[i,ii], split = " ")[[1]])){
        if(is.na(domain_order[i,ii]) == TRUE){break}
        test <- (strsplit(domain_order[i,ii], split = " ")[[1]][[iii]])
        DO_list <- c(DO_list, test)
      }
    }
  }
  
  #create xed, a copy of x, the parsed .gbk file. 
    #rename the domain entries in xed with the module formatting in DO_list#
  xed <- x
  countr <- 1
  modnumb <- NA
  for (i in dfx$FramePosition) {
    modnumb <- sum(grepl("nrpspks", x = xed[["FEATURES"]][[i]])) + 5
    for(ii in 6:modnumb){
      print(xed[["FEATURES"]][[i]][[ii]])
      
      xed[["FEATURES"]][[i]][[ii]]<- paste0(xed[["FEATURES"]][[i]][[ii]],"_rID_", DO_list[countr])
      countr = countr + 1
    }  
  }
  
  #Loop populating dfy with amino acid coordinates of domains. 
  for (i in 1:nrow(dfx)) 
  {
    for (ii in 3:ncol(dfx)) 
    {
      
      if (grepl(pattern = "_ACP", x = dfz[i,ii]) == TRUE) {
        p1 <- scan(text = xed[["FEATURES"]][[dfx[i,2]]][[dfx[i,ii]]], what = "")
        p2 <- sub(".*nrpspksdomains_","", p1[3])
        p2 <- gsub("\\(|\\).", "", p2)
        p3 <- strsplit(p1[length(p1)], split= "\\.")
        dfycol <- paste("ACP", p3[[1]][[length(p3[[1]])]] ,sep = ".")
        if(nchar(dfycol) >= 5){dfycol <- substring(dfycol, 1, 5)} else{print("x")}
        dfy[i,dfycol] <- dfz2[i,ii]
      } 
      if (grepl(pattern = "_PP", x = dfz[i,ii]) == TRUE) {
        p1 <- scan(text = xed[["FEATURES"]][[dfx[i,2]]][[dfx[i,ii]]], what = "")
        p2 <- sub(".*nrpspksdomains_","", p1[3])
        p2 <- gsub("\\(|\\).", "", p2)
        p3 <- strsplit(p1[length(p1)], split= "\\.")
        dfycol <- paste("ACP", p3[[1]][[length(p3[[1]])]] ,sep = ".")
        if(nchar(dfycol) >= 5){dfycol <- substring(dfycol, 1, 5)} else{print("x")}
        dfy[i,dfycol] <- dfz2[i,ii]
      }
      if (grepl(pattern = "_KS", x = dfz[i,ii]) == TRUE) {
        p1 <- scan(text = xed[["FEATURES"]][[dfx[i,2]]][[dfx[i,ii]]], what = "")
        p2 <- sub(".*nrpspksdomains_","", p1[3])
        p2 <- gsub("\\(|\\).", "", p2)
        p3 <- strsplit(p1[length(p1)], split= "\\.")
        dfycol <- paste("KS", p3[[1]][[length(p3[[1]])]] ,sep = ".")
        if(nchar(dfycol) >= 5){dfycol <- substring(dfycol, 1, 4)} else{print("x")}
        dfy[i,dfycol] <- dfz2[i,ii]
      } 
      if (grepl(pattern = "_AT", x = dfz[i,ii]) == TRUE) {
        p1 <- scan(text = xed[["FEATURES"]][[dfx[i,2]]][[dfx[i,ii]]], what = "")
        p2 <- sub(".*nrpspksdomains_","", p1[3])
        p2 <- gsub("\\(|\\).", "", p2)
        p3 <- strsplit(p1[length(p1)], split= "\\.")
        dfycol <- paste("AT", p3[[1]][[length(p3[[1]])]] ,sep = ".")
        if(nchar(dfycol) >= 5){dfycol <- substring(dfycol, 1, 4)} else{print("x")}
        dfy[i,dfycol] <- dfz2[i,ii]
      }
      if (grepl(pattern = "_KR", x = dfz[i,ii]) == TRUE) {
        p1 <- scan(text = xed[["FEATURES"]][[dfx[i,2]]][[dfx[i,ii]]], what = "")
        p2 <- sub(".*nrpspksdomains_","", p1[3])
        p2 <- gsub("\\(|\\).", "", p2)
        p3 <- strsplit(p1[length(p1)], split= "\\.")
        dfycol <- paste("KR", p3[[1]][[length(p3[[1]])]] ,sep = ".")
        if(nchar(dfycol) >= 5){dfycol <- substring(dfycol, 1, 4)} else{print("x")}
        dfy[i,dfycol] <- dfz2[i,ii]
      } 
      if (grepl(pattern = "_DH", x = dfz[i,ii]) == TRUE) {
        p1 <- scan(text = xed[["FEATURES"]][[dfx[i,2]]][[dfx[i,ii]]], what = "")
        p2 <- sub(".*nrpspksdomains_","", p1[3])
        p2 <- gsub("\\(|\\).", "", p2)
        p3 <- strsplit(p1[length(p1)], split= "\\.")
        dfycol <- paste("DH", p3[[1]][[length(p3[[1]])]] ,sep = ".")
        if(nchar(dfycol) >= 5){dfycol <- substring(dfycol, 1, 4)} else{print("x")}
        dfy[i,dfycol] <- dfz2[i,ii]
      }
      if (grepl(pattern = "_ER", x = dfz[i,ii]) == TRUE) {
        p1 <- scan(text = xed[["FEATURES"]][[dfx[i,2]]][[dfx[i,ii]]], what = "")
        p2 <- sub(".*nrpspksdomains_","", p1[3])
        p2 <- gsub("\\(|\\).", "", p2)
        p3 <- strsplit(p1[length(p1)], split= "\\.")
        dfycol <- paste("ER", p3[[1]][[length(p3[[1]])]] ,sep = ".")
        if(nchar(dfycol) >= 5){dfycol <- substring(dfycol, 1, 4)} else{print("x")}
        dfy[i,dfycol] <- dfz2[i,ii]
      } 
      if (grepl(pattern = "_Docking_Nterm", x = dfz[i,ii]) == TRUE) {
        p1 <- scan(text = xed[["FEATURES"]][[dfx[i,2]]][[dfx[i,ii]]], what = "")
        p2 <- sub(".*nrpspksdomains_","", p1[3])
        p2 <- gsub("\\(|\\).", "", p2)
        p3 <- strsplit(p1[length(p1)], split= "\\.")
        dfycol <- "N.term.dock"
        dfy[i,dfycol] <- dfz2[i,ii]
      }
      if (grepl(pattern = "_Docking_Cterm", x = dfz[i,ii]) == TRUE) {
        p1 <- scan(text = xed[["FEATURES"]][[dfx[i,2]]][[dfx[i,ii]]], what = "")
        p2 <- sub(".*nrpspksdomains_","", p1[3])
        p2 <- gsub("\\(|\\).", "", p2)
        p3 <- strsplit(p1[length(p1)], split= "\\.")
        dfycol <- "C.term.dock"
        dfy[i,dfycol] <- dfz2[i,ii]
      } 
      if (grepl(pattern = "_Thioesterase", x = dfz[i,ii]) == TRUE) {
        p1 <- scan(text = xed[["FEATURES"]][[dfx[i,2]]][[dfx[i,ii]]], what = "")
        p2 <- sub(".*nrpspksdomains_","", p1[3])
        p2 <- gsub("\\(|\\).", "", p2)
        p3 <- strsplit(p1[length(p1)], split= "\\.")
        dfycol <- paste("TE", p3[[1]][[length(p3[[1]])]] ,sep = ".")
        if(nchar(dfycol) >= 5){dfycol <- substring(dfycol, 1, 4)} else{print("x")}
        dfy[i,dfycol] <- dfz2[i,ii]
      }    
      if (is.na(dfz[i,ii])) {NULL}
    }
  }
  print("WARNING: annotations for PP are being treated as ACPs in dfy")
  
  
  #convert coordinates in dfy to amino acid sequences#
  #Add full AA sequences to dfy
  dfy$gene_AA <- NA
  countr <- 0
  for(i in dfx$FramePosition)
  {
    countr <- countr + 1
    dfy[countr, "gene_AA"] <- x[["FEATURES"]][[i]][["translation"]]
    print(i)
    print(countr)
  }
  
  #loop for taking aminoacids from dfy and reapplying to dfy2
  dfy2 <- dfy
  
  for(i in 1:nrow(dfy)){
    for(ii in 1:(length(dfy)-1))
    {
      if(is.na(dfy[i,ii]) == TRUE){next}
      test <- AAString(x = dfy$gene_AA[i])
      domain_range_lower <- strsplit(x = dfy[i,ii], split = "-")[[1]][[1]]
      domain_range_upper <- strsplit(x = dfy[i,ii], split = "-")[[1]][[2]]
      test <- letter(test, domain_range_lower:domain_range_upper)
      dfy2[i,ii] <- test
    }
  }
  
  write.csv(dfy2, file = paste(date, "_", bgc_name, ".csv", sep = ""))
  framesize <- c(framesize, paste(bgc_name, "columns:" , ncol(dfy2), sep = "_"))
 
})

##Post .gbk extraction analysis##
print(framesize)
length(framesize)
length(df3)

#list ommited BGCs#
missed_bgc <- c()
for(i in 1:length(df3[,1])){
  if(grepl(df3[i,1], paste(framesize, collapse = " ")) == FALSE){missed_bgc <- c(missed_bgc, df3[i,1])}
}
write.csv(missed_bgc, file= paste(date, "_missed_bgcs.csv"))

#List bgcs with non-standard files for manual correction.#
#if a frame is larger than 83 then something has gone wrong with the domain assignment.#
column_error_bgcs <- c()
for(i in 1:length(framesize)){
  if(grepl("_83", paste(framesize[i], collapse = " ")) == FALSE){column_error_bgcs <- c(column_error_bgcs, framesize[i])}
}
write.csv(column_error_bgcs, file= paste(date, "_column_errors.csv"))

#loop for creating a mastersheet#
bgc_outputs <- list.files()
bgc_outputs <- bgc_outputs[grepl(paste0(date, "_BGC"), bgc_outputs)]
mastersheet <- data.frame()
for(i in bgc_outputs){
  ii <- read.csv(i)
  if (ncol(ii) != 84) {bgc_outputs <- bgc_outputs[-grep(pattern = i, x = bgc_outputs)]  
  print(paste0("WARNING! Manual curation required for: ", i))} 

}
  
for(i in bgc_outputs){
  mastersheet <- rbind(mastersheet,ii)
  }
write.csv(mastersheet, file = paste0(date, "_mastersheet.csv"))


