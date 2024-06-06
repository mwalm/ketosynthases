### This script scrapes clusterCAD for SMILES and AT-types ###

library(rvest)
library(polite)
library(tidyverse)

## Inputs ##
date <- "YYYYMMDD"

# For reviewed only #
x1 <- bow("https://clustercad.jbei.org/pks/")
# Reviewed and unreviewed # 
x3 <- bow("https://clustercad.jbei.org/pks/all/")

# Return BGC numbers #
result <- scrape(x1,)%>%
  html_nodes('tr')%>%
  html_node('td')%>%
  html_text2()
head(result)
df1 <- result[2:length(result)]
# Remove decimal on bgc number # 
df3 <- data.frame(substr(df1,1,nchar(df1)-2))
colnames(df3) <- c('MiBIG BGC number')



## Smiles extraction test script ##
x2 <- bow("https://clustercad.jbei.org/pks/BGC0000001.1")
r2 <- scrape(x2, )%>%
  html_nodes('ul')%>%
  html_nodes('li > div > div > img')
kl <- data.frame(as.character(r2))

# r3 targets buttons tooltips #
# this can be refitted to target other domain button tips #
r3 <- scrape(x2, )%>%
  html_nodes('ul')%>%
  html_nodes('li > div > div > div')
substrate_frame <- data.frame(as.character(r3))

for(ii in kl){
  kl2 <- scan(text = ii, what = "")
}

kl3 <- data.frame(kl2[grepl("data-smiles=", kl2 )])
colnames(kl3) <- "Smiles"
rownames(kl3) <- gsub(":", "", 
                      paste(
                        kl2[c(FALSE, grepl("subunit", kl2))], 
                        kl2[c(FALSE, grepl("module", kl2))], 
                        sep = '_Mod.'
                      )
)
#Tidy Smiles#
kl3$Smiles <- gsub("\"", "", gsub("data-smiles=", "", as.character(kl3$Smiles)))

#Add compound name to dataframe#
kl3$BGC_Compound <- gsub("alt=\"", "", kl2[grepl("alt=", kl2)])


for(i in substrate_frame){
  sf1 <- scan(text = i, what = "")
}
sf2 <- data.frame(sf1[grep("substrate", sf1) + 1])
colnames(sf2) <- c("AT-type")
sf2$`AT-type` <- gsub(",$", "", sf2$`AT-type`)

kl3 <- data.frame(kl3, sf2)


test <- c(df1[1:length(df1)])
outputframe <- data.frame()
##smiles extraction loop##
for(bgcpage in test)try({
  address <- paste("https://clustercad.jbei.org/pks/", bgcpage, sep = "")
  x2 <- bow(noquote(address))
  r2 <- scrape(x2, )%>%
    html_nodes('ul')%>%
    html_nodes('li > div > div > img')
  kl <- data.frame(as.character(r2))
  
  r3 <- scrape(x2, )%>%
    html_nodes('ul')%>%
    html_nodes('li > div > div > div')
  substrate_frame <- data.frame(as.character(r3))
  
  for(ii in kl){
    kl2 <- scan(text = ii, what = "")
  }
  
  for(i in substrate_frame){
    sf1 <- scan(text = i, what = "")
  }
  
  kl3 <- data.frame(kl2[grepl("data-smiles=", kl2 )])
  colnames(kl3) <- "Smiles"
  rownames(kl3) <- gsub(":", "", 
                        paste(
                          kl2[c(FALSE, grepl("subunit", kl2))], 
                          kl2[c(FALSE, grepl("module", kl2))], 
                          sep = '_Mod.'
                        )
  )
  
  sf2 <- data.frame(sf1[grep("substrate", sf1) + 1])
  colnames(sf2) <- c("AT-type")
  sf2$`AT-type` <- gsub(",$", "", sf2$`AT-type`)
  
  
  kl3$Smiles <- gsub("\"", "", gsub("data-smiles=", "", as.character(kl3$Smiles)))
  kl3$BGC_Compound <- gsub("alt=\"", "", kl2[grepl("alt=", kl2)])
  
  if(nrow(kl3) != nrow(sf2)){sf2 <- rbind(NA,sf2)}
  if(nrow(kl3) != nrow(sf2)){print(paste0("error" , bgcpage))}
  if(nrow(kl3) != nrow(sf2)){next}
  kl3 <- data.frame(kl3, sf2)
  
  outputframe<- rbind.data.frame(outputframe,kl3)
})

write.csv(outputframe, file = paste0(date, "_Smiles_and_AT_reviewed.csv"))

# Checks: # 
final<- read.csv(paste0(date, "_Smiles_and_AT_reviewed.csv"))
length(unique(unlist(strsplit(final$BGC_Compound, " "))))



## Module unification script ##

uniques <- unique(gsub(pattern = "\\..*", "", final$X))  
final$renamed_modules <- NA
for (i in uniques) {
  genename <- final[grep(pattern = i, x = final$X),1]
  rownumbers <- grep(pattern = i, x = final$X)
 
  if(sub(".*\\.", "", (genename[1])) == 0){countr <- 0}
  else{countr <- 1}
  countr2 <- 1
  for(ii in final[rownumbers, 1]){
    final[rownumbers[countr2],5] <- paste0(gsub(pattern = "\\..*", ".", ii), countr)
    countr <- countr + 1
    countr2 <- countr2 + 1 
  }
 }

write.csv(final, file = paste0(date, "_Smiles_and_AT_reviewed.csv"))




