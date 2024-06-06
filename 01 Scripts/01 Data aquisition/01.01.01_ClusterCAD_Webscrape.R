library(rvest)
library(polite)
library(tidyverse)


# For reviewed only #
x1 <- bow("https://clustercad.jbei.org/pks/")
# Reviewed and unreviewed #
#x1 <- bow("https://clustercad.jbei.org/pks/all/")

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

## write df3 as csv for later ##




## 1. Smiles extraction test ##
  x2 <- bow("https://clustercad.jbei.org/pks/BGC0000001.1")
  r2 <- scrape(x2, )%>%
    html_nodes('ul')%>%
    html_nodes('li > div > div > img')
  print.AsIs(r2)
  kl <- data.frame(as.character(r2))
  
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
  # Tidy Smiles #
  kl3$Smiles <- gsub("\"", "", gsub("data-smiles=", "", as.character(kl3$Smiles)))
  
  #Add compound name to dataframe#
  kl3$BGC_Compound <- gsub("alt=\"", "", kl2[grepl("alt=", kl2)])
  test <- c(df1)
  outputframe <- data.frame()


## 2. smiles extraction loop ##
for(bgcpage in test){
  address <- paste("https://clustercad.jbei.org/pks/", bgcpage, sep = "")
  x2 <- bow(noquote(address))
  r2 <- scrape(x2, )%>%
    html_nodes('ul')%>%
    html_nodes('li > div > div > img')
  kl <- data.frame(as.character(r2))
  
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
  kl3$Smiles <- gsub("\"", "", gsub("data-smiles=", "", as.character(kl3$Smiles)))
  kl3$BGC_Compound <- gsub("alt=\"", "", kl2[grepl("alt=", kl2)])
  outputframe<- rbind.data.frame(outputframe,kl3)
}

### Write outputframe as csv with write.csv() ###
  
## 3. Output checks ##
final<- read.csv("20231107_Smiles_reviewed.csv")
length(unique(unlist(strsplit(final$BGC_Compound, " "))))

