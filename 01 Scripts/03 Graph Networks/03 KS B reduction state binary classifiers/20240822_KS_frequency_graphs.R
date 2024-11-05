file.choose()
route <-   "C:\\Users\\q31032mw\\The University of Manchester Dropbox\\Max Walmsley\\Max\\17_ML_Project\\"

DHvsKR_DH <- paste0(route, "Figures scripts\\20240913_KRvsDH_GCN_corrected_labels\\","20240913_KRvsDH_binary_GCN_corrected_labels_DH_netx_score_nysC_Mod.2.pdb")
DHvsKR_KR <- paste0(route, "Figures scripts\\20240913_KRvsDH_GCN_corrected_labels\\","20240913_KRvsDH_binary_GCN_corrected_labels_KR_netx_score_SeEryAI_Mod.2.pdb")

DHvsER_DH <- paste0(route, "Figures scripts\\20240914_DHvsER_GCN_corrected_labels\\", "20240917_DHvsER_binary_GCN_corrected_labels_DH_netx_score_NysC_Mod.2.pdb")
DHvsER_ER <- paste0(route, "Figures scripts\\20240914_DHvsER_GCN_corrected_labels\\","20240917_DHvsER_binary_GCN_corrected_labels_ER_netx_score_SeEryAIII_Mod.2.pdb")
  
DHvsNR_DH <- paste0(route, "Figures scripts\\20240917_DHvsNR_GCN_corrected_labels\\","20240919_DHvsNR_binary_GCN_corrected_labels_DH_netx_score_NysC_Mod.2.pdb")
DHvsNR_NR <- paste0(route, "Figures scripts\\20240917_DHvsNR_GCN_corrected_labels\\","20240919_DHvsNR_binary_GCN_corrected_labels_NR_netX_score_AceP4_Mod.6.pdb")

KRvsNR_KR <- paste0(route, "Figures scripts\\20240923_KRvsNR_corrected_labels\\","20240923_KRvsNR_binary_GCN_corrected_labels_KR_netx_score_SeEryAII_Mod.2.pdb") 
KRvsNR_NR <- paste0(route, "Figures scripts\\20240923_KRvsNR_corrected_labels\\","20240923_KRvsNR_binary_GCN_corrected_labels_NR_netX_score_AceP4_Mod.6.pdb")


l <- c(DHvsKR_DH, DHvsKR_KR, DHvsNR_DH, DHvsNR_NR, DHvsER_DH, DHvsER_ER, KRvsNR_KR, KRvsNR_NR)
l2 <- c("DHvsKR_DH", "DHvsKR_KR", "DHvsNR_DH", "DHvsNR_NR", "DHvsER_DH", "DHvsER_ER", "KRvsNR_KR", "KRvsNR_NR")
library(bio3d)
countr <- 0L
for (i in l) {
  countr <- countr + 1L
  pdb <- read.pdb(i)
  CA <- which(pdb$atom$elety == "CA")
  b <- pdb$atom$b[CA]
  atomnum <- pdb$atom$resno[CA]
  df <- data.frame(atomnum, b)
  assign( l2[countr] ,df[1:max(df[,1]),])
}

library(plotly)

l3 <- list(DHvsKR_DH, DHvsKR_KR, DHvsNR_DH, DHvsNR_NR, DHvsER_DH, DHvsER_ER, KRvsNR_KR, KRvsNR_NR)

countr <- 0L
for (i in l3) {
  countr <- countr + 1
  assign(paste0(l2[countr], "_g" ), plot_ly(x = i[,1], y = i[,2],  type = 'scatter', mode = 'lines', name = l2[countr])) 
}


# DHvsKR #
fig <- plot_ly(x = DHvsKR_DH[,1], y = DHvsKR_DH[,2],  type = 'scatter', mode = 'lines', name = "DH") 
fig <- fig %>% add_trace(x = DHvsKR_KR[,1], y = DHvsKR_KR[,2],  type = 'scatter', mode = 'lines', name = "KR") 
fig <- fig %>% layout(
  yaxis = list(
    title = list(text = "Frequency %", font = list(size = 22, family = "Arial", color = "black")),
    tick0 = 0, 
    dtick = 10, 
    tickfont = list(size = 18),
    range = c(0, 100)
  ),
  xaxis = list(
    title = list(text = "Residue number", font = list(size = 22, family = "Arial", color = "black")),
    tick0 = 0, 
    dtick = 20, 
    tickfont = list(size = 18),
    tickangle = 90, 
    range = c(0, 440)
  )
)
fig

# DHvsER #
fig <- plot_ly(x = DHvsER_DH[,1], y = DHvsER_DH[,2],  type = 'scatter', mode = 'lines', name = "DH") 
fig <- fig %>% add_trace(x = DHvsER_ER[,1], y = DHvsER_ER[,2],  type = 'scatter', mode = 'lines', name = "ER") 
fig <- fig %>% layout(
  yaxis = list(
    title = list(text = "Frequency %", font = list(size = 22, family = "Arial", color = "black")),
    tick0 = 0, 
    dtick = 10, 
    tickfont = list(size = 18),
    range = c(0, 100)
  ),
  xaxis = list(
    title = list(text = "Residue number", font = list(size = 22, family = "Arial", color = "black")),
    tick0 = 0, 
    dtick = 20, 
    tickfont = list(size = 18),
    tickangle = 90, 
    range = c(0, 440)
  )
)
fig

# DHvsNR #
fig <- plot_ly(x = DHvsNR_DH[,1], y = DHvsNR_DH[,2],  type = 'scatter', mode = 'lines', name = "DH") 
fig <- fig %>% add_trace(x = DHvsNR_NR[,1], y = DHvsNR_NR[,2],  type = 'scatter', mode = 'lines', name = "NR") 
fig <- fig %>% layout(
  yaxis = list(
    title = list(text = "Frequency %", font = list(size = 22, family = "Arial", color = "black")),
    tick0 = 0, 
    dtick = 10, 
    tickfont = list(size = 18),
    range = c(0, 100)
  ),
  xaxis = list(
    title = list(text = "Residue number", font = list(size = 22, family = "Arial", color = "black")),
    tick0 = 0, 
    dtick = 20, 
    tickfont = list(size = 18),
    tickangle = 90, 
    range = c(0, 440)
  )
)
fig

# KRvsNR #
fig <- plot_ly(x = KRvsNR_KR[,1], y = KRvsNR_KR[,2],  type = 'scatter', mode = 'lines', name = "KR") 
fig <- fig %>% add_trace(x = KRvsNR_NR[,1], y = KRvsNR_NR[,2],  type = 'scatter', mode = 'lines', name = "NR") 
fig <- fig %>% layout(
  yaxis = list(
    title = list(text = "Frequency %", font = list(size = 22, family = "Arial", color = "black")),
    tick0 = 0, 
    dtick = 10, 
    tickfont = list(size = 18),
    range = c(0, 100)
  ),
  xaxis = list(
    title = list(text = "Residue number", font = list(size = 22, family = "Arial", color = "black")),
    tick0 = 0, 
    dtick = 20, 
    tickfont = list(size = 18),
    tickangle = 90, 
    range = c(0, 440)
  )
)
fig

# Combined #
fig <- plot_ly(x = DHvsKR_DH[,1], y = DHvsKR_DH[,2],  type = 'scatter', mode = 'lines', name = "DH-KR") 
#fig <- fig %>% add_trace(x = DHvsKR_KR[,1], y = DHvsKR_KR[,2],  type = 'scatter', mode = 'lines', name = "KR") 
fig <- fig %>% add_trace(x = DHvsER_DH[,1], y = DHvsER_DH[,2],  type = 'scatter', mode = 'lines', name = "DH-ER") 
#fig <- fig %>% add_trace(x = DHvsER_ER[,1], y = DHvsER_ER[,2],  type = 'scatter', mode = 'lines', name = "ER") 
fig <- fig %>% add_trace(x = DHvsNR_DH[,1], y = DHvsNR_DH[,2],  type = 'scatter', mode = 'lines', name = "DH-NR") 
#fig <- fig %>% add_trace(x = DHvsNR_NR[,1], y = DHvsNR_NR[,2],  type = 'scatter', mode = 'lines', name = "NR")
fig <- fig %>% layout(
  yaxis = list(
    title = list(text = "Frequency %", font = list(size = 22, family = "Arial", color = "black")),
    tick0 = 0, 
    dtick = 10, 
    tickfont = list(size = 18),
    range = c(0, 100)
  ),
  xaxis = list(
    title = list(text = "Residue number", font = list(size = 22, family = "Arial", color = "black")),
    tick0 = 0, 
    dtick = 20, 
    tickfont = list(size = 18),
    tickangle = 90, 
    range = c(0, 440)
  )
)
fig




