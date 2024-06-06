### AT explainer histogram ###

## bin graphs ##
file.choose()
file <- read.csv("C:\\Users\\q31032mw\\The University of Manchester Dropbox\\Max Walmsley\\Max\\17_ML_Project\\Figures scripts\\Methylation_KS\\20240229_methylation_state_graphViz.csv")
substrates_sheet <- read.csv("C:\\Users\\q31032mw\\The University of Manchester Dropbox\\Max Walmsley\\Max\\17_ML_Project\\Webscraping\\ClusterCAD Webscraping\\20230322_AT_SMILES_unified.csv")

mmal <- substrates_sheet[which(substrates_sheet$AT_substrate == "mmal"),]
mal <- substrates_sheet[which(substrates_sheet$AT_substrate == "mal"),]

mmal_list <- c()
for(i in mmal$Gene){
  mmal_list <- c(mmal_list, which(file$Gene..AT.== i))
}

mal_list <- c()
for(i in mal$Gene){
  mal_list <- c(mal_list, which(file$Gene..AT.== i))
}

mal_nodes <- file[mal_list,]
mmal_nodes <- file[mmal_list,]

# Sort nodes for Mal #
mal_nodes_string <- c()
for (i in mal_nodes$Residues.strong) {
  mal_nodes_string <- c(mal_nodes_string, i)
}
for (i in mal_nodes$Residues.weak) {
  mal_nodes_string <- c(mal_nodes_string, i)
}

all_elements <- c()
for (str in mal_nodes_string) {
  elements <- as.numeric(unlist(strsplit(str, ",")))
  all_elements <- c(all_elements, elements)
}
mal_sorted_elements <- sort(all_elements)
print(mal_sorted_elements)




# sort nodes for mmal #
mmal_nodes_string <- c()
for (i in mmal_nodes$Residues.strong) {
  mmal_nodes_string <- c(mmal_nodes_string, i)
}
for (i in mmal_nodes$Residues.weak) {
  mmal_nodes_string <- c(mmal_nodes_string, i)
}

all_elements <- c()
for (str in mmal_nodes_string) {
  elements <- as.numeric(unlist(strsplit(str, ",")))
  all_elements <- c(all_elements, elements)
}
mmal_sorted_elements <- sort(all_elements)
print(mmal_sorted_elements)


hist(mal_sorted_elements,breaks = max(mal_sorted_elements), xlab = "Residue number", col = "red")
box(which = "plot", bty = "l")
axis(pos = 0)

hist(mmal_sorted_elements[1:448],breaks = max(mmal_sorted_elements[1:448]), xlab = "Residue number", col = "blue")
box(which = "plot", bty = "l")
axis(pos = 0)

library(plotly)

fig <- plot_ly(x = mal_sorted_elements, type = "histogram", name = "Mal", nbinsx = 900) 
fig <- fig %>% add_histogram(x = mmal_sorted_elements[1:448], alpha = 0.8, type = "histogram", name = "Mmal")
fig <- fig %>% layout(barmode = "overlay")
fig <- fig %>% layout(
  yaxis = list(title = "Frequency", tick0 = 0, dtick = 1, tickfont = list(size = 16)),
  xaxis = list(title = "Residue number", tick0 = 0, dtick = 20, tickfont = list(size = 16))
)

fig

