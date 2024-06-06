### Graphical readout of b scores ###
library(msa)
library(Biostrings)
library(bio3d)
library(ggplot2)

# file to make graph from #
route <- "/home/q31032mw/Dropbox (The University of Manchester)/Max/17_ML_Project/Figures scripts/20240408_SeEryAIII_Mod.1_GCN_binary_DH_netx_score_SeEryAIII_Mod.2.pdb"


# read as pdb #
pdb <- read.pdb(route)

# get sequence # 
pdb_seq <- pdbseq(pdb)
pdb_seq2 <- paste( unlist(pdb_seq), collapse='')

CAs <- which(pdb$atom$elety == "CA")
scores <- pdb$atom$b[CAs]

# Create a data frame with your data
data <- data.frame(letters = strsplit(pdb_seq2, "")[[1]], values = scores)


# Use if dimer only! #
data <- data[1:(nrow(data)/2),]
print("WARNING: DIMER MODE")

# Create a new variable to indicate the position of each letter
data$position <- rep(1:20, length.out = nrow(data))

# Reorder the levels of the factor variable for the y-axis
data$y_factor <- factor(rep(1:ceiling(nrow(data)/20), each = 20, length.out = nrow(data)))
data$y_factor <- factor(data$y_factor, levels = rev(levels(data$y_factor)))

# Create the heatmap with letters forming more square tiles and increased Y-axis tick width using geom_raster() to add black borders
plot <- ggplot(data, aes(x = position, y = y_factor, fill = values)) +
  geom_raster() +
  geom_rect(aes(x = as.numeric(position) - 0.5, y = as.numeric(y_factor) - 0.5, 
                xmin = as.numeric(position) - 0.5, xmax = as.numeric(position) + 0.5, 
                ymin = as.numeric(y_factor) - 0.5, ymax = as.numeric(y_factor) + 0.5), 
            color = "black", fill = NA) +  # Add black borders
  geom_text(aes(label = letters), color = "black", size = 4) +  # Adjust text size
  scale_fill_gradient(low = "ivory", high = "red" ) +  # Adjust the color gradient as needed and add a legend title
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),  # Remove x-axis label
        axis.title.y = element_blank(),  # Remove y-axis label
        legend.position = "bottom",  # Move legend to bottom
        legend.key = element_rect(color = "black", size = 1),  # Add black outline around legend key
        plot.title = element_text(hjust = 0.5)) +  # Center title
  coord_fixed(ratio = 1) +
  labs(title = "xAI average scores: DH-type Ketosynthases (GCN binary). 
       Model: KR vs. DH. 
       Template: SeEryAIII_Mod.1")
plot



ggsave("/home/q31032mw/Dropbox (The University of Manchester)/Max/17_ML_Project/Figures scripts/20240408_KRvsDH_DH_SeEryAIIIMod2_GCN_binary_netx_crossword.jpeg", plot = plot, width = 10, height = 10, units = "in")


