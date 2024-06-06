## Cube faces with z depths ## 

## Cube face script ##
# idea is to feed "photos" of the cube, where each photo is a cube face. Values = depth dimension##
# 0. Read the pdb files into a large list/array
{
  ## Dependencies ## 
  library(bio3d)
  library(tensorflow)
  library(keras)
  
  ## Inputs ## 
  #route <- "/home/q31032mw/Dropbox (The University of Manchester)/Max/17_ML_Project/Rank_1/Aligned/"
  #route <- "/home/q31032mw/Dropbox (The University of Manchester)/Max/17_ML_Project/Rank_1/Aligned_big_set3/"
  route <- "/home/q31032mw/Dropbox (The University of Manchester)/Max/17_ML_Project/Rank_1/Aligned_KS_2/"
  #route <- "/home/q31032mw/Dropbox (The University of Manchester)/Max/17_ML_Project/Rank_1/Unaligned_pdb_files/"
  
  
  ### 1. Collate aligned domains into a unified data frame ###
  # Extract file names in route # 
  files_list <- list.files(path = route)
  
  ## Generate total_data file structure and populate ## 
  # Data architecture:
  #     1. List: genes
  #          2. Data frame: Atom coordinates of protein structure
  
  total_data <- list()
  dcs <- c(4,6, 8, 14, 16)
  #dcs are empty columns
  
  for (i in files_list){
    name_in_list <- paste0(route, i)
    pdb <- read.pdb(name_in_list)
    pdb2 <- cbind(pdb[1], as.data.frame(pdb[3]))
    colnames(pdb2) <- gsub(x = colnames(pdb2),pattern = "atom.", replacement = "") 
    total_data[[i]] <- pdb2[-dcs]
  }
  
  
  # Define labels (substrate) # 
  Substrates <-  read.csv("/home/q31032mw/Dropbox (The University of Manchester)/Max/17_ML_Project/Webscraping/ClusterCAD Webscraping/20230322_AT_SMILES_unified.csv")
  #training_names <- gsub(x = names(total_data), pattern = "\\.pdb_aligned.pdb", replacement = "")
  training_names <- gsub(x = names(total_data), pattern = "\\.pdb_aligned.pdb", replacement = "")
  
  refined_substrates <- c()
  for (i in training_names) {
    if (sum(grepl(pattern = i, x = Substrates[,2])) >=2) {
      print(i)
      print(sum(grepl(pattern = i, x = Substrates[,2])))
      i2 <- grep(pattern = i, x = Substrates[,2])
      refined_substrates <-  c(refined_substrates,i2[1])
    } 
    if (sum(grepl(pattern = i, x = Substrates[,2])) ==0) {
      print(paste0("missing: ", i))
    } 
    else{
      refined_substrates <-  c(refined_substrates, grep(pattern = i, x = Substrates[,2]))
    }
    
  }
  y <- Substrates[refined_substrates, c(2,5)]
  
  # Report Mal-Mmal bias # 
  print(paste0("mmal: ", sum(grepl(pattern = "mmal", y$AT_substrate))))
  print(paste0("mal: ", (nrow(y) - sum(grepl(pattern = "mmal", y$AT_substrate)))))
  
  # Dump excess mal from y#
  if(sum(grepl(pattern = "mmal", y$AT_substrate)) <= (nrow(y) - sum(grepl(pattern = "mmal", y$AT_substrate)))){
    dump_number <- (nrow(y) - sum(grepl(pattern = "mmal", y$AT_substrate))) - sum(grepl(pattern = "mmal", y$AT_substrate))
    mal_list <- y[!grepl(pattern = "mmal", y$AT_substrate),]
    dump_list <- paste0(sample(x = mal_list[,1], size = dump_number))
  }
  rm_list <- c()
  for(i in dump_list){
    rm_list <- c(rm_list, grep(pattern = i, y$Gene))
  }
  y <- y[-rm_list,]
  
  
  print(paste0("mmal: ", sum(grepl(pattern = "mmal", y$AT_substrate))))
  print(paste0("mal: ", (nrow(y) - sum(grepl(pattern = "mmal", y$AT_substrate)))))
  
  # Remove the dumped mal from training data ##
  for (i in dump_list) {
    total_data <- total_data[-grep(pattern = i, x = names(total_data))]
  }
  training_names <- gsub(x = names(total_data), pattern = "\\.pdb", replacement = "")
  
  ## Data exploration ##
  # You need graphs that display the ranges of values before and after pre-processing. This will pick up outliers.# 
  
  
  ## Data pre-processing ##
  # Assumptions: atoms coordinates can have both positive and negative values. 
  
  #assess range for atom coordinates
  max_x <- 0
  max_y <- 0
  max_z <- 0
  
  metrics_x <- data.frame(matrix(ncol = 2, nrow = 0))
  metrics_y <- data.frame(matrix(ncol = 2, nrow = 0))
  metrics_z <- data.frame(matrix(ncol = 2, nrow = 0))
  
  for(i in 1:length(total_data)){
    if(max_x <= max(total_data[[i]][[6]]))
    {max_x <- max(total_data[[i]][[6]])}
    
    if(max_y <= max(total_data[[i]][[7]]))
    {max_y <- max(total_data[[i]][[7]])}
    
    if(max_z <= max(total_data[[i]][[8]]))
    {max_z <- max(total_data[[i]][[8]])}
    
    metrics_x <- rbind(metrics_x, c(range(total_data[[i]][[6]])))
    metrics_y <- rbind(metrics_y, c(range(total_data[[i]][[7]])))
    metrics_z <- rbind(metrics_z, c(range(total_data[[i]][[8]])))
    
  }
  colnames(metrics_x) <- c("LowerX", "UpperX")
  colnames(metrics_y) <- c("LowerY", "UpperY")
  colnames(metrics_z) <- c("LowerZ", "UpperZ")
  metrics_x$RowName <- cbind(training_names)
  metrics_all <- cbind(metrics_x, metrics_y, metrics_z)
  
  # Plotting
  # Load the required libraries
  library(ggplot2)
  library(plotly)
  
  p <- ggplot(metrics_all, aes(x = RowName)) +
    geom_point(aes(y = LowerX, fill = "LowerX"), position =  "dodge", stat = "identity", color = "black") +
    geom_point(aes(y = UpperX, fill = "UpperX"), position =  "dodge", stat = "identity", color = "black") +
    geom_point(aes(y = LowerY, fill = "LowerY"), position =  "dodge", stat = "identity", color = "black") +
    geom_point(aes(y = UpperY, fill = "UpperY"), position =  "dodge", stat = "identity", color = "black") +
    geom_point(aes(y = LowerZ, fill = "LowerZ"), position =  "dodge", stat = "identity", color = "black") +
    geom_point(aes(y = UpperZ, fill = "UpperZ"), position =  "dodge", stat = "identity", color = "black") +
    scale_fill_manual(values = c("LowerX" = "red", "UpperX" = "red", "LowerY" = "green", "UpperY" = "green", "LowerZ" = "blue", "UpperZ" = "blue")) +
    labs(x = "Row Name", y = "Value")
  
  # Convert ggplot to plotly
  coordinate_range_plot <- ggplotly(p, tooltip = c("x", "y"))
  
  # Show the plot
  coordinate_range_plot
  
  min_x <- 0
  min_y <- 0
  min_z <- 0
  
  for(i in 1:length(total_data)){
    if(min_x >= min(total_data[[i]][[6]]))
    {min_x <- min(total_data[[i]][[6]])}
    
    if(min_y >= min(total_data[[i]][[7]]))
    {min_y <- min(total_data[[i]][[7]])}
    
    if(min_z >= min(total_data[[i]][[8]]))
    {min_z <- min(total_data[[i]][[8]])}
  }
  
  print(paste("Maximum XYX:", max_x, max_y, max_z))
  print(paste("Minimum XYX:", min_x, min_y, min_z))
}
# 1. move numbers into a positive only space. 
# Find the minimum of X, Y and Z. Add it to all pdbs. Decide on rounding number. 
rounding_number <- 2L
faces <- 3L
{
  for (i in 1:length(total_data)) {
    total_data[[i]][[6]] <-  (total_data[[i]][[6]]) - min_x
    total_data[[i]][[7]] <-  (total_data[[i]][[7]]) - min_y
    total_data[[i]][[8]] <-  (total_data[[i]][[8]]) - min_z
  }
  
  
  
  max_x <- 0
  max_y <- 0
  max_z <- 0
  
  metrics_x <- data.frame(matrix(ncol = 2, nrow = 0))
  metrics_y <- data.frame(matrix(ncol = 2, nrow = 0))
  metrics_z <- data.frame(matrix(ncol = 2, nrow = 0))
  
  for(i in 1:length(total_data)){
    if(max_x <= max(total_data[[i]][[6]]))
    {max_x <- max(total_data[[i]][[6]])}
    
    if(max_y <= max(total_data[[i]][[7]]))
    {max_y <- max(total_data[[i]][[7]])}
    
    if(max_z <= max(total_data[[i]][[8]]))
    {max_z <- max(total_data[[i]][[8]])}
    
    metrics_x <- rbind(metrics_x, c(range(total_data[[i]][[6]])))
    metrics_y <- rbind(metrics_y, c(range(total_data[[i]][[7]])))
    metrics_z <- rbind(metrics_z, c(range(total_data[[i]][[8]])))
    
  }
  colnames(metrics_x) <- c("LowerX", "UpperX")
  colnames(metrics_y) <- c("LowerY", "UpperY")
  colnames(metrics_z) <- c("LowerZ", "UpperZ")
  metrics_x$RowName <- cbind(training_names)
  metrics_all <- cbind(metrics_x, metrics_y, metrics_z)
  
  min_x <- 0
  min_y <- 0
  min_z <- 0
  
  for(i in 1:length(total_data)){
    if(min_x >= min(total_data[[i]][[6]]))
    {min_x <- min(total_data[[i]][[6]])}
    
    if(min_y >= min(total_data[[i]][[7]]))
    {min_y <- min(total_data[[i]][[7]])}
    
    if(min_z >= min(total_data[[i]][[8]]))
    {min_z <- min(total_data[[i]][[8]])}
  }
  
  print(paste("Maximum XYX:", max_x, max_y, max_z))
  print(paste("Minimum XYX:", min_x, min_y, min_z))
  
  
  # convert coordinates to 0-1 range
  atom_refactor <-max(c(max_x, max_y, max_z, (min_x*-1), (min_y*-1),(min_z*-1))) 
  
  refined_total_data <- total_data
  # Apply atom_refactor to coordinates. 
  for (i in 1:length(total_data)) {
    refined_total_data[[i]][[6]] <-  round(digits = rounding_number, x = ((total_data[[i]][[6]])/atom_refactor))
    refined_total_data[[i]][[7]] <-  round(digits = rounding_number, x = ((total_data[[i]][[7]])/atom_refactor))
    refined_total_data[[i]][[8]] <-  round(digits = rounding_number, x = ((total_data[[i]][[8]])/atom_refactor)) 
    refined_total_data[[i]][[10]] <-  (total_data[[i]][[10]])/100
  }
  
  # Partition data
  #######################################################################################################
  # Partition training and test data ## 
  # Align y with refined_total_data # 
  # Creates 3rd column in y that corresponds the ordering in total_data # 
  for (i in 1:nrow(y)) {
    y[i,3] <- grep(pattern = y[i,1], names(refined_total_data))
  }
  
  # Define partition factor # 
  partition_ratio <- 0.1
  
  # Partition the training data from test data # 
  # Partition evenly distributes between mmal and mal # 
  test_data_length <- length(refined_total_data)*partition_ratio
  y_mal <- y[!grepl("mmal", y[,2]),][sample(1:nrow(y[!grepl("mmal", y[,2]),]), size = test_data_length/2),]
  y_mmal <- y[grepl("mmal", y[,2]),][sample(1:nrow(y[grepl("mmal", y[,2]),]), size = test_data_length/2),]
  partition <- as.numeric(c(y_mal[,3], y_mmal[,3]))
  test_data <- refined_total_data[partition]
  training_data <- refined_total_data[-partition]
  
  # Labels (y) lists and encoding # 
  y_train <- rbind(y_mal, y_mmal)
  removal_list <- c()
  for (i in y_train$Gene) {
    removal_list <- c(removal_list, grep(pattern = i, x = y$Gene))
  }
  y_test <- y[-removal_list,]
  # Convert labels to numeric # 
  y_train <- to_categorical(as.numeric(factor(y_train[1:nrow(y_train),2])))
  y_test <- to_categorical(as.numeric(factor(y_test$AT_substrate)))
  # to_categorical seems to be adding an extra column (V1), this is removed below but is a bit of a bodge #
  # Possible this is a disconnect between where python and R start counting? # 
  y_train <- y_train[,-1]
  y_test <- y_test[,-1]
  # You fucked up and mixed the training and test data for y. This fudge repairs that.   
  y_train2 <- y_test
  y_test2 <- y_train
  y_train <- y_train2
  y_test <- y_test2
  
  # TIDY UP # 
  columns_to_ignore <- c(1:5, 9:70)
  final_test_data <- test_data
  for (i in 1:length(test_data)) {
    final_test_data[[i]] <- final_test_data[[i]][,-columns_to_ignore]
  }
  final_training_data <- training_data
  for (i in 1:length(training_data)) {
    final_training_data[[i]] <- final_training_data[[i]][,-columns_to_ignore]
  }  
  x_test <- final_test_data
  x_train <- final_training_data
  
  #############################################################################################################
}

# 2.  Create and populate a 100x100x100 cube 

# Define the atomic coordinate cube boundaries and row density. 
# this is a product of the rounding_number. Lower rounding needs bigger cubes. 
max_atom <- 0L
for(i in 1:length(final_test_data)){
  if(max_atom <= max(final_test_data[[i]])){
    max_atom <- max(final_test_data[[i]]) 
  }
}
cube_boundary <- 10^rounding_number
surfaces <- 6L
x_train1 <- x_train
x_test1 <- x_test
## XY front face ## 
for (i in 1:length(names(x_train))) {
  data <- as.matrix(x_train[[i]][,1:3])
  coords <- t(apply(data, 1, function(row) as.integer(row * 10^rounding_number)))
  empty_array <- array(data= 0, dim = c(faces, cube_boundary, cube_boundary)) 
  
  for(ii in 1:nrow(coords)){
    empty_array[1,coords[ii, 1],coords[ii, 2]] <- (coords[ii, 3]/(10^rounding_number))
  }
  x_train1[[i]] <- empty_array
}

for (i in 1:length(names(x_test))) {
  data <- as.matrix(x_test[[i]][,1:3])
  coords <- t(apply(data, 1, function(row) as.integer(row * 10^rounding_number)))
  empty_array <- array(data= 0, dim = c(faces, cube_boundary, cube_boundary)) 
  
  for(ii in 1:nrow(coords)){
    empty_array[1, coords[ii, 1],coords[ii, 2]] <- (coords[ii, 3]/(10^rounding_number))
  }
  x_test1[[i]] <- empty_array
}

## XZ face ##
for (i in 1:length(names(x_train))) {
  data <- as.matrix(x_train[[i]][,1:3])
  coords <- t(apply(data, 1, function(row) as.integer(row * 10^rounding_number)))
  empty_array <- array(data= 0, dim = c(cube_boundary, cube_boundary)) 
  
  for(ii in 1:nrow(coords)){
    empty_array[coords[ii, 1],coords[ii, 3]] <- (coords[ii, 2]/(10^rounding_number))
  }
  x_train1[[i]][2,,] <- empty_array
}

for (i in 1:length(names(x_test))) {
  data <- as.matrix(x_test[[i]][,1:3])
  coords <- t(apply(data, 1, function(row) as.integer(row * 10^rounding_number)))
  empty_array <- array(data= 0, dim = c(cube_boundary, cube_boundary)) 
  
  for(ii in 1:nrow(coords)){
    empty_array[coords[ii, 1],coords[ii, 3]] <- (coords[ii, 2]/(10^rounding_number))
  }
  x_test1[[i]][2,,] <- empty_array
}

## YZ face ###
for (i in 1:length(names(x_train))) {
  data <- as.matrix(x_train[[i]][,1:3])
  coords <- t(apply(data, 1, function(row) as.integer(row * 10^rounding_number)))
  empty_array <- array(data= 0, dim = c(cube_boundary, cube_boundary)) 
  
  for(ii in 1:nrow(coords)){
    empty_array[coords[ii, 2],coords[ii, 3]] <- (coords[ii, 1]/(10^rounding_number))
  }
  x_train1[[i]][3,,] <- empty_array
}

for (i in 1:length(names(x_test))) {
  data <- as.matrix(x_test[[i]][,1:3])
  coords <- t(apply(data, 1, function(row) as.integer(row * 10^rounding_number)))
  empty_array <- array(data= 0, dim = c(cube_boundary, cube_boundary)) 
  
  for(ii in 1:nrow(coords)){
    empty_array[coords[ii, 1],coords[ii, 3]] <- (coords[ii, 1]/(10^rounding_number))
  }
  x_test1[[i]][3,,] <- empty_array
}

#CNN
library(keras)

# Convert the list of genes to a single array
x_train2 <- array_reshape(unlist(x_train1), dim = c(length(x_train1), faces, cube_boundary, cube_boundary))
x_test2 <- array_reshape(unlist(x_test1), dim = c(length(x_test1), faces, cube_boundary, cube_boundary))

# Define the CNN model
library(keras)

# Rearrange dimensions of x_train2
x_train2 <- array_reshape(x_train2, c(dim(x_train2)[1], faces, cube_boundary, cube_boundary))
x_train2 <- aperm(x_train2, c(1, 3, 4, 2))

# Create the model
model <- keras_model_sequential()
set.seed(42)

# Add the layers
model %>%
  layer_conv_2d(filters = 32, kernel_size = c(3, 3), activation = "relu", input_shape = c(cube_boundary, cube_boundary, faces)) %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  layer_conv_2d(filters = 64, kernel_size = c(3, 3), activation = "relu") %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  layer_flatten() %>%
  layer_dense(units = 64, activation = "relu") %>%
  #layer_dropout(rate = 0.2) %>%
  layer_dense(units = 2, activation = "softmax")

# Compile the model
model %>% compile(
  loss = "categorical_crossentropy",
  optimizer = "adam",
  metrics = c("accuracy")
)

# Train the model
model %>% fit(
  x_train2, y_train,
  epochs = 30,
  batch_size = 35,
  validation_split = 0.2
)



x_test3 <- array_reshape(x_test2, c(dim(x_test2)[1], faces, cube_boundary, cube_boundary))
x_test3 <- aperm(x_test2, c(1, 3, 4, 2))


# Evaluate the model
model %>% evaluate(x_test3, y_test)

