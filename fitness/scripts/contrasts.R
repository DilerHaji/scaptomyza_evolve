# Setting the reference levels to "M" for source and test factors 
tab$PopSpecies <- relevel(tab$PopSpecies, ref = "M")
tab$PhenotypeCageSpecies <- relevel(tab$PhenotypeCageSpecies, ref = "M")

# Making a contrast matrix for hypothesis 1 (Difference between M and average of TB)
contrast_matrix1 <- matrix(c(2, -1, -1), nrow=3, byrow=TRUE)
rownames(contrast_matrix1) <- c("M", "T", "B")
colnames(contrast_matrix1) <- c("M_vs_TB")

# Making a contrast matrix for hypothesis 2 (Difference between T and B)
contrast_matrix2 <- matrix(c(0, -1, 1), nrow=3, byrow=TRUE)
rownames(contrast_matrix2) <- c("M", "T", "B")
colnames(contrast_matrix2) <- c("T_vs_B")

# Combining into one contrast matrix 
contrast_matrix <- cbind(contrast_matrix1, contrast_matrix2)

# Mixture
contrast_matrix_mixture <- matrix(c(
  1,  1,    # M: gets positive coefficients for both contrasts
  -1,  0,    # B: compared to M in first contrast  
  0, -1     # T: compared to M in second contrast
), nrow = 3, byrow = TRUE)

colnames(contrast_matrix_mixture) <- c("M_vs_B", "M_vs_T")
rownames(contrast_matrix_mixture) <- c("M", "B", "T")



