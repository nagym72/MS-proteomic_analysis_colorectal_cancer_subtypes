#Part 2 Statistical analysis of protein counts in 3 subtypes of colorectal tumors 

install.packages("gplots")
library("gplots")
library("ibb")
#data_path <- "./example-3groups.txt"

df <- read.csv(data_path, sep="\t")
df

#function used to normalize both columns(total proteom in the cell) and rows(same protein within patients)

normalize_func <- function(data){

  #First normalize columns:
  #Normalization to account for total count of proteins
  
  #1) Total counts of spectra per sample (col)
  total_counts <- colSums(data)
  #2) Compute mean over all Sample totals 
  avg_count_over_cols <- mean(total_counts)
  #3) Generate multipliers for each sample depending on which col they belong
  ratios <- avg_count_over_cols/total_counts

  #double transpose in order to multiply the vector (ratios) with transposed matrix and then back transpose to get original matrixshape
  norm_cols <- t((t(data)*(ratios)))
  
  
 #ColSums are identical after normalization. sum of ColSums before and after are both 153272(identical)
  
  
  #Second normalization of rows:
  
  #mean = 0 sd = 1 -> z-score  (x-mean) / sd 
  #mean of row -> mean for protein X between the 8 samples
  #Margin = 1 for rows
  #transpose afterwards to retrieve a 1786 x 8 matrix instead of 8 x 1786
  norm_col_and_row <- t(apply(norm_cols, 1, FUN=function(x)((x-mean(x))/sd(x))))
  
  return(norm_col_and_row)
}


normalized_matrix <- normalize_func(df)

#?heatmap.2
heatmap.2 (t(normalized_matrix), col=topo.colors(50),
           key=TRUE, keysize=1,cexRow=1,labCol=NA, scale="none",
           symkey=FALSE, density.info="none", trace="none")


#?ibb.test
out <- ibb.test(x=df[,1:6], tx=colSums(df[,1:6]),c(rep("pre", 3), rep("post", 3)), n.threads=1)
#wrting results to outfile

#just for me to inspect
write.table(out, file="beta_binom_p_grouped_vales.tsv", sep="\t")  

#converting the list result from the ibb containing the required probabilities to a vector
results <- as.vector(out[[1]])
#filtering out the relevant proteins (only for col 1:6 because 7:8 are already established as subtype c) 
relevant_p_values <- normalized_matrix[results<0.01, 1:6]

#plot 2 with only those values that are p < 0.01 and only 6 samples (a1:3 and b1:3)
heatmap.2 (t(relevant_p_values), col=topo.colors(50),
           key=TRUE, keysize=1,cexRow=1,labCol=NA, scale="none",
           symkey=FALSE, density.info="none", trace="none")


