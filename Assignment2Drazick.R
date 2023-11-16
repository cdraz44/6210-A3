##load packages to be used in this analysis

library(tidyverse)
library(dbplyr)
library(ape)
library(muscle)
library(DECIPHER)
library(rentrez)
library(Biostrings)
library(ggplot2)
library(ggtext)
library(cluster)
library(factoextra)
library(viridis)

# Listing global variables to be searched in NCBI Genback nucleotide database. In this case, NOTCH3 and BRCA1 genes from the Cetacea family were searched against a certain length range to isolate gene sequences rather than genomes.
family = "Cetacea"
gene1 = "NOTCH3"
gene1_min = 6000
gene1_max = 9000
gene2 = "BRCA1"
gene2_min = 5000
gene2_max = 8000

# Function to search and fetch genes using the rentrez package, and create a dataframe containing gene sequences. Also creates global variables for quality control purposes.
create_dfGene <- function(family_name, gene, gene_min, gene_max) {
  # Search NCBI nucleotide database for genes of a given length
  gene_search <- entrez_search(db = "nuccore", term = paste0(family_name,"[ORGN] AND ", gene, "[Gene] AND ", gene_min, ":", gene_max, "[SLEN]"), retmax = 1000)
  
  # Fetch fasta files corresponding to search parameters
  gene_fetch <- entrez_fetch(db = "nuccore", id = gene_search$ids, rettype = "fasta")

  # write to file, separate and remove Ns from sequence data
  write(gene_fetch, paste0(gene, "_fetch.fasta"), sep = "\n")
  # Rewriting to a DNA string set from the .fasta file without any Ns in the sequences
  gene_string <- readDNAStringSet(paste0(gene, "_fetch.fasta"))
  
  # Creating a data frame from our stringset taking values(names) from stringset and plugging them into data frame
  dfGene <- data.frame(title = names(gene_string), sequence = paste(gene_string))
  
  # Creating the column names for the dataframe
  gene_title <- paste0(gene, "_Title")
  gene_sequence <- paste0(gene, "_Sequence")
  names(dfGene)[1:2] <- c(gene_title, gene_sequence)
  
  ## Creating global variables for each gene for quality control purposes
  
  # Entrez fetch variable
  assign(paste0(gene, "_fetch"), gene_fetch, parent.frame())
  # Gene DNAStringSet Object using Biostrings package
  assign(paste0(gene, "_string"), gene_string, parent.frame())
  # Create Gene Data frame
  assign(paste0("df", gene), dfGene, parent.frame())
}

# Fetch NOTCH3 sequences and create data frame
create_dfGene(family_name = family, gene = gene1,  gene_min = gene1_min, gene_max = gene1_max)
# Fetch BRCA1 sequences and create data frame
create_dfGene(family_name = family, gene = gene2,  gene_min = gene2_min, gene_max = gene2_max)

##Check class to ensure it is a character vector

class(NOTCH3_fetch)

class(BRCA1_fetch)

##Having a look at the data to ensure we have the right sequence data

head(BRCA1_fetch)

head(NOTCH3_fetch)

##Checking data, ensuring we have proper class and viewing our data
class(NOTCH3_string)
head(names(NOTCH3_string))

class(BRCA1_string)
head(names(BRCA1_string))

##This is some name editing to clean up sequence names and remove unwanted words from species names, and rearranging data frame columns

clean_df <- function(dfGene, gene_name) {
  
  dfGene$Species_Name <- word(dfGene[,1], 3L, 4L)
  colnames <- c(paste0(gene_name, "_Title"), "Species_Name", paste0(gene_name, "_Sequence"))
  return(dfGene[, colnames])
  
}

dfNOTCH3 <- clean_df(dfNOTCH3, gene1)
dfBRCA1 <- clean_df(dfBRCA1, gene2)

##Checking dimensions of the dataframe as quality control

dim(dfBRCA1)

dim(dfNOTCH3)

##Checking names of dataframe as quality control
names(dfBRCA1)

names(dfNOTCH3)

##Getting a summary of sequence lengths to ensure sequence length is appropriate

summary(nchar(dfBRCA1$BRCA1_Sequence))

summary(nchar(dfNOTCH3$NOTCH3_Sequence))

##coding for a histogram displaying the distribution of sequence lengths, to ensure quality control and identify any errors or outliers

create_plot <- function(dfGene, gene_name, family_name, gene_min, gene_max) {
  
  plot_title <- paste0("Frequency of Sequence Length of ", gene_name, " in ", family_name)
  
  hist_plot <- ggplot(data = dfGene,
                      mapping = aes(x = nchar(dfGene[,3]))) + 
    geom_histogram(breaks= seq(gene_min,gene_max, by = 75), 
                   col="blue", aes(fill=..count..)) + 
    scale_fill_distiller(palette= "Spectral") + 
    labs(title = plot_title, x = "Sequence Length", y = "Number of Species")  +
    theme_bw()
  
  return(hist_plot)
}

NOTCH3PLOT <- create_plot(dfNOTCH3, gene1, family, gene1_min, gene1_max)

BRCA1PLOT <- create_plot(dfBRCA1, gene2, family, gene2_min, gene2_max)

##Plotting the plot to visualize

plot(NOTCH3PLOT)

plot(BRCA1PLOT)

# Alignment Function
run_alignment <- function(dfGene, gene) {
  #Create new column name
  colname <- paste0(gene, "_Sequence2")
  #Add DNAStringSet sequences to new column
  sequences <- dfGene[,3]
  dfGene[[colname]] <- DNAStringSet(sequences)
  #Assign species name to each sequence for alignment
  names(dfGene[[colname]]) <- dfGene$Species_Name
  #Align using MUSCLE
  dfGene.alignment <- DNAStringSet(muscle::muscle(dfGene[[colname]]), use.names = TRUE)
  # Update data frame global variable
  assign(paste0("df", gene), dfGene, parent.frame())
  # Return sequence alignment
  return(dfGene.alignment)
}

##Create alignments for NOTCH3 and BRCA1
dfNOTCH3.alignment <- run_alignment(dfNOTCH3, gene1)
dfBRCA1.alignment <- run_alignment(dfBRCA1, gene2)

##checking class -> string set
class(dfNOTCH3$NOTCH3_Sequence2)

class(dfBRCA1$BRCA1_Sequence2)

##Double checking that the names have transferred in and checking they are cleaned up
names(dfNOTCH3$NOTCH3_Sequence2)

names(dfBRCA1$BRCA1_Sequence2)

##Checking that alignment worked, with species names 

dfBRCA1.alignment 

dfNOTCH3.alignment

##Opening an alignment in the browser to double check alignment quality

BrowseSeqs(dfBRCA1.alignment)

BrowseSeqs(dfNOTCH3.alignment)

##Changing alignment class to DNA bin to be used in a distance matrix for subsequent cluster analysis

dnaBin.NOTCH3 <- as.DNAbin(dfNOTCH3.alignment)

dnaBin.BRCA1 <- as.DNAbin(dfBRCA1.alignment)

##Checking that it worked and is now DNA bin class

class(dnaBin.NOTCH3)

class(dnaBin.BRCA1)

##Setting values for cluster analysis

missing.data <- 0.01

length.var <- 10

chosen.model <- "TN93"

clustering.threshold <- 0.03

clustering.method <- "UPGMA"

##Creating distance matrixes

distanceMatrixBRCA <- dist.dna(dnaBin.BRCA1, 
                               model = chosen.model, 
                               as.matrix = TRUE, 
                               pairwise.deletion = TRUE)

distanceMatrixNOTCH3 <- dist.dna(dnaBin.NOTCH3, 
                                 model = chosen.model, 
                                 as.matrix = TRUE, 
                                 pairwise.deletion = TRUE)

##Ensuring that it worked 

head(distanceMatrixBRCA)

head(distanceMatrixNOTCH3)


##!!please note using IdClusters here instead of TreeLine!! -Clustering my data using previously chosen cluster values
# clusters.NOTCH3 <- DECIPHER::IdClusters(myDistMatrix = distanceMatrixNOTCH3, 
#                                         method = clustering.method,
#                                         cutoff = clustering.threshold,
#                                         showPlot = TRUE,
#                                         type = "both",
#                                         verbose= TRUE)
# 
# clusters.BRCA1 <- DECIPHER::IdClusters(myDistMatrix = distanceMatrixBRCA,
#                                        method = clustering.method,
#                                        cutoff = clustering.threshold,
#                                        type = "both",
#                                        showPlot = TRUE,
#                                        verbose = TRUE)
clusters.NOTCH3 <- DECIPHER::TreeLine(myDistMatrix = distanceMatrixNOTCH3, 
                                        method = clustering.method,
                                        cutoff = clustering.threshold,
                                        showPlot = TRUE,
                                        type = "both",
                                        verbose= TRUE)

clusters.BRCA1 <- DECIPHER::TreeLine(myDistMatrix = distanceMatrixBRCA,
                                       method = clustering.method,
                                       cutoff = clustering.threshold,
                                       type = "both",
                                       showPlot = TRUE,
                                       verbose = TRUE)

##Checking class, viewing and checking length of clusters for quality control

class(clusters.BRCA1)

clusters.BRCA1

length(clusters.BRCA1)


class(clusters.NOTCH3)

clusters.NOTCH3

length(clusters.NOTCH3)

##numbering species names making each a unique name to use as a unique data point on our plot
rownames(distanceMatrixNOTCH3) = make.names(rownames(distanceMatrixNOTCH3), unique = TRUE)

rownames(distanceMatrixBRCA) = make.names(rownames(distanceMatrixBRCA), unique = TRUE)

##Setting the seed for reproduction
set.seed(123)

##Finding kmeans values and indicating number of clusters.(I decided on these numbers for k-means clusters through trial and error and looking at kmeans results to decide which values to use for plotting)

kmNOTCH3 <- kmeans(distanceMatrixNOTCH3,3)

kmBRCA1 <- kmeans(distanceMatrixBRCA, 6)

##Viewing kmeans outcomes and evaluating whether this is the right number of clusters

kmNOTCH3

kmBRCA1

##Creating kmeans cluster plot in ggplot style, colouring and tweaking default parameters to match the goal of these plots
kmeans_cluster <- function(kmGene, distanceMatrix, gene_name, labelsize = 5, viridis_option = "H") {
  figure_title <- paste("K-means Analysis of", gene_name)
  gene_cluster <- fviz_cluster(kmGene,
                               data = distanceMatrix, 
                               geom = c("point","text"), 
                               repel = TRUE, 
                               show.clust.cent = TRUE, 
                               ellipse.type = "convex", 
                               labelsize = 5) +
    labs(title= figure_title) + 
    theme(panel.background = element_rect(fill = "darkgrey")) + 
    scale_color_viridis_d(option = viridis_option)
  return(gene_cluster)
}

NOTCH3_cluster <- kmeans_cluster(kmNOTCH3, distanceMatrixNOTCH3, gene1, 5, "H")

BRCA1_cluster <- kmeans_cluster(kmBRCA1, distanceMatrixBRCA, gene2, 8, "A")

##Plotting and viewing of plots
plot(BRCA1_cluster)

plot(NOTCH3_cluster)

##To ensure we have the right amount of clusters for our k-means analysis, I am creating a plot to tell us the optimal number of clusters for each using silhouette index calculations
OptNOTCH3 <- fviz_nbclust(distanceMatrixNOTCH3, FUNcluster = kmeans, method = "silhouette", linecolor = "lightpink2")

OptBRCA1 <- fviz_nbclust(distanceMatrixBRCA, FUNcluster = kmeans, method = "silhouette", linecolor = "lightblue2")

##Plotting optimal cluster plot
plot(OptNOTCH3)

plot(OptBRCA1)

##Changing the number of K-mean clusters

kmNOTCH3 <- kmeans(distanceMatrixNOTCH3,5)

kmBRCA1 <- kmeans(distanceMatrixBRCA, 4)

##Re-making the graphs using this new parameter

NOTCH3_cluster <- fviz_cluster(kmNOTCH3, data = distanceMatrixNOTCH3, geom = c("point","text"), 
                               repel = TRUE, show.clust.cent = TRUE, ellipse.type = "convex", labelsize = 5) + 
  labs(title= "K-means Analysis of NOTCH3") + theme(panel.background = element_rect(fill = "darkgrey")) + scale_color_viridis_d(option = "H")

BRCA1_cluster <- fviz_cluster(kmBRCA1, data = distanceMatrixBRCA, geom = c("point","text"), 
                              repel = TRUE, ellipse.type = "convex", labelsize = 8) + 
  labs(title= "K-means Analysis of BRCA1") + theme(panel.background = element_rect(fill = "darkgrey"))   + scale_color_viridis_d(option = "A") 

##Plotting and viewing new plot with updated kmeans parameter

plot(NOTCH3_cluster)

plot(BRCA1_cluster)





