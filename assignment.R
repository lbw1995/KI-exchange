# Task 1 - Literature
  
## 1. Read the research article of the hands-on working group you are assigned to (see file "Student Groups.pdf" in shared folder General course material).
  
## 2. Answer the following questions
  
### a. What is the medically relevant insight from the article?
  
Answer: In this study, they identified the cells with the characteristics of lacrimal gland primordium, which appeared in the two-dimensional eye-like organs cultured from human pluripotent stem cells. SEAM region 3 containing epithelioid cells on the surface of the eye produces precursors of lacrimal glands. They used flow cytometry to analyze the cells obtained from the isolated lacrimal gland clusters, and used previously determined antibody combinations (anti-ITGB4, anti-SSEA-4, and anti-CD200) to identify the eye surface epithelial stem cells. The study found that the lacrimal gland-like cell cluster was included in the eye surface epithelial stem cell group in SEAM, and the cell group (CD200+− ITGB4SSEA-4) was identified as the origin of lacrimal gland-like organs produced by HIPS cells. The researchers screened the best medium for the formation of lacrimal gland like organs, namely EGF/Y-27632 medium, and formed lacrimal gland like organs on this medium. When separated by cell sorting and growing under specified conditions, these cells form a three-dimensional lacrimal gland-like tissue organ with ducts and tips, which can be realized by budding and branching. Clonal colony analysis showed that these organs originated from multipotential epithelial stem cells. Organs have obvious similarities with local lacrimal gland in terms of morphology, immune marker characteristics and gene expression patterns. When transplanted to the eyes of recipient rats, functional maturation will occur, forming lumen and producing tear film protein.
### b. Which genomics technology/ technologies were used?

Answer: Gene expression analysis

## 3. Further related research questions

### a. List and explain at least three questions/ hypotheses you can think of that extend the analysis presented in the paper.

Sorry, I am a computer science background student and I can not list questions about this article.

### b. [Optional] Devise a computational analysis strategy for (some of) the listed questions under 3a.

Sorry, I can not devise any computational analysis strategy for it.

# Task 4 - R basic operations

## 1. What is the square root of 10?


cat(sqrt(10))


## 2. What is the logarithm of 32 to the base 2?

cat(log2(32))


## 3. What is the sum of the numbers from 1 to 1000?


cat(sum(1:1000))


## 4. What is the sum of all even numbers from 2 to 1000?

sum <- 0 

for (i in seq(2, 1000, by=2)) { 
  sum <- sum + i  
}

cat(sum_even)


## 5. How many pairwise comparisons are there for 100 genes?

sum_even <- sum(seq(from = 2, to = 1000, by = 2))
cat(sum_even)

## 6. And how many ways to arrange 100 genes in triples?

n <- 100
r <- 3

num_triplets <- factorial(n) / (factorial(r) * factorial(n - r))

print(num_triplets)


# Task 5 - Using R example datasets

## 1. Use the R internal CO2 dataset ("data(CO2)").

data(CO2)

## 2. Describe briefly the content of the CO2 dataset using the help function.

A: The CO2 data set is a time series object, including the atmospheric CO2 concentration data collected by the Monaloa Observatory in Hawaii from 1959 to 1997.

The dataset has two variables:
  
CO2: concentration of carbon dioxide in the atmosphere.

Plant: plant type used for experiment.

There are 468 observations in the data set.

## 3. What is the average and median CO2 uptake of the plants from Quebec and Mississippi?

install.packages("dplyr")
library(dplyr)
CO2 %>%
  # select only Quebec and Mississippi
  filter(Type %in% c("Quebec", "Mississippi")) %>% 
  # group by Type
  group_by(Type) %>%
  # compute mean and median
  summarize(mean_uptake = mean(uptake), median_uptake = median(uptake))  


## 4. [Optional] In the "airway" example data from Bioconductor, how many genes are expressed in each sample? How many genes are not expressed in any sample?

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pasilla", version = "3.16")
if (!require("airway"))
  BiocManager::install("airway")

library(pasilla)

data(airway)

airway_mat <- as.matrix(assay(airway))

airway_mat <- na.omit(airway_mat)

expressed_genes <- sum(rowSums(airway_mat) > 0)
not_expressed_genes <- sum(rowSums(airway_mat) == 0)

cat("Number of expressed genes:", expressed_genes, "\n")
cat("Number of not expressed genes:", not_expressed_genes, "\n")


# Task 6 - R Functions

## 1. Write a function that calculates the ratio of the mean and the median of a given vector.


mean_median_ratio <- function(vector) {
  mean <- mean(vector)
  median <- median(vector)

  ratio <- mean / median

  return(ratio)
}


## 2. Write a function that ignores the lowest and the highest value from a given vector and calculate the mean.

trimmed_mean <- function(vector) {
  trimmed_vector <- vector[-c(which.min(vector), which.max(vector))]
  mean(trimmed_vector)
}


## 3. Read about piping from here:<https://r4ds.had.co.nz/pipes.html#pipes> (you don't have to learn everything, a basic understanding of the usage is enough). Write a short (max. 300 characters, no spaces) explanation of why, how, and when not to use pipes.

The pipes in R are used to link multiple operations together, making the code more readable and efficient.
Pipelines allow data to flow from one operation to the next, thus reducing the need for intermediate variables.
However, when the pipeline becomes too complex or nested too deep, the pipeline may be difficult to read.
In addition, some operations may not be applicable to pipelines, such as functions that require multiple parameters or functions that need to group or sort data in a specific way.
Therefore, pipelines should be used wisely, not at the expense of code readability or functionality.

## 4. Familiarize yourself with the apply-family of functions (apply, lapply, sapply etc.) <http://uc-r.github.io/apply_family> Write a short explanation (max. 300 characters, no spaces) of why they could be useful in your work.

The application function series (apply, apply, apply, etc.) in R are very useful in my work because they allow efficient and simplified operations on data in arrays and lists.
These functions provide a simpler and more concise way to apply functions to subsets of data sets or across multiple data sets, thus reducing the number of duplicate code.
The apply-family function also allows you to return output in various formats, such as lists, vectors, or matrices, according to your analysis needs.

# Task 7 - Basic visualization with R

## 1. Compare the distributions of the body heights of the two species from the 'magic_guys.csv' dataset graphically

### a. using the basic 'hist' function as well as 'ggplot' and 'geom_histogram' functions from the ggplot2 package. Optimize the plots for example by trying several different 'breaks'. Note that ggplot2-based functions give you many more options for changing the visualization parameters, try some of them.

magic_guys <- read.csv("/home/bowen/下载/magic_guys.csv")
library(ggplot2)
body_height <- magic_guys$weight / magic_guys$length^2


species1 <- subset(body_height, magic_guys$species == "jedi")
species2 <- subset(body_height, magic_guys$species == "sith")

type(species1)
hist(species1, breaks = 20, col = "blue", main = "Species 1")
hist(species2, breaks = 20, col = "red", add = TRUE)

# Use ggplot2 and geom_histogram() to plot histograms of the body height data for both species
# Load the necessary library
library(ggplot2)

# Create a new data frame with the body heights of both species
new_df <- data.frame(Body_Height = c(species1, species2),
                     Species = factor(rep(c("Species 1", "Species 2"), c(length(species1), length(species2)))))

# Create the histogram plot
ggplot(new_df, aes(x=Body_Height, fill=Species)) +
  geom_histogram(alpha=0.5, position="identity", bins=30) +
  labs(title="Body Height Distribution by Species", x="Body Height", y="Count") +
  scale_fill_manual(values=c("red", "blue"))

```

### b. Do the same comparison as in a. but with boxplots. If you want to use the ggplot2-package, use the functions 'ggplot' and 'geom_boxplot'.

ggplot(new_df, aes(x=Species, y=Body_Height, fill=Species)) +
  geom_boxplot(alpha=0.5, position="dodge") +
  labs(title="Body Height Distribution by Species", x="Species", y="Body Height") +
  scale_fill_manual(values=c("red", "blue"))

```

### c. Save the plots with the 'png', 'pdf', and 'svg' formats. In which situation would you use which file format?



ggsave("temp_save.png", plot = last_plot(), width = 6, height = 4, dpi = 300)
ggsave("temp_save.pdf", plot = last_plot(), width = 6, height = 4, dpi = 300)
# install.packages("svglite") first
library(svglite)
ggsave("temp_save.svg", plot = last_plot(), width = 6, height = 4, dpi = 300, device = "svg")

If we need a transparent, lossless, Web-compatible format, we may choose PNG.
If we need a high-quality vector format that supports CMYK colors and can be embedded in documents, we may choose PDF.
If we need a scalable, network-compatible vector format, we may choose SVG.
## 2.
Load the gene expression data matrix from the 'microarray_data.tab' dataset provided in the shared folder, it is a big tabular separated matrix.
### a.
How big is the matrix in terms of rows and columns?

data <- read.table("/home/bowen/下载/microarray_data.tab", header = TRUE, sep = "\t")
dim(data)



### b.
Count the missing values per gene and visualize this result.


missing_values <- apply(data, 1, function(x) sum(is.na(x)))
missing_values
# Visualize the result
hist(missing_values, main = "Missing Values per Gene", xlab = "Number of Missing Values")


### c. Find the genes for which there are more than X% (X=10%, 20%, 50%) missing values.


# Calculate the percentage of missing values for each gene
percent_missing <- apply(is.na(data), 1, mean) * 100

X_values <- c(10, 20, 50)
for (X in X_values) {
  missing_genes <- rownames(data)[percent_missing > X]
  cat("Genes with more than", X, "% missing values:\n")
  print(missing_genes)
}

### d. Replace the missing values by the average expression value for the particular gene. (Note: Imputing data has to be used with caution!)

# Replace missing values with the average expression value for the particular gene
for (gene in colnames(data)) {
  gene_values <- data[,gene]
  missing_indices <- is.na(gene_values)
  if (any(missing_indices)) {
    avg_value <- mean(gene_values, na.rm=TRUE)
    gene_values[missing_indices] <- avg_value
    data[,gene] <- gene_values
  }
}
head(data)

# Task 8

## 1. Install the Tidybiology package, which includes the data 'chromosome' and 'proteins' devtools::install_github("hirscheylab/tidybiology")

### a. Extract summary statistics (mean, median and maximum) for the following variables from the 'chromosome' data: variations, protein coding genes, and miRNAs. Utilize the tidyverse functions to make this as simply as possible.

library(ggplot2)
data(chromosome)
library(tidybiology)
# Extract summary statistics for variations, protein coding genes, and miRNAs
chromosome %>%
  summarize(
    mean_variations = mean(variations),
    median_variations = median(variations),
    max_variations = max(variations),
    mean_protein_coding_genes = mean(protein_codinggenes),
    median_protein_coding_genes = median(protein_codinggenes),
    max_protein_coding_genes = max(protein_codinggenes),
    mean_miRNAs = mean(mi_rna),
    median_miRNAs = median(mi_rna),
    max_miRNAs = max(mi_rna)
  )


### b. How does the chromosome size distribute? Plot a graph that helps to visualize this by using ggplot2 package functions.

# Create a histogram of chromosome sizes
ggplot(chromosome, aes(x = basepairs)) +
  geom_histogram(binwidth = 10000000, fill = "blue", color = "white") +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "Distribution of Chromosome Sizes",
       x = "Size (base pairs)", y = "Frequency")

### c. Does the number of protein coding genes or miRNAs correlate with the length of the chromosome? Make two separate plots to visualize these relationships.

# Create a scatterplot for protein coding genes
ggplot(chromosome, aes(x = basepairs, y = protein_codinggenes)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "Correlation between Chromosome Size and Protein Coding Genes",
       x = "Size (base pairs)", y = "Protein Coding Genes")

# Create a scatterplot for miRNAs
ggplot(chromosome, aes(x = basepairs, y = mi_rna)) +
  geom_point(color = "green") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "Correlation between Chromosome Size and miRNAs",
       x = "Size (base pairs)", y = "miRNAs")

### d. Calculate the same summary statistics for the 'proteins' data variables length and mass. Create a meaningful visualization of the relationship between these two variables by utilizing the ggplot2 package functions. Play with the colors, theme- and other visualization parameters to create a plot that pleases you.

# To calculate the summary statistics for the 'proteins' data variables length and mass, we can use the summarize() function from the dplyr package.
proteins_summary <- proteins %>%
  summarize(mean_length = mean(length),
            median_length = median(length),
            max_length = max(length),
            mean_mass = mean(mass),
            median_mass = median(mass),
            max_mass = max(mass))

proteins_summary

# Use the ggplot() function to create the plot, and add geom_point() to add the points.
ggplot(proteins, aes(x = length, y = mass)) +
  geom_point() +
  labs(x = "Length", y = "Mass") +
  theme_classic()
