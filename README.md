# R-based data analysis (gene expression & tree growth) and comparative sequence analysis

## Project summary
This repository contains two linked analyses:

- **Part 1**: Data import, wrangling, summary statistics, visualization and simple hypothesis testing using:
  - `gene_expression.tsv` (RNA-seq counts)
  - `growth_data.csv` (tree circumference, control vs treatment over 20 years)

- **Part 2**: Comparative sequence analysis of coding DNA and translated proteins for:
  - *Oxalobacter formigenes* OXCC13 (assigned organism)
  - *Escherichia coli* K-12 MG1655 (reference)

All scripts are written in **R**, designed for reproducibility, and generate outputs (tables and plots) in the `outputs/` folder.

## Repository structure

# Repository structure
Assessment/

├── data/  

├── output

├── script.R

└── README.md


## Requirements and installation

### Software
- R ≥ 4.0 (recommended)
- RStudio (recommended for development and RMarkdown rendering)

### R packages
Install once:

install.packages(c("dplyr", "ggplot2", "tidyr", "seqinr", "R.utils"))

## How to run the analysis (step-by-step)
1.	Clone the repository and set the working directory to the repository root.
2.	Confirm that input data files are present in data/. For Part 2, the Ensembl-provided .fa.gz files may be downloaded via browser and placed in data/ (recommended if download.file() fails).
3.	Install required packages (see above).
4.	Run
source("script.R")


## Part 1 — Gene expression & growth analysis
Summary of approach
-	Data import verification is performed prior to any modification.
-	Row-wise operations (e.g., means) are added as new columns to preserve raw counts.
-	Visualizations use ggplot2 and are saved to outputs/part1/ for report inclusion.
-	Statistical tests check preconditions (non-missing data, adequate sample counts). If t-test requirements are not met, results are reported and alternative visualizations are provided.


## Q1–Q3: Data import, mean calculation, top genes
Objective: import gene_expression.tsv, make gene identifiers row names, compute gene-wise mean expression, and list the 10 genes with highest mean expression.

Method:
-	read.table("data/gene_expression.tsv", header = TRUE, sep = "\t", row.names = 1) imports gene counts and sets gene IDs as the data frame row names. Using row names simplifies row-based aggregation and ensures that gene identifiers are preserved when adding new columns.
-	rowMeans() calculates the mean expression for each gene across all samples. Adding the mean as a separate column allows subsequent sorting and filtering without altering the original counts.
-	Sorting is performed with order(-gene_data$mean) to rank genes by mean descending; the top 10 are selected using head(..., 10).

<img width="975" height="187" alt="image" src="https://github.com/user-attachments/assets/0fdfbff4-c329-4e8e-b974-69cdd63015d3" />


Top 10 genes by mean:

<img width="975" height="329" alt="image" src="https://github.com/user-attachments/assets/1c0434dd-4763-47de-9a89-204615888c41" />


The top expressed genes include mitochondrial genes and canonical high-expression genes. The mean-based ranking reduces influence from single-sample extremes and highlights consistently highly expressed genes across samples.

## Q4–Q6: Growth dataset columns, mean & SD, boxplot
Objective: Understand the growth data layout, compute site-specific summary statistics at start (2005) and end (2020), and visualize distributions.
Key steps:
-	Load with read.csv("data/growth_data.csv").
-	Inspect columns via colnames() to ensure the expected fields are present.
-	Compute means and standard deviations grouped by Site.
-	Reshape to long format and plot boxplots with ggplot2::geom_boxplot(). Y-axis ticks have been adjusted to show 10 cm increments for clarity.

<img width="975" height="589" alt="image" src="https://github.com/user-attachments/assets/9e1e395a-5874-4cba-be29-6aa18a727271" />

 
Interpretation guidance:

Examine differences in median values and spread (IQR). Larger median increase at treatment site suggests treatment-related growth advantages; verify with Q7 and Q8.

## Q7–Q10: Mean growth over last 10 years & t-test
Objective: calculate per-tree growth 2010→2020, report site means, attempt statistical test (t-test), and provide alternative visualization when t-test cannot be performed due to insufficient data.
Method:
-	Growth per tree: Growth_10yr = Circumf_2020_cm - Circumf_2010_cm.
-	Group summary: group_by(Site) %>% summarise(Mean_Growth = mean(Growth_10yr, na.rm=TRUE)).
-	T-test:
  - Remove NAs, check group counts.
  -	If both groups have ≥ 2 observations, run t.test(); otherwise, print a clear message and produce a boxplot of Growth_10yr by site for visual comparison.

Note:

The t-test could not be performed because one or both sites had insufficient non-missing observations for Growth_10yr. In such cases, the code conservatively reports the limitation and includes boxplot visual assessment.

<img width="975" height="559" alt="image" src="https://github.com/user-attachments/assets/b76841c2-7650-4885-ae4b-ee6066fcd61a" />

 
# Part 2 — Comparative sequence diversity analysis

Overview of approach
-	CDS FASTA files were loaded with seqinr::read.fasta() (DNA sequences). If .fa.gz files were downloaded, R.utils::gunzip() was used to decompress in R while keeping .gz.
-	All sequence-based operations used the raw CDS sequences; translation to protein sequences was performed via seqinr::translate() when protein analyses were required.
-	Codon usage was calculated using seqinr::uco() (codon counts) and RSCU values were computed by grouping by amino acid and normalizing observed codon counts against equal-synonym expectation.
-	K-mer analysis for proteins used sliding-window counting across all translated proteins; counts were normalized by total k-mers to obtain relative frequencies. Over/under representation was assessed by ranking and by comparison of species frequencies (observed vs. reference).

## Q1: Number of coding sequences
Method: length(read.fasta(...)) was used to count CDS.
Observed results:

<img width="975" height="237" alt="image" src="https://github.com/user-attachments/assets/9974afcd-d327-4fcc-9e17-2f26848931cd" />


Interpretation: E. coli has roughly twice the number of coding sequences, a finding consistent with larger genome size and broader metabolic repertoire.

## Q2: Total coding DNA length
Method summary:
-	For each CDS list: compute per-sequence length with nchar(getSequence(..., as.string=TRUE)).
-	Sum lengths with sum() to obtain total coding DNA length in base pairs.
 Interpretation guideline: larger total coding length indicates more protein-coding content and greater coding capacity.

## Q3: CDS length distribution, mean and median
Method
-	Per-CDS lengths computed with nchar() and summarized using mean() and median().
-	Distributions visualized with a ggplot2 boxplot. A histogram may also be used for a more granular view.

| Organism                      | Mean_Length | Median_Length |
|-------------------------------|-------------|---------------|
| Escherichia coli              | 938.5534    |  831          | 
| Oxalobacter formigenes        | 969.2820    |  801          | 


 <img width="975" height="600" alt="image" src="https://github.com/user-attachments/assets/a3242c3d-8c78-448b-a1c4-339cc267a0be" />


Interpretation: Although Oxalobacter has fewer CDS, the mean length is slightly higher, indicating fewer but somewhat longer coding sequences relative to E. coli. The median differences indicate distributional differences (e.g., more short genes in E. coli).

## Q4: Nucleotide and amino acid frequence
Method:
-	Concatenate all CDS sequences for a species (unlisting the FASTA object) then compute nucleotide counts with table() or alphabetFrequency() and normalize to relative frequencies by dividing by the total base count.
-	Translate CDS to protein using translate() and flatten the protein list to count amino acid occurrences via table().
-	Produce bar charts for nucleotide composition (A, T, G, C) and for the 20 amino acids. Use ggplot2 with position = "dodge" to compare species side-by-side on the same plot.

<img width="975" height="589" alt="image" src="https://github.com/user-attachments/assets/2f77168b-fa6a-4205-a50e-54b2901539d2" />

 
Interpretation guidance:
-	Differences in GC content are immediately visible in nucleotide bars (higher G/C proportion in one species vs. A/T in the other).
-	Amino acid frequencies reflect protein-coding preferences and may relate to codon bias and genome GC content.

## Q5: Codon usage & RSCU — method and top-codon results
Computation method:
1.	Use uco() (from seqinr) on concatenated CDS to obtain codon counts across all CDS for each organism.
2.	For each amino acid, compute expected codon count under equal synonymous usage: Expected = sum(counts_for_AA) / n_codons_for_AA.
3.	Compute RSCU: RSCU = Observed_count / Expected for each codon.
4.	Rank codons by RSCU to identify overrepresented codons.
5.	Save full RSCU tables to CSV (outputs/part2/codon_rscu_ecoli.csv, codon_rscu_oxalo.csv). 

Interpretation guidance:
-	A tendency to favor G/C-ending codons in E. coli indicates higher genomic GC bias whereas Oxalobacter tends to favor A/T ending codons.
-	These biases affect translational efficiency and may be adaptive.


## Q6: Protein k-mer (3–5 mers) analysis 
Method:
1.	Translate all CDS to protein sequences with translate().
2.	For each protein, extract overlapping k-mers (length k = 3, 4, 5) using a sliding window:
o	For a protein of length L, total possible k-mers per protein = L - k + 1 (skip sequences shorter than k).
3.	Aggregate k-mer counts across all proteins for each organism.
4.	Convert raw counts to relative frequencies by dividing by total number of k-mers of that length.
5.	Identify the top 10 over-represented (highest relative frequency) and top 10 under-represented (lowest non-zero relative frequency) k-mers in Oxalobacter.
6.	For each of these k-mers, report the relative frequency in E. coli for direct comparison (saved to CSV).
7.	Visualize the comparison with side-by-side bar plots or dot plot.

Interpretation guidance:
- Over-represented k-mers often correspond to common protein motifs, signal sequences or structural patterns.
-	Under-represented k-mers may avoid unfavorable physico-chemical properties or motifs that incur deleterious effects.
-	Differences between species reflect divergent functional repertoire and evolutionary pressures.


## End of README. ##
