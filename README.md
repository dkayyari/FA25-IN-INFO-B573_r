# B573_Assignment 2  
**Name**: Deeksha Kayyari  
**Programming Language**: R  
**Date**: 10/12/2025  

---

## **Description**

This project demonstrates the use of **fundamental Base R commands and flow control** through the analysis of a DNA sequence file.  
The program reads, processes, and analyzes the nucleotide sequence of `chr1_GL383518v1_alt.fa` to perform several bioinformatics operations including string manipulation, sequence reversal, base counting, and structured data summarization.  

The assignment emphasizes proficiency with **loops**, **conditionals**, **lists (dictionaries)**, and **data frames** — building toward a complete, well-documented R workflow.  

---

## **Required Files**
- `chr1_GL383518v1_alt.fa` – FASTA file containing the DNA sequence  
- `Assignment_R_Base_Deeksha_Kayyari.R` – R script implementing all analysis steps  

---

## **Required Packages**
- R version 4.4.0 or later  
- `stringr` (for base counting operations)  

You can install the package in R using:
```r
install.packages("stringr")
```

## **Steps for Execution**
1.Run in RStudio or R Console

2.Ensure both the R script and FASTA file are located in the same folder.

3.Open RStudio (or R GUI).

4.Set the working directory to the folder containing both files:


```
setwd("path_to_your_folder")
```
5.Run the script:

```
source("Assignment_R_Base_Deeksha_Kayyari.R")
```
6.Review the console output for results from each part.

 ## output summary

-Successfully read and cleaned the DNA sequence from chr1_GL383518v1_alt.fa.

-Displayed the 10th and 758th bases from the original sequence.

-Generated the reverse complement of the DNA sequence using Watson–Crick–Franklin base pairing (A↔T, C↔G).

-Displayed the 79th base and the 500th–800th bases of the reverse complement.

-Created a list (dictionary) showing the number of A, C, G, and T bases in each kilobase (1000 bases) of the sequence.

-Displayed the first five kilobase entries to verify correct list structure.

-Converted the nested list into a data frame containing columns for A, C, G, T, and the calculated Sum per kilobase.

-Verified that each kilobase’s total equals the expected 1000 bases, with exceptions explained (partial final segment or removed invalid bases).

-Provided an interpretation section explaining observed discrepancies in total counts and confirming the overall data integrity.

-Output displays are clearly labeled for each section (1a–4e) to match assignment requirements.

-Program concludes with a message indicating successful execution.




