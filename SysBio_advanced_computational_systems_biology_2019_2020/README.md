# SysBio - Advanced Computational Systems Biology course

Here can be found the materials used in the Advanced Computational Systems
Biology course imparted by Julio Saez-Rodriguez's group.

The course covers both theoretical and computational systems biology.
Students learn about data analysis, dimensionality reduction, data modeling and
high-throughput data. This part mainly covers pre-processing of omics data and
some downstream analyses.

Prerequesites: install latest git, Rstudio and R (optionally Cytoscape)

**Slides [here](https://www.dropbox.com/sh/3leev4mgdd231li/AABTPGr1WI1eHiJXGy4tergCa?dl=0)**

**To be updated?**

## Day 1 - Motivation (SaezLab research) & Data types, resources, and tools overview
0. General talk of the research in the lab - as motivation of what they will learn

1. Different -omics- (they have done this already, just reminder)
    1. genomics (SNPs, WES, WGS)
    2. transcriptomics (microarrays, RNA-seq, mention scRNA)
    3. proteomics (antibody based, mass spec based)
    4. metabolomics (mass spec based)

2. Different resources for different -omics, and tools
     1. Resources: EGA, arrayexpress/Expression atlas, GEO, PRIDE, metabolights,...
     2. Tools:
        1. online tools to download/explore data
        2. R/Bioconductor, Python (mention biojupies)
        3. Mention good practices of coding - GitHub/repo etc.
     3. Identifiers etc - get the data I want to work with

3. Summary statistics on data
    1. count, median, mean, deviation,..
    2. Correlation
    3. Density plot, scatter-plots,.. (do not use pies!)
    4. ggplot?

#### Day 1 - practical
0. git and version control
1. Installing packages from Bioconductor
2. Explore resources - find data sets,...
3. Download the NCI-60 from GEO directly to R + metadata
4. Manipulate, BioMart - e.g. probe names to uniprot
5. Saving data to file
6. Retrieve pathway annotation from MSigDB (*)
7. Obtaining PPI network from OmniPath (*)
8. Mapping expression data to network in Cytoscape (*)

(*) Optional/extra, depending on time and students interest

## Day 2 - Exploratory Analysis and normalization
1. Visualization
    1. PCA, tSNE
    2. Clustering
    3. ....
2. Normalisation, scaling,  

#### Day 2 - practical
1. Installing packages from CRAN
2. Loading data from file (from previous day)
3. Subsetting data
4. Generating basic plots + basic stats
5. Heatmap + hierarhical clustering
6. Scatter plots + correlation
7. Dimensionality reduction (PCA and t-SNE)

## Day 3 - Basic statistical tests
1. Experimental design,
2. t-test, ANOVA, linear models, …
3. Multiple hypothesis corection
4. contrasts, differential exrpression Limma, edger, deseq (just basic ideas, not in depth, focus on limma)

#### Day 3 - practical
Building on the NCI-60 data, try on this data the different tests: t- tests/ANOVAs on genomic-drug response. Linear model/limma..

## Day 4 - Functional analysis
1. Basis enrichment, hypergeometric test, GSEA,
    1.2. Mention more advanced enrichment methods  (AREA), also PIANO as aggregated
2. Footprint methods: DOROTHEA, PROGENY and CARNIVAL
3. Online-tools (ENRICHR,DAVID, shiny apps like ours)

#### Day 4 - practical
Building on the NCI-60, run on the expression of the cell lines GSEA, PROGENY,
DOROTHEA, CARNIVAL.


## References
##### General reference on statistical analysis for omics data (free online)
> [Modern Statistics for Modern Biology. Susan Holmes, Wolfgang Huber](https://www.huber.embl.de/msmb/index.html)
##### Other good introductory books on statistics
> [Intuitive biostatistics](https://www.amazon.de/Intuitive-Biostatistics-Nonmathematical-Statistical-Thinking/dp/0190643560/ref=dp_ob_title_bk)
>[An introduction to statistical learning](http://www-bcf.usc.edu/~gareth/ISL/)
