Functional analysis using footprint approaches
================
Rosa Hernansaiz Ballesteros
22/11/2021

## Installing packages

The required packages to run the different steps of this analysis are in
the CRAN repository, Bioconductor or gitHub.

``` r
library(dplyr)
```

    ## Warning: package 'dplyr' was built under R version 4.1.3

``` r
library(ggplot2)

# source("assignPROGENyScores.r")
# source("generateTFList.r")
```

## Data description and loading

The dataset we are using for this sesion is a subset of the GEO entry
`GSE116436`. This is a GEO series of the whole data set. It comprehends
transcriptome profiling upon **60 cell lines using 15 anti-cancer
agents** with different time points (2, 6, 24 hrs) and drug
concentrations (incl. null concentration, i.e. control). A short
description and experimental overview are shown. For further details,
you can read the original article:

> Monks A et al. The NCI Transcriptional Pharmacodynamics Workbench: A
> Tool to Examine Dynamic Expression Profiling of Therapeutic Response
> in the NCI-60 Cell Line Panel. Cancer Res 2018 Dec
> 15;78(24):6807-6817. PMID:
> [30355619](https://www.ncbi.nlm.nih.gov/pubmed/30355619)

or you can also visit the GEO entry `GSE116436` at

> <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116436>

As this **dataset is large (7.5k samples)**, it is stratified in
subseries. Thus, we choose to work with one subserie: lapatinib
(`GSE116445`).

As sometimes Ensembl servers are down, we provide here directly the
filter matrices:

-   expression data matrix: data.frame of 13226 genes (gene symbol) vs
    485 conditions

-   metada matrix: data.frame containing extra information about the
    samples

``` r
metadata = readRDS("../data/NC60_metadata.rds") %>%
  dplyr::mutate_if(is.factor, as.character) %>% 
  dplyr::rename(tissue = "tissue:ch1") %>%
  tibble::tibble() 

expression = readRDS("../data/NC60_expression.rds")
```

## Getting ready for downstream analysis

First of all, we need to twick the `metadata` matrix. The information in
the `title` column is provided in an inconvinient format. So we are
going to get it in a more easy-access way.

``` r
#split title column, get it with the correct data type and remove unecesary information
metadata = metadata %>%
  tidyr::separate(col = title, sep = "_",
                  into = c("sample", "drug", "concentration_nM", "time_h")) %>%
  dplyr::mutate(concentration_nM = gsub("nM", "", concentration_nM)) %>%
  dplyr::mutate(time_h = gsub("h", "", time_h)) %>%
  dplyr::mutate(dplyr::across(c(concentration_nM, time_h), 
            as.numeric)) %>%
  dplyr::select (-c(platform_id))
```

It is usually a good idea to know the data we are working with before
starting any analysis…

``` r
# Create a palette

## List of unique tissues in our data set
tissues = unique(metadata$tissue)

## List of unique colors for each tissue
tissue_colors = RColorBrewer::brewer.pal(length(tissues), 'Set1')

## List of tissue-assigned colors for each sample
sample_colors = tissue_colors[as.integer(factor(metadata$tissue,
levels=tissues))]

barplot(table(metadata$tissue), cex.names=.6, col=tissue_colors)
```

And compute some exploratory statistics in order to get a preliminary
idea of our expression data.

``` r
gg_xplr = expression %>%
  tibble::rownames_to_column(var = "genes") %>%
  reshape2::melt(id.vars="genes", 
                 variable.name='geo_accession', 
                 value.name="exp") %>%
  merge.data.frame(., metadata, by = 'geo_accession')

# Density curves
ggplot(gg_xplr, aes(x=exp, colour=time_h)) + 
  geom_density(alpha=.2) +
  facet_grid(concentration_nM~tissue) + 
  xlim(0, 15) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
```

### Differential Expression of NCI-60

Our collaborators are working on drug reporpusing. They need help to
analyse their data. Specifically, they want to know if there are any
differencies between treated and untreated samples when using lapatinib
in different cell lines derived from Lung cancer.

Now that we know how our data looks like, we will perform a
*differential expression analysis* to try to answer this specific
question. Thus, we are going to select all cell lines derived from Lung
to check if the drug treatment has an effect on its expression.
Concretely, we will choose basal and maximal drug concentration (`0nM`
and `10000nM`, respectively) after `24h`.

How would you interpret the MDS?

``` r
meta_lung = metadata %>% 
  dplyr::filter(concentration_nM %in% c(0, 10000) & 
                time_h == 24 & 
                tissue == 'Lung') %>%
  dplyr::select(-c(time_h, tissue)) %>%
  dplyr::mutate(concentration_nM = replace(concentration_nM, concentration_nM== '0', 'basal'), 
                concentration_nM = replace(concentration_nM, concentration_nM == '10000', 'drug'))

gex_lung = expression %>%
  dplyr::select(meta_lung$geo_accession)

#write.csv(gex_lung, "data/NC60_expression.csv")

limma::plotMDS(gex_lung, labels = paste(meta_lung$sample, meta_lung$concentration_nM, sep = '_'))
```

There are different methods that we can use to conduct a differential
expression analysis, such as linear models, t-test or anovas. In our
case, we will use the package `limma`.

LIMMA is a library for the analysis of gene expression microarray data,
especially the use of linear models for analysing designed experiments
and the assessment of differential expression. The linear model and
differential expression functions apply to all gene expression
technologies, including microarrays, RNA-seq and quantitative PCR.

``` r
f = factor(meta_lung$concentration_nM, levels = c("basal","drug"))
design <- model.matrix(~0+f)
colnames(design) <- c("basal","drug")
fit <- limma::lmFit(gex_lung, design) #Fit linear model (intercept 0, so each coefficient is the mean)
cont.matrix <- limma::makeContrasts(treatment_dif = basal - drug,
levels=design)
fit2 <- limma::contrasts.fit(fit, cont.matrix)
fit2 <- limma::eBayes(fit2)
DEA_results = as.data.frame(limma::topTable(fit2,adjust.method = "BH",number = Inf)) %>% 
  tibble::rownames_to_column(var = "Gene") %>% 
  dplyr::arrange(desc(abs(t))) %>% as_tibble()
limma_genes = filter(DEA_results, adj.P.Val < 0.05)[[1]]

DEA_results %>%
  dplyr::rename(ID = Gene) %>%
  write.csv("../data/DEA_lung_NC60.csv", row.names = F)

DEA_results %>% filter(Gene%in%limma_genes)
```

By selecting the differential expressed genes, we can come back to our
collaborators. How would you explain these results to them?

We may conclude our analysis here, or we can go a step further and
conduct a functional analysis on the data…

## Session Info Details

    ## R version 4.1.1 (2021-08-10)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 22621)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United Kingdom.1252 
    ## [2] LC_CTYPE=English_United Kingdom.1252   
    ## [3] LC_MONETARY=English_United Kingdom.1252
    ## [4] LC_NUMERIC=C                           
    ## [5] LC_TIME=English_United Kingdom.1252    
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] ggplot2_3.3.5 dplyr_1.0.8  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] rstudioapi_0.13  knitr_1.39       magrittr_2.0.3   munsell_0.5.0   
    ##  [5] tidyselect_1.1.2 colorspace_2.0-3 R6_2.5.1         rlang_1.0.2     
    ##  [9] fastmap_1.1.0    fansi_1.0.3      stringr_1.4.0    tools_4.1.1     
    ## [13] grid_4.1.1       gtable_0.3.0     xfun_0.30        utf8_1.2.2      
    ## [17] DBI_1.1.3        cli_3.2.0        withr_2.5.0      htmltools_0.5.2 
    ## [21] ellipsis_0.3.2   assertthat_0.2.1 yaml_2.3.5       digest_0.6.29   
    ## [25] tibble_3.1.6     lifecycle_1.0.1  purrr_0.3.4      vctrs_0.4.0     
    ## [29] glue_1.6.2       evaluate_0.15    rmarkdown_2.14   stringi_1.7.6   
    ## [33] compiler_4.1.1   pillar_1.8.0     scales_1.2.0     generics_0.1.3  
    ## [37] pkgconfig_2.0.3
