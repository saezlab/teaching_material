Functional analysis using footprint approaches
================
Rosa Hernansaiz Ballesteros
09/12/2020

## Installing packages

The required packages to run the different steps of this analysis are in
the CRAN repository, Bioconductor or gitHub.

*NOTE: For running CARNIVAL analysis, an interactive version of **IBM
Cplex solver** is required as the network optimiser. The IBM ILOG Cplex
is freely available through Academic Initiative
[here](https://www.ibm.com/products/ilog-cplex-optimization-studio?S_PKG=CoG&cm_mmc=Search_Google-_-Data+Science_Data+Science-_-WW_IDA-_-+IBM++CPLEX_Broad_CoG&cm_mmca1=000000RE&cm_mmca2=10000668&cm_mmca7=9041989&cm_mmca8=kwd-412296208719&cm_mmca9=_k_Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB_k_&cm_mmca10=267798126431&cm_mmca11=b&mkwid=_k_Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB_k_%7C470%7C135655&cvosrc=ppc.google.%2Bibm%20%2Bcplex&cvo_campaign=000000RE&cvo_crid=267798126431&Matchtype=b&gclid=Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB).*

``` r
library(tidyverse)

source("assignPROGENyScores.r")
source("generateTFList.r")
```

## Data description and loading

The dataset we are using for this sesion is a subset of the GEO entry
`GSE116436`. This is a GEO series of the whole data set. It comprehends
transcriptome profiling upon **60 cell lines using 15 anti-cancer
agents** with different time points (2, 6, 24 hrs) and drug
concentrations (incl. null concentration, i.e. control). A short
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

  - expression data matrix: data.frame of 13226 genes (gene symbol) vs
    485 conditions

  - metada matrix: data.frame containing extra information about the
    samples

<!-- end list -->

``` r
metadata = readRDS("data/NC60_metadata.rds") %>%
  dplyr::mutate_if(is.factor, as.character) %>% 
  dplyr::rename(tissue = "tissue:ch1") %>%
  tibble::tibble() 

expression = readRDS("data/NC60_expression.rds")
```

## Getting ready for downstream analysis

First of all, we need to twick the `metadata` matrix. The information in
the `title` column is provided in an inconvinient format. So we are
going to get it in a more easy-access
way.

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

DEA_results %>% filter(Gene%in%limma_genes)
```

By selecting the differential expressed genes, we can come back to our
collaborators. How would you explain these results to them?

We may conclude our analysis here, or we can go a step further and
conduct a functional analysis on the data…

## Footprint analysis

With the diferential gene expression, or even using directly the gene
expression, we can conduct our functional analysis: Can we found deeper
differences between lapatinib treated/untreated cell lines derived from
lung cancer?

There are many ways of understand or giving answers for our biological
questions, here we are going to explore several methods of functional
analysis based on footprint approaches:

  - Activity of Transcription Factors using DoRothEA
  - Activity of cancer-related pathways using PROGENy
  - Network’s reconstruction using CARNIVAL

### Transcription factor analysis using DoRothEA

DoRothEA (Discriminant Regulon Expression Analysis) is a framework to
estimate single sample TF activities from gene expression data and
consensus TF-target DNA binding networks. The approach assumes that the
activity of a TF can be estimated from the mRNA levels of its direct
target genes (regulons). TF regulons are *signed* (to account for
activation/repression), when possible, and accompanied by a confidence
score.

> [Garcia-Alonso et al. Benchmark and integration of resources for the
> estimation of human transcription factor activities. Genome Research
> 2019 July
> 24;29(8):1363-1375.](https://genome.cshlp.org/content/early/2019/07/24/gr.240663.118.abstract)

> [Garcia-Alonso et al. Transcription Factor Activities Enhance Markers
> of Drug Sensitivity in Cancer. Cancer Research 2018 Feb
> 1;78(3):769-780.](https://www.ncbi.nlm.nih.gov/pubmed/30355619)

First, we need to **create the regulons**. We are going to use the
OmniPath’s R package to retreve the information about the TF’s targets.
Omnipath is a comprehensive collection of literature curated human and
rodent signaling pathways.

> [Turei et al. OmniPath: guidelines and gateway for literature-curated
> signaling pathway resources. Nature Methods 2016;
> 13(12):966.](https://www.nature.com/articles/nmeth.4077.epdf?author_access_token=H15Yqi8F5sADwj5dpJ0kCdRgN0jAjWel9jnR3ZoTv0Pg4abXTKgFXxt5OxdZm-zVfjCooSFBEc3ZjUdSSJ-VJQBgCAAtKMadqrq6Ck2eJVmj7x6PjR6FCEj269Fi7h6E)

``` r
#Dorothea
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B", "C"))

tf_activities_stat <- DEA_results %>%
    dplyr::select(Gene,t) %>%
    tibble::column_to_rownames(var = "Gene") %>%
    dorothea::run_viper(., 
                        regulons,
                        options = list(minsize = 5, eset.filter = FALSE, 
                        cores = 1, verbose = FALSE, nes = TRUE)) %>%
  as.data.frame() %>%
  dplyr::rename(activity = t)

# Check TF activities
ggplot(as.data.frame(tf_activities_stat), aes(x=activity)) + geom_density()
```

We have now the TF’s activity based on the expression of their targets…

  - How do you interpret the distribution of the TF’s activities?

  - Where would you put a “threshold” to select the most up/down
    regulated TF?

<!-- end list -->

``` r
tf_activities_stat %>%
  dplyr::filter( abs(activity) > 3 ) %>%
  tibble::rownames_to_column(var = "TF") %>%
  ggplot(., aes(x = activity, y = reorder(TF, activity))) + 
    geom_bar(aes(fill = activity), stat="identity") +
    theme_minimal() +
    scale_fill_gradient2(low = "#99004C", high = "#0859A2",
                         mid = "whitesmoke", midpoint = 0) +
    ylab("TFs") +
    xlab("NES")
```

How do you interpret these results? What would you say to the
collaborators?

### Pathway analysis using PROGENy

PROGENy is a method that uses a large compendium of publicly available
perturbation experiments to yield a common core of Pathway RespOnsive
GENes.

Unlike existing methods, PROGENy can:

  - recover the effect of known driver mutations

  - provide or improve strong markers for drug indications

  - distinguish between oncogenic and tumor suppressor pathways for
    patient survival

Collectively, these results show that PROGENy more accurately infers
pathway activity from gene expression than other methods.

> [Schubert et al. Perturbation-response genes reveal signaling
> footprints in cancer gene expression. Nature Communications 2018;
> 9(1):20.](https://doi.org/10.1038/s41467-017-02391-6)

We can use PROGENy as we did with DoRothEA: using the differential gene
expression to see differences between pathways’ activities. However, we
can also apply these methods directly to gene expression matrices\!

``` r
data(model_human_full, package = "progeny")

pscores <- gex_lung %>%
  as.matrix() %>%
  progeny::progeny(., z_scores = FALSE, 
                           organism = "Human", top = 100, perm = 100)
```

We can now visualise the results…

``` r
paletteLength = 100
myColor <- colorRampPalette(c("#99004C", "whitesmoke", "#0859A2"))(paletteLength)

pmap_annotation = meta_lung %>%
  tibble::column_to_rownames(var = "geo_accession") %>%
  dplyr::select(sample, concentration_nM) %>%
  dplyr::rename(drug_status = concentration_nM)

pheatmap::pheatmap(pscores, 
                   annotation_row = pmap_annotation,
                   fontsize = 14,
                   fontsize_row = 10, fontsize_col = 10,
                   color = myColor,
                   angle_col = 45, treeheight_col = 0, 
                   border_color = NA)
```

How do you interpret the pathway activities? What would you say to the
collaborators?

### Network reconstruction analysis using CARNIVAL

CARNIVAL is a framework to perform causal reasoning to infer a subset of
signalling network from transcriptomics data. It uses the Transcription
factors’ (TFs) activities from DoRothEA and the pathways scores from
PROGENy to drive the ILP algorithm.

> [Liu A\*, Trairatphisan P\*, Gjerga E\*, et al. From expression
> footprints to causal pathways: contextualizing large signaling
> networks with CARNIVAL. *bioRxiv*
> 2019](https://doi.org/10.1101/541888) (\*equal contributions).

The most *important* requirment for CARNIVAL is the network. The network
is the one that will constrain the results, and also our interpretation.

``` r
# signed and directed
omniR <-OmnipathR::import_Omnipath_Interactions() %>%
  dplyr::filter(consensus_direction == 1 & (consensus_stimulation == 1 | consensus_inhibition == 1 ))

# changing 0/1 criteria to -1/1
omniR$consensus_stimulation[which( omniR$consensus_stimulation == 0)] = -1
omniR$consensus_inhibition[which( omniR$consensus_inhibition == 1)] = -1
omniR$consensus_inhibition[which( omniR$consensus_inhibition == 0)] = 1

# check consistency on consensus sign and select only those in a SIF format
sif <- omniR[,c('source_genesymbol', 'consensus_stimulation', 'consensus_inhibition', 'target_genesymbol')] %>%
      dplyr::filter(consensus_stimulation==consensus_inhibition) %>%
      unique.data.frame() %>%
      dplyr::select(source_genesymbol, consensus_stimulation, target_genesymbol) %>%
      dplyr::rename(source = source_genesymbol , interaction = consensus_stimulation, target = target_genesymbol) %>%
      dplyr::mutate(source = gsub(":", "_", source), target = gsub(":", "_", target))

#save SIF
write_tsv(sif, "data/omnipath_carnival.tsv")
```

Once we have the network, we can use DoRothEA and PROGENy estimations to
run CARNIVAL. DoRothEA gives status of the TFs that will be reached,
while PROGENy helps the IPL algorithm to navegate. As we want to analyse
the differential expression network, we calculate PROGENy pathway scores
using the t-statistic.

*NOTE: Although DoRothEA is a must, PROGENy is optional but specially
useful when we want to reconstruct a whole network without perturbation
node to reach.*

``` r
# dorothea for CARNIVAL
tfList = generateTFList(tf_activities_stat, top=50, access_idx = 1)

# calculate progeny from DE
pathway_scores <- DEA_results %>%
  dplyr::select(Gene, t) %>%
  tibble::column_to_rownames(var = "Gene") %>%
  as.matrix() %>%
  progeny::progeny(., z_scores = FALSE, 
                           organism = "Human", top = 100, perm = 100)

# format for CARNIVAL
load(file = system.file("progenyMembers.RData",package="CARNIVAL"))

progenylist = assignPROGENyScores(progeny = pathway_scores, 
                                            progenyMembers = progenyMembers, 
                                            id = "gene", 
                                            access_idx = 1)
names(progenylist) = "t"
```

We can now reconstruct the network with CARNIVAL. CARNIVAL provides two
options: - Reconstruct the network up to the inital nodes - Select
target(s) to reach.

As we know the targets of lapatinib, we will reconstruct the network
looking for this targets

The default mode for CARNIVAL is the second opction, but we can easily
toggel this with the option `inverseCR=T`.

We can analyse the network result in different ways; in order to explore
the different options, check notebook 6 of our
[transcriptutorial](https://github.com/saezlab/transcriptutorial). Here
we will manly analyse visually the network obtainedfrom CARNIVAL.

``` r
g = carnival_result$weightedSIF %>%
    dplyr::select(Node1, Node2, Sign) %>%
    dplyr::rename(from = Node1, to = Node2, weight = Sign) %>%
    as.matrix() %>%
    igraph::graph.data.frame(., directed = T)

igraph::E(g)$color <- as.factor(carnival_result$weightedSIF$Sign)

plot(g, 
     vertex.size=1, 
     edge.arrow.size=0.1,
     layout= igraph::layout.davidson.harel,
     vertex.label.font=0.5,
     vertex.label.color='black')
legend('topleft')
```

## Session Info Details

    ## R version 4.0.2 (2020-06-22)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS  10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] forcats_0.5.0   stringr_1.4.0   dplyr_1.0.2     purrr_0.3.4    
    ## [5] readr_1.4.0     tidyr_1.1.2     tibble_3.0.4    ggplot2_3.3.2  
    ## [9] tidyverse_1.3.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.5        cellranger_1.1.0  pillar_1.4.7      compiler_4.0.2   
    ##  [5] dbplyr_2.0.0      tools_4.0.2       digest_0.6.27     lubridate_1.7.9.2
    ##  [9] jsonlite_1.7.1    evaluate_0.14     lifecycle_0.2.0   gtable_0.3.0     
    ## [13] pkgconfig_2.0.3   rlang_0.4.9       reprex_0.3.0      cli_2.2.0        
    ## [17] DBI_1.1.0         rstudioapi_0.13   yaml_2.2.1        haven_2.3.1      
    ## [21] xfun_0.19         withr_2.3.0       xml2_1.3.2        httr_1.4.2       
    ## [25] knitr_1.30        fs_1.5.0          hms_0.5.3         generics_0.1.0   
    ## [29] vctrs_0.3.5       grid_4.0.2        tidyselect_1.1.0  glue_1.4.2       
    ## [33] R6_2.5.0          fansi_0.4.1       readxl_1.3.1      rmarkdown_2.5    
    ## [37] modelr_0.1.8      magrittr_2.0.1    backports_1.2.0   scales_1.1.1     
    ## [41] ellipsis_0.3.1    htmltools_0.5.0   rvest_0.3.6       assertthat_0.2.1 
    ## [45] colorspace_2.0-0  stringi_1.5.3     munsell_0.5.0     broom_0.7.2      
    ## [49] crayon_1.3.4
