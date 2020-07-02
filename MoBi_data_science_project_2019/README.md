# Project 2 : Cellular response to drug perturbations

*Project overview and guidelines*

* [Introduction](#Introduction)
* [Objectives](#Objectives)
    - [Broad analysis](#Broad-analysis)
    - [Specific analysis](#Specific-analysis)
* [Description of data sets](#Description-of-data-sets)
* [Literature review](#Literature-review)
* [How to structure your project](#How-to-structure-your-project)
    - [Project proposal](#Project-proposal)
    - [Project](#Project)
* [Glosary](#Glosary)

Supervisors:
* Nicolàs Palacio-Escat (nicolas.palacio@bioquant.uni-heidelberg.de)
* Javier Perales-Patón (javier.perales@bioquant.uni-heidelberg.de)

Tutor:
* Julia Rühle ([julia.ruehle@web.de](mailto:julia.ruehle@web.de))

## Introduction
Cancer is a generic term for a heterogeneous group of diseases that arise in
different parts of the body. Cancer stems from genetic alterations that
transform normal cells into malignant cells. These malignant cells are
characterized by an uncontrolled growth and resistance to apoptosis (programmed
cell death). This malignant behavior is mainly achieved through genetic
alterations that affect the activity and wiring of many molecular mechanisms
like signaling and metabolic pathways. Therefore cancer cells hijack many
cellular processes enabling them to proliferate uncontrollably, avoid apoptosis
and/or immune response, migrate to other tissues (metastasis) or even obtain
drug resistance. Moreover, large-scale studies from international cancer
consortia have described that patients' tumors acquire this malignancy through
diverse individual molecular alterations. In other words, individual cancers
are different each other.

Cancer therapies aim to develop drugs that exploit these diverse aberrant
characteristics in a targeted manner, potentially resulting in an effective
inhibition of the cancer cells' growth with lesser toxicity for the patients.
These include **chemotherapy agents** that generally deplete dividing cells
through the disruption of DNA replication or chromosomal segregation. On the
other hand, **targeted therapies** hit particular key molecules (usually
proteins) involved in deregulation of a particular biological pathway that
provides benefit for the cancer cells.

This data set is indented to show how cancer cells respond to well established
anti-cancer agents in terms of **transcriptome modulation** and **cell growth**
after their exposure to the drug.

## Objectives
The main objective is that the students look into the different cellular
responses to drug perturbations in the comprehensive data set of the
NCI-60<sup>[1](#Original-paper-of-the-data-set)</sup>. For this,
transcriptomics and drug sensitivity data are provided already formatted and
subsetted. The student should explore the heterogeneous drug response to 15
anti-cancer agents, followed by a more comprehensive understanding of a more
specific case-study.

<img src="objectives.svg" alt="Summary of objectives" width="750"/>

To address the main objective, the following individual objectives are proposed
(note that the following steps are a suggestion, you can deviate from them or
use other approaches depending on your interests or findings):

### Broad analysis
* Explore gene expression profiles in treated and untreated profiles (e.g.
density plots, violin or box plots, profile/correlation heatmap...).
> Do they look already normalized?<br>
> Could you identify number of batches from the data if any?

* Explore the space of gene expression using **dimensionality reduction**
methods such as PCA. Relate the **main latent factors** to the sample metadata
(e.g. cell lines, tissue, drug mechanism of action, dose, time point, etc).
> Is there any co-variate that correlates with the main aspects of
transcriptional heterogeneity in the data set?

* Extract **individual gene expression signatures** by ranking the most
differentially expressed genes using the log2 fold-change (FC) between treated
and untreated conditions for each cell line at each time point and dose.

    ```
    log2(FC) = log2(gene expr in treated / gene expr in untreated)
             = log2(gene expr in treated) - log2(gene expr in untreated)
    ```

    > **Note:** These *individual gene expression signatures* are individual
    for each cell line. Since there are no technical replicates, these must be
    taken with caution.

* Repeat the exploratory analysis on the matrix of log2 FC using
**dimensionality reduction** methods. Relate with co-variates.
> Do they still preserved the original patterns observed at first glance?

### Specific analysis
From now on, focus your analysis in one particular anti-cancer agent from the
data set: **Choose one particular targeted therapy or a chemotherapy agent**.

#### Enrich your metadata annotation
* **Review literature** to understand what are the main factors driving the
drug sensitivity (e.g. inhibition of cancer cell growth) of the chosen drug
like for instance: mutational status or amplification of certain proteins, gene
over-expression, microsatellite instability, cell division rate, impaired DNA
repair, etc
* **Collect functional annotation of the drug response biomarkers**. Use the
supplementary data of **basal molecular profiles** of the cell lines (from CCLE
data portal) to extend the cell line annotation with those biomarkers of drug
response (mutations, gene copy-number, over-expression, ...).

#### Dissect transcriptional changes of drug perturbation in cancer cell lines
* **Perform a genome-wide t-test between treated and untreated cell lines**.
For this, you have to perform multiple individual t-tests along the
transcriptome (e.g. one t-test per gene), in order to assess the differential
gene expression between treated and untreated cell lines.
* Use the drug response **biomarkers to dissect patterns of transcriptome
modulation (clusters) using PCA**. You can stratify cell lines based on the
mutational status of the drug's biomarkers.

#### Model drug sensitivity based on the collected biomarkers
* **Perform an exploratory analysis** similarly to the one done for the broad
analysis of the perturbation data.
* Model the drug response based on the insight obtained in previous steps (e.g.
which biomarkers or covariates explain better the drug response). Use linear
models (remember that the data can be log-transformed).
* **[Optional]** Explore other advanced methods for gene expression analysis
such as the inference pathway activities. Relate to [Further reading: omics
data analysis](#Advanced-methods). Test if these pathway activities better
correlate with drug sensitivity rather than the specific biomarkers found in
the literature.

## Description of data sets
Data matrices stored in individual RDS (R data format) and TSV files are
provided at:

* **Drug perturbation data from cell lines** (transcriptome modulation and
inhibition of cell growth)
  - Gene Expression profiles: https://figshare.com/s/db1411debcfbe2618d2f
  - Drug sensitivity assay: https://figshare.com/s/074e0120fe5e682f3d14
* **Basal molecular profiles** of cancer cell lines (mutations + gene
copy number alterations + basal gene expression):
https://figshare.com/s/fc0c71246dc192982a3c
* **Feature annotation** (cell lines and drugs):
https://figshare.com/s/efb6a529eaf2d4dc6432

Please **download** every single file described in the three repositories
above. The files are summarized bellow:

| Type | Filename | R object | Description |
| --- | ---- | --- | --- |
| Drug perturbation | `NCI_TPW_gep_treated.rds`  | Matrix  | Gene expression profiling of transcriptome modulation after drug exposure - **treated cell lines**. Rows are genes, and columns are the samples. Gene expression values are log2-transformed. <br> Source: NCI-TPW |
| Drug perturbation | `NCI_TPW_gep_untreated.rds`  | Matrix  | Gene expression profiling of transcriptome modulation after drug exposure - **untreated cell lines** (null concentration, controls). Rows are genes, and columns are the samples. Gene expression values are log2-transformed. <br> Source: NCI-TPW |
| Drug perturbation | `NCI_TPW_metadata.tsv`  | data.frame | Experiment metadata for the gene expression profiling above. <br> Source: NCI-TPW |
| Drug perturbation | `NegLogGI50.rds` | Matrix | Drug sensitivity profiling (50% growth inhibition). Rows are drugs, and columns are the cancer cell lines. Values are -log10 transformed of drug concentration required for 50% growth inhibition (EC50), thus higher values indicate more sensitivity. <br> Source: NCI-60 |
| Basal molecular profile  | `CCLE_mutations.rds`  | data.frame  | Somatic mutations (Single-Nucleotide-Variants and INDELS) from cancer cell lines. The data is formatted as [MAF](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/#protected-maf-file-structure) <br> Source: CCLE |
| Basal molecular profile | `CCLE_copynumber.rds` | data.frame | Gene Copy-number alterations. Values are log2-transformed gene copy-number ratio (i.e. log2(gene CN / 2); it is divided by 2 because in diploid genomes usually there are 2 copies of a gene. <br> Source: CCLE |
| Basal molecular profiles  | `CCLE_basalexpression.rds` | data.frame  | Basal gene expression of cell lines. This might be useful to discern which cell lines over-expressed a gene as compared to the others. Values are log2-transformed. <br> Source: CCLE |
| Feature annotation | `cellline_annotation.tsv` | data.frame | Initial cell line metadata. Please, expand with the omics profiles provided.  |
| Feature annotation | `drug_annotation.tsv`  | data.frame | Drug annotation: mechanism of action, targets, etc. |

You can load this data using the following function in `R`:

```
mat <- readRDS("filename.rds")
```

The **data types** are described as follows:
* **Gene Expression** : These values reflect the level of gene expression,
higher values suggests over expression of genes and vice-versa.
* **Drug sensitivity** : It reflects the drug concentration required to inhibit
the growth of a cancer cell by 50% (usually termed effective concentration to
achieve half of the maximum response, shortened EC50). The values are in -log10
scale of the concentration required for the 50% inhibition. Therefore, higher
values indicate that the cell line is more sensitive to the drug than the
others with lower values because it requires lower drug concentration to reach
the 50% of growth inhibition.
* **Gene Copy-Number Alterations** : Values are log2(CN ratio), thus log2(CN ratio)=0 indicate a neutral number of copies (2 copies in a diploid genome). Higher values indicates amplifications (e.g. log2(CN ratio) > 1) and lower values indicates deletions (e.g. log2(CN ratio) < -1).
* **Somatic point mutations** : DNA mutation consequences of single-nucleotide variants and indels affecting protein coding genes.

## Literature Review
#### Original paper of the data set
> Monks A *et al*. The NCI Transcriptional Pharmacodynamics Workbench: a tool
to examine dynamic expression profiling of therapeutic response in the NCI-60
cell line panel. **Cancer Res**. 2018 Dec 15;78(24):6807-6817. doi:
[10.1158/0008-5472.CAN-18-0989](https://doi.org/10.1158/0008-5472.CAN-18-0989).<br>
[Link to pubmed](https://www.ncbi.nlm.nih.gov/pubmed/30355619).

#### Conceptual approach
> Lamb J *et al*. The Connectivity Map: using gene-expression signatures to
connect small molecules, genes, and disease. **Science**. 2006 Sep
29;313(5795):1929-35. doi:
[10.1126/science.1132939](https://doi.org/10.1126/science.1132939).

#### Further reading: getting biological insights derived from this kind of data
* Chen B *et al*. Reversal of cancer gene expression correlates with drug
efficacy and reveals therapeutic targets. **Nat Commun**. 2017 Jul 12;8:16022.
doi: [10.1038/ncomms16022](https://doi.org/10.1038/ncomms16022).
* Szalai B *et al*. Signatures of cell death and proliferation in perturbation transcriptomics data - from confounding factor to effective prediction. **bioRxiv** 2018. doi: [10.1101/454348](https://doi.org/10.1101/454348)
* Iorio F *et al*. A Landscape of Pharmacogenomic Interactions in Cancer. **Cell**. 2016 Jul 28;166(3):740-754. doi: [10.1016/j.cell.2016.06.017](https://doi.org/10.1016/j.cell.2016.06.017)
* Iorio F *et al*. Discovery of drug mode of action and drug repositioning from transcriptional responses. **Proc Natl Acad Sci U S A**. 2010 Aug 17;107(33):14621-6. doi: [10.1073/pnas.1000138107](https://doi.org/10.1073/pnas.1000138107)


#### Advanced methods
These methods unlock the possibility to infer the transcription factor
activities and pathway activities.
* [PROGENy](https://saezlab.github.io/progeny/)
* [DoRothEA](https://saezlab.github.io/DoRothEA/)
* [GSVA](https://bioconductor.org/packages/release/bioc/html/GSVA.html)

## How to structure your project
### Project proposal
You first task will be to define a project proposal, which should include:
* List of planned analysis
* Choose a drug of interest for the specific analysis (**MUST** be different
for each group).
* Milestones (important achievements)
* Deliverables (what kind of result(s) will we produce for each milestone?)
* Approximate timetable

You will present this project proposal together with a literature review on the
subject 3 weeks after the beginning of the semester (10 minutes presentation +
5 minutes discussion).

### Project
Your project **MUST** contain the following elements:
* Descriptive statistics about the data sets, including graphical
representations
* A dimension reduction analysis (PCA, clustering or k-means, ...)
* Statistical tests (t-test, proportion tests, etc)
* A linear regression analysis, either uni- or multivariate

#### General plan
1. **Broad analysis** - All groups must perform the following (specific
approaches may differ across groups):
    * Familiarize with the data (descriptive statistics, basic plots, ...).
    * Assess the need for normalization, batch correction and/or presence of
    outliers.
    * Apply the changes if necessary and observe the results on the data
    distribution.
    * Data reduction and pattern/cluster identification (may require
    literature and/or database research).
2. **Specific analysis** - Each group must choose a drug on which the
following analysis will be focused:
    * Metadata collection (literature, databases...).
    * Statistical analysis of transcriptional changes upon drug treatment.
    * Modeling drug-response and performance assessment.

---

## Glosary
* **Cancer cell-line**: is an individual cell lineage of a cancer patient which
has established during years. Cancer cell lines are fully profiled in terms of
genomic aberrations, doubling time (cell division time), etc.
* **Gene Expression Profile** (GEP): represents the gene-level measurements of
a sample in terms of transcription levels (quantified transcripts).
* **Drug sensitivity**: refers to the magnitude of sensitiveness in terms of
cellular growth e.g. in the context of cancer cell-lines, usually is regarded
as the level of growth inhibition.
* **Targeted therapy**: refers to a treatment (usually a drug) whose mechanism
of action is to target a specific molecule in a particular pathway, cellular
process and/or cell type.
* **Chemotherapy agent**: in contrast to targeted therapy, chemotherapy affects
broader range of cellular processes or cells such as DNA replication, or cell
stress. Since are citotoxic agents with a wider spectrum of action than
targeted therapies, these are more toxic in terms of whole tissue/individual.
* **PCA** : Principal Component Analysis.
