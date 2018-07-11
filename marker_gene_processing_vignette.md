
# Filtering marker gene datasets

Here, we'll walk through some filtering strategies applicable to compositional marker gene studies. We will employ several custom R scripts designed to manipulate phyloseq-class datasets that aid in the filtering and reproducability of marker gene studies. We apply this workflow to a subset of a larger 16S community dataset, though these scripts are scalable to any marker gene and large datasets with thousands of taxa and samples. 

<div class="alert alert-block alert-info">
<b>A preface on filtering:</b> Filtering is an _essential_ step in the analyis of marker gene datasets. This workflow uses a two-pronged approach to filtering contaminants from read datasets. That is, we attempt to remove two types of contamination: <br>
<br>
1. Real biological sequences that have cross contaminated unassociated samples. <br>
<br>
2. Exogenous contaminants introduced during sample preparation. <br>
<br>
To achieve this, we'll first filter out taxa below an abundance threshold individually within each sample. We estimate this threshold empirically by measuring the number of amplicon sequence variants (ASVs) within known mock communities under various thresholds and select the threshold that returns the desired results. Next, we'll remove control samples from the dataset and then apply a dataset-wide prevalence filter to remove exogenous contaminants present as low abundance ASVs. 
</div>

***

## Before we begin...
This workflow assumes that your sequencing project has been:
-  processed with __[dada2](https://benjjneb.github.io/dada2/index.html)__ (namely, that your ASV names are the ASV nucleotide sequences)
-  converted into a __[phyloseq](https://joey711.github.io/phyloseq/)__ object


```R
library(RCurl); packageVersion("RCurl")
```

    Loading required package: bitops



    [1] ‘1.95.4.10’



```R
library(phyloseq); packageVersion("phyloseq")
```


    [1] ‘1.19.1’



```R
library(ggplot2); packageVersion("ggplot2")
```


    [1] ‘2.2.1’


## Getting started
First, we'll need to source the functions into our R environment. We can do that directly from Github with some handy functions provided by HW. 


```R
library(devtools); packageVersion("devtools")
mylink <- "https://raw.githubusercontent.com/itsmisterbrown/marker_gene_processing_scrips/master/marker_gene_processing_functions_070918.R"
source_url(mylink) #compliments of Hadley
```


    [1] ‘1.13.6’


    SHA-1 hash of file is 1e2dc1f20fa49ebb0635850acf0431ed5cac33bf


For this vignette, we'll work with a subset of a recent study of stool samples from Nigerian infants. Let's take a look at the raw dataset


```R
NI
```


    Error in eval(expr, envir, enclos): object 'NI' not found
    Traceback:




```R
(ws.threshold <- estimate.threshold(ps = NI, WSmin = 1e-4, WSmax = 2e-4, WSstep = 1e-5, controlID = "pc1"))
```

And let's take a peek at a quick visual summary


```R
ggplot(data = ws.threshold, aes(x = threshold.value, y = control.taxa.count)) + geom_line(size=2, color="steelblue2") + 
  theme(axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"))

```


```R
ggplot(data = ws.threshold, aes(x = threshold.value, y = read.percent)) + geom_line(size=2, color="orangered1") + 
  theme(axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"))

```

***
It looks like a within-sample threshold of 2e-4 is sufficient for recovering the desired taxa count for our mock communities, without sacrificing too many reads. Now that we've determined our threshold, let's apply this cutoff to our dataset with the __filter.dataset()__ function


```R
NI.f <- filter.dataset(ps = NI, controlID = "pc1", controlCAT = "timepoint", 
                      controlFACTOR = "control", WSF = 1e-3, PF = 1e-5, return.all = TRUE)
```


```R
NI.f$control.taxa.sequences
```


```R
NI.f$taxonomy.of.control.taxa
```


```R
NI.f$read.count.table[142:149, ]
```


```R
NI.f$prevalence.filter.threshold
```


```R
(NIF <- NI.f$filtered.phyloseq)
```

Next, we'll export all of the data associated with our final, filtered dataset using the __write.dataset()__ function. Unfortunately, I've never been able to get the biomformat package (v1.2.0) to play nice with my version of phyloseq, but y'all may have more luck with newer versions. In case the kinks have been worked out, I'd recommend using the __write.dataset.biom()__ function to export your dataset directly to biom format. Until then, we'll use __write.dataset()__ to perform some final minor modifications and export our dataset.

<div class="alert alert-block alert-info">
<b>Output:</b> The output of these scripts are as follows: <br>
<br>
<b>Both:</b>
-  FASTA of all ASV sequences. dada2 does not assign names to each ASV, but rather uses the nucleotide sequence as the identifier; we can leverage this feature to write an output in FASTA format of each ASV before we assign names.  <br>
<br>
<b>write.dataset.biom()</b>
-   A biom file with observation and sample metadata integrated. I've run into problems attempting to import the resulting biom back into phyloseq, but this might warrant a closer look with a newer version. <br>
<br>
<b>write.dataset()</b>
-  A standard count table (eg. OTU table) of read counts of taxa (ASVs) within each sample. <br>
-  A taxonomy table with the sequence identifier (eg. ASV1) in the first column and semicolon-separated taxonomy in the second column.<br>
-  A sample data table of all associated metadata. 
<br>

</div>

Before writing, all files written by __write.dataset()__ are formatted for easy merging into a biom file using the __[standalone biom package](http://biom-format.org/documentation/biom_conversion.html)__ without any further manipulation. Additionally, both functions rename all taxa as ASV1, ASV2, ..., ASVN, for ease of visualization and handling in downstream analysis. 


```R
fp <- paste(getwd(), "/", sep = "")
NIF.f <- write.dataset(ps = NIF, filepath = fp, fileprefix = "Nif.f")
```


```R
phyloseq::otu_table(NIF.f)[1:10, 1:10]
```
