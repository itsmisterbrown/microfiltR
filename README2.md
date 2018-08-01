
# Processing and filtering marker gene datasets

Here, we'll walk through some filtering strategies applicable to compositional marker gene studies. We will employ several custom R functions designed to manipulate phyloseq-class datasets that aid in the processing, filtering, and reproducability of marker gene studies. We apply this workflow to a moderately sized 16S community dataset generated on an Illumina HiSeq, though these strategies should be widely useful to any marker gene and easily scalable to large datasets with thousands of taxa and samples. _For reference, this entire notebook can be run in about 30 seconds on a 13" 2013 Macbook Pro._<br>
<br>
Briefly, the aims of this vignette and associated functions are to:<br>
1. Estimate and apply dataset-specific filtering parameters.<br>
2. Provide a set of robust and highly reproducible processing, filtering, and exporting functions for marker gene studies.<br>
3. Faciltate transparency and reproducability of marker gene data analysis, regardless of user skill level.<br>

<div class="alert alert-block alert-info">
<b>A preface on filtering:</b> Filtering is an <i>essential</i> step in the analyis of marker gene datasets. This workflow uses a two-pronged approach (within sample and across sample) to filtering contaminants from read datasets. That is, we attempt to remove two types of contamination: <br>
<br>
1. Real biological sequences that have cross contaminated unassociated samples. <br>
<br>
2. Exogenous contaminants introduced during sample preparation. <br>
<br>
To achieve this, we first estimate and apply <b>within-sample (WS)</b> filters to remove taxa below an abundance threshold individually within each sample. We estimate this threshold empirically by measuring the number of amplicon sequence variants (ASVs) within known mock communities under various thresholds and select the threshold that returns the appropriate ASV count. Next, we (optionally) remove control samples from the dataset and then apply several <b>across-sample (AS)</b> filters (relative abundance, prevalence, and/or coefficient of variation (CV)) to remove exogenous contaminants present as low abundance ASVs. 
</div>

***

## Before we begin...
This workflow assumes that your sequencing project: 
-  has been converted into a __[phyloseq](https://joey711.github.io/phyloseq/)__ object.

and ideally:<br>
-  includes a mock community or other positive control of known composition. <br>
-  has been processed with __[dada2](https://benjjneb.github.io/dada2/index.html)__ (this is not necessary, but we leverage the output of dada2 to align control sequences and write filtered FASTA files).<br>
<br>
Now, let's load some required packages and begin!


```R
suppressMessages(library(RCurl)); packageVersion("RCurl")
```


    [1] ‘1.95.4.10’



```R
library(phyloseq); packageVersion("phyloseq")
```


    [1] ‘1.19.1’



```R
suppressMessages(library(cowplot)); packageVersion("cowplot")
```


    [1] ‘0.9.2’



```R
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw()) #gray is a bummer
```


    [1] ‘2.2.1’



```R
suppressMessages(library(tidyr)); packageVersion("tidyr")
```


    [1] ‘0.8.1’



```R
library(devtools); packageVersion("devtools")

```


    [1] ‘1.13.6’



```R
suppressMessages(library(ShortRead)); packageVersion("ShortRead")

```


    [1] ‘1.32.1’


## Getting started
First, we'll need to source these functions into our R environment. We can do that directly from Github with some handy tools provided by HW. 


```R
myfunctions <- "https://github.com/itsmisterbrown/marker_gene_processing_scripts/raw/master/marker_gene_processing_functions.R"
source_url(myfunctions) #compliments of Hadley
```

    SHA-1 hash of file is 418450e427b6f3a3661153735841ebd5d7c7933e


Now, let's load the dataset and subset it to remove taxa not annotated at the Phylum level and identified as Chloroplast sequences.


```R
#load the dataset
NI <- readRDS(file = "/Users/bpb/Desktop/marker_gene_filtering_functions/marker_gene_processing_vignette/NI.dada2.RDS")

NI %>%
  subset_taxa(Phylum != "NA" & 
                Phylum != "Cyanobacteria/Chloroplast") -> NI

table(tax_table(NI)[, "Phylum"], exclude = NULL)
```


    
                  Acidobacteria              Actinobacteria 
                             22                         629 
                Armatimonadetes               Bacteroidetes 
                              6                         359 
       candidate_division_WPS-1 Candidatus_Saccharibacteria 
                              1                           6 
                    Chloroflexi             Deferribacteres 
                             17                           1 
            Deinococcus-Thermus                  Firmicutes 
                             10                        1211 
                   Fusobacteria            Gemmatimonadetes 
                             25                           1 
                    Nitrospirae              Planctomycetes 
                              3                          59 
                 Proteobacteria                Spirochaetes 
                            815                           3 
                  Synergistetes                 Tenericutes 
                              1                           4 
                Verrucomicrobia                        <NA> 
                              3                           0 


For this vignette, we'll work with a subset of a recent study of stool samples from Nigerian infants. Let's take a look at the raw dataset


```R
NI
```


    phyloseq-class experiment-level object
    otu_table()   OTU Table:         [ 3176 taxa and 162 samples ]
    sample_data() Sample Data:       [ 162 samples by 9 sample variables ]
    tax_table()   Taxonomy Table:    [ 3176 taxa by 7 taxonomic ranks ]


Our unfiltered object has 3,176 taxa (that's a lot!), though some of those bugs are likely contaminants or irrelevent stochastic associates. We'll use some filtering functions to pare that number down to a smaller set of  bugs that we can interrogate with downstream analyses. <br>
<br>
## Filtering
<br>

The first step in our workflow is to estimate the threshold at which our data suggest cross contamination exists. Because cross contamination cannot be corrected with standard global filters, we need to estimate the extent of cross contamination and apply our filter individually within each sample. We estimate this filtering threshold against a known positive control library (HMP Mock; HM-782D) of 20 taxa using the __estimate.WSthreshold()__ function. However, due to operon variation along 16S, __we expect 23 unique ASVs along the region sequenced (V6)__. In addition to threshold estimation, __ws.threshold()__ also performs a local alignment of each ASV sequence within your postiive control against a fasta file to determine the number of exact matches between the actual composition of your sequenced control and its theoretical composition. <br>
<br>
__Note:__ sequence alignment is only supported for datasets processed by DADA2. DADA2 does not assign names to each ASV, but rather uses the nucleotide sequence as the identifier; we leverage this feature to perform sequence alignment and write output in FASTA format.
<br>
<br>

<div class="alert alert-block alert-info">
<b><font size="5">estimate.WSthreshold</font></b><br>
Performs WS filter threshold estiamtion by iterative application over a range of thresholds. <br>
<br>
<b>Input:</b> <br>
<b>ps:</b> a phyloseq-class dataset<br>
<br>
<b>Parameters:</b><br>
<b>WSrange:</b> The range of WS filter thresholds to iterate over, followed by the estimation interval.<br>
<b>controlID:</b> The sample ID of a positive control library, if applicable. Used for reporting ASV read count, taxonomy, and sequence alignment of taxa in the positive control.<br>
<b>controlFASTA:</b> the path to a FASTA file of nucleotide sequences to be used for alignment to positive control sequences, if applicable.<br>
<br>
<b>Output:</b> <br>
A dataframe of:<br>
<b>threshold.value:</b> The WS threshold applied in that iteration. <br>
<b>control.taxa.count:</b> The count of taxa in the given positive control. <br>
<b>control.taxa.matches:</b> The number of perfect alignments between taxa in the given control and the controlFASTA file. <br>
<b>read.count:</b> The read count in the positive control at the listed threshold value. <br>
<b>read.percent:</b> The read percent (vs. pre-filter read count) in the positive control at the listed threshold value. <br>

</div>

<br>
<div class="alert alert-block alert-success">
<b>TIP:</b> All of these functions support <b>vectorized input</b> for both convenience and ease of visualization! The format for the vectorized input used here is: <b>lower threshold:higher threshold, estimation interval</b>. For example, below, we will estimate over the range from 1e-5 to 2e-4 with an estimation interval of 1e-5.
</div>


```R
#load a FASTA file of full length 16S sequences for each operon of each community member
mockf <- "/Users/bpb/Desktop/HMP_MOCK.v35.fasta"

(ws.threshold <- estimate.WSthreshold(ps = NI, WSrange = c(1e-5:2e-4, 1e-5), controlID = "pc1", controlFASTA = mockf))
```

    Estimating filtering statistics from WS thresholds 1e-05 to 2e-04 by 1e-05



<table>
<thead><tr><th scope=col>threshold.value</th><th scope=col>control.taxa.count</th><th scope=col>control.taxa.matches</th><th scope=col>read.count</th><th scope=col>read.percent</th></tr></thead>
<tbody>
	<tr><td>0.00001 </td><td>36      </td><td>23      </td><td>92978550</td><td>99.98106</td></tr>
	<tr><td>0.00002 </td><td>35      </td><td>23      </td><td>92949748</td><td>99.95009</td></tr>
	<tr><td>0.00003 </td><td>30      </td><td>22      </td><td>92923580</td><td>99.92195</td></tr>
	<tr><td>0.00004 </td><td>30      </td><td>22      </td><td>92899476</td><td>99.89603</td></tr>
	<tr><td>0.00005 </td><td>28      </td><td>22      </td><td>92876359</td><td>99.87117</td></tr>
	<tr><td>0.00006 </td><td>27      </td><td>22      </td><td>92856037</td><td>99.84932</td></tr>
	<tr><td>0.00007 </td><td>27      </td><td>22      </td><td>92835116</td><td>99.82682</td></tr>
	<tr><td>0.00008 </td><td>26      </td><td>22      </td><td>92815880</td><td>99.80614</td></tr>
	<tr><td>0.00009 </td><td>25      </td><td>22      </td><td>92797926</td><td>99.78683</td></tr>
	<tr><td>0.00010 </td><td>25      </td><td>22      </td><td>92779812</td><td>99.76735</td></tr>
	<tr><td>0.00011 </td><td>25      </td><td>22      </td><td>92764199</td><td>99.75056</td></tr>
	<tr><td>0.00012 </td><td>25      </td><td>22      </td><td>92748930</td><td>99.73415</td></tr>
	<tr><td>0.00013 </td><td>25      </td><td>22      </td><td>92732313</td><td>99.71628</td></tr>
	<tr><td>0.00014 </td><td>25      </td><td>22      </td><td>92713406</td><td>99.69595</td></tr>
	<tr><td>0.00015 </td><td>24      </td><td>21      </td><td>92698673</td><td>99.68010</td></tr>
	<tr><td>0.00016 </td><td>24      </td><td>21      </td><td>92682842</td><td>99.66308</td></tr>
	<tr><td>0.00017 </td><td>24      </td><td>21      </td><td>92668320</td><td>99.64747</td></tr>
	<tr><td>0.00018 </td><td>24      </td><td>21      </td><td>92653169</td><td>99.63117</td></tr>
	<tr><td>0.00019 </td><td>22      </td><td>20      </td><td>92641098</td><td>99.61819</td></tr>
	<tr><td>0.00020 </td><td>22      </td><td>20      </td><td>92629339</td><td>99.60555</td></tr>
</tbody>
</table>



And let's take a peek at a quick visual summary


```R
plot1 <- ggplot(data = ws.threshold, aes(x = threshold.value, y = control.taxa.count)) + geom_line(size=2, color="steelblue2") + 
  theme(axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 10, colour = "black", angle = 315, vjust = 0.7),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"))

plot2 <- ggplot(data = ws.threshold, aes(x = threshold.value, y = control.taxa.matches)) + geom_line(size=2, color="orange") + 
  theme(axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 10, colour = "black", angle = 315, vjust = 0.7),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"))

plot3 <- ggplot(data = ws.threshold, aes(x = threshold.value, y = read.percent)) + geom_line(size=2, color="orangered1") + 
  theme(axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 10, colour = "black", angle = 315, vjust = 0.7),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black")) 

plot_grid(plot1, plot2, plot3, labels = "AUTO")


```




![png](output_19_1.png)


***
Though it's not a perfect match, we'll use a WS threshold of 1.4e-4, as it yields 22/23 mock community matches and significantly pares down contaminant ASVs without sacrificing too many reads. Now that we've determined our WS threshold, let's estimate our AS filter thresholds with this value as a fixed parameter. There are several global, AS filtering strategies for removing exogenous contaminants from marker gene datasets. Here, we focus on three popular approaches that can all be estimated and applied using the suite of functions sourced above. <br>
<br>

<div class="alert alert-block alert-info">
<b><font size="5">estimate.ASthreshold</font></b><br>
Performs AS filter threshold estiamtion by iterative application over a range of threshold values and filtering strategies. Supports the (optional) application of fixed thresholds if you want to optimize one or multiple threshold with the other(s) fixed. <br>
<br>
<b>Input:</b> <br>
<b>ps:</b> a phyloseq-class dataset.<br>
<br>
<b>Parameters (all are optional):</b><br>
<b>WST:</b> The WS threshold to apply before estimation of AS thresholds.<br>
<b>mdCAT:</b> The sample data category used in conjunction with mdFACTOR for metadata based filtering.<br>
<b>mdFACTOR:</b> The factor within mdCAT used for metadata based filtering.<br>
<b>mdNEGATIVE:</b> A logical flag indicating whether metadata-based filtering should be positive (remove matching samples) or negative (remove non matching samples).<br>
<b>minLIB:</b> The minimum sequence depth for a sample to be retained.<br>
<b>RArange:</b> The range of relative abundance filter thresholds to iterate over, followed by the estimation interval. <br>
<b>CVrange:</b> The range of CV filter thresholds to iterate over, followed by the estimation interval.<br>
<b>PFrange:</b> The range of prevalence filter thresholds to iterate over, followed by the estimation interval. <br>
<b>RAT:</b> Relative Abundance Threshold; the fixed percentage of reads across the dataset for which an ASV must be present for that ASV to be retained. <br>
<b>CVT:</b> Coefficient of Variation Threshold; the fixed lowest CV value of an ASV for that ASV to be retained.<br>
<b>PFT:</b> Prevalence Filter Threshold; the fixed percentage of samples across the dataset for which an ASV must be present for that ASV to be retained. <br>

<br>
<b>Output:</b> <br>
<b>relative.abundance.filtering.stats:</b> The ASV count remaining in the dataset at the given relative abundance threshold. <br>
<b>CV.filtering.stats:</b> The ASV count remaining in the dataset at the given CV threshold. <br>
<b>prevalence.filtering.stats:</b> The ASV count remaining in the dataset at the given prevalence threshold. <br>
<b>ASV.filtering.stats:</b> A dataframe listing the read count, read percent, prevalence count, prevalence percent, CV, taxonomic ranks, and sequence identity of each ASV.<br>

</div> 


We'll use the __estimate.ASthreshold()__ function to simultaneously estimate all three filter thresholds.


```R
as.threshold <- estimate.ASthreshold(ps = NI, WST = 1.4e-4, minLIB = 4000, mdCAT = "timepoint", 
                                      mdFACTOR = "control", Prange = c(0.01:0.1, 0.005), CVrange = c(0:5, 0.1), 
                                      RArange = c(1e-7:1e-5, 1e-7))
```

    Removing 13 samples with read count < 4000
    Applying WS filter threshold of 0.00014
    Removing 6 samples matching metadata identifiers timepoint:control
    Estimating filtering statistics from relative abundance thresholds 1e-07 to 1e-05 by 1e-07
    Estimating filtering statistics from CV thresholds 0 to 5 by 0.1
    Estimating filtering statistics from prevalence thresholds 0.01 to 0.1 by 0.005



```R
as.threshold$relative.abundance.filtering.stats[1:10,]
```


<table>
<thead><tr><th scope=col>relative.abundance.filter</th><th scope=col>ASV.count</th></tr></thead>
<tbody>
	<tr><td>1e-07</td><td>887  </td></tr>
	<tr><td>2e-07</td><td>837  </td></tr>
	<tr><td>3e-07</td><td>794  </td></tr>
	<tr><td>4e-07</td><td>770  </td></tr>
	<tr><td>5e-07</td><td>761  </td></tr>
	<tr><td>6e-07</td><td>743  </td></tr>
	<tr><td>7e-07</td><td>732  </td></tr>
	<tr><td>8e-07</td><td>723  </td></tr>
	<tr><td>9e-07</td><td>708  </td></tr>
	<tr><td>1e-06</td><td>695  </td></tr>
</tbody>
</table>




```R
as.threshold$CV.filtering.stats[1:10,]
```


<table>
<thead><tr><th scope=col>CV.filter</th><th scope=col>ASV.count</th></tr></thead>
<tbody>
	<tr><td>0.0</td><td>978</td></tr>
	<tr><td>0.1</td><td>978</td></tr>
	<tr><td>0.2</td><td>978</td></tr>
	<tr><td>0.3</td><td>978</td></tr>
	<tr><td>0.4</td><td>978</td></tr>
	<tr><td>0.5</td><td>978</td></tr>
	<tr><td>0.6</td><td>978</td></tr>
	<tr><td>0.7</td><td>978</td></tr>
	<tr><td>0.8</td><td>978</td></tr>
	<tr><td>0.9</td><td>978</td></tr>
</tbody>
</table>




```R
as.threshold$prevalence.filtering.stats[1:10,]
```


<table>
<thead><tr><th scope=col>prevalence.filter</th><th scope=col>ASV.count</th></tr></thead>
<tbody>
	<tr><td>0.010</td><td>620  </td></tr>
	<tr><td>0.015</td><td>416  </td></tr>
	<tr><td>0.020</td><td>416  </td></tr>
	<tr><td>0.025</td><td>328  </td></tr>
	<tr><td>0.030</td><td>280  </td></tr>
	<tr><td>0.035</td><td>248  </td></tr>
	<tr><td>0.040</td><td>248  </td></tr>
	<tr><td>0.045</td><td>217  </td></tr>
	<tr><td>0.050</td><td>197  </td></tr>
	<tr><td>0.055</td><td>197  </td></tr>
</tbody>
</table>




```R
as.threshold$ASV.filtering.stats[1:10,]
```


<table>
<thead><tr><th scope=col>ASV.read.count</th><th scope=col>ASV.read.percent</th><th scope=col>ASV.prevalence</th><th scope=col>ASV.prevalence.percent</th><th scope=col>ASV.CV</th><th scope=col>Kingdom</th><th scope=col>Phylum</th><th scope=col>Class</th><th scope=col>Order</th><th scope=col>Family</th><th scope=col>Genus</th><th scope=col>Species</th><th scope=col>ASV.ID</th></tr></thead>
<tbody>
	<tr><td>16756456                                                                                                                  </td><td>18.075627                                                                                                                 </td><td>144                                                                                                                       </td><td>96.64430                                                                                                                  </td><td>1.337880                                                                                                                  </td><td>Bacteria                                                                                                                  </td><td>Firmicutes                                                                                                                </td><td>Bacilli                                                                                                                   </td><td>Lactobacillales                                                                                                           </td><td>Streptococcaceae                                                                                                          </td><td>Streptococcus                                                                                                             </td><td>lactarius/salivarius                                                                                                      </td><td>GGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCCGATGCTATTTCTAGAGATAGGAAGTTACTTCGGTACATCGGTGACAGG</td></tr>
	<tr><td> 9524897                                                                                                                  </td><td>10.274755                                                                                                                 </td><td>145                                                                                                                       </td><td>97.31544                                                                                                                  </td><td>1.915157                                                                                                                  </td><td>Bacteria                                                                                                                  </td><td>Firmicutes                                                                                                                </td><td>Bacilli                                                                                                                   </td><td>Lactobacillales                                                                                                           </td><td>Enterococcaceae                                                                                                           </td><td>Enterococcus                                                                                                              </td><td>NA                                                                                                                        </td><td>GGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCTTTGACCACTCTAGAGATAGAGCTTCCCCTTCGGGGGCAAAGTGACAGG</td></tr>
	<tr><td> 5018847                                                                                                                  </td><td> 5.413961                                                                                                                 </td><td>139                                                                                                                       </td><td>93.28859                                                                                                                  </td><td>1.598313                                                                                                                  </td><td>Bacteria                                                                                                                  </td><td>Actinobacteria                                                                                                            </td><td>Actinobacteria                                                                                                            </td><td>Bifidobacteriales                                                                                                         </td><td>Bifidobacteriaceae                                                                                                        </td><td>Bifidobacterium                                                                                                           </td><td>breve/longum                                                                                                              </td><td>GGGCCCGCACAAGCGGCGGAGCATGCGGATTAATTCGATGCAACGCGAAGAACCTTACCTGGGCTTGACATGTTCCCGACGATCCCAGAGATGGGGTTTCCCTTCGGGGCGGGTTCACAGGT</td></tr>
	<tr><td> 4489657                                                                                                                  </td><td> 4.843110                                                                                                                 </td><td>121                                                                                                                       </td><td>81.20805                                                                                                                  </td><td>2.306163                                                                                                                  </td><td>Bacteria                                                                                                                  </td><td>Actinobacteria                                                                                                            </td><td>Actinobacteria                                                                                                            </td><td>Coriobacteriales                                                                                                          </td><td>Coriobacteriaceae                                                                                                         </td><td>Collinsella                                                                                                               </td><td>aerofaciens                                                                                                               </td><td>GGGCCCGCACAAGCAGCGGAGCATGTGGCTTAATTCGAAGCAACGCGAAGAACCTTACCAGGGCTTGACATATGGGTGAAGCGGGGGAGACCCCGTGGCCGAGAGGAGCCCATACAGGTGGT</td></tr>
	<tr><td> 3923735                                                                                                                  </td><td> 4.232635                                                                                                                 </td><td>132                                                                                                                       </td><td>88.59060                                                                                                                  </td><td>2.912638                                                                                                                  </td><td>Bacteria                                                                                                                  </td><td>Proteobacteria                                                                                                            </td><td>Gammaproteobacteria                                                                                                       </td><td>Enterobacteriales                                                                                                         </td><td>Enterobacteriaceae                                                                                                        </td><td>Salmonella                                                                                                                </td><td>NA                                                                                                                        </td><td>GGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACAGAACTTTCCAGAGATGGATTGGTGCCTTCGGGAACTGTGAGACAGG</td></tr>
	<tr><td> 3708829                                                                                                                  </td><td> 4.000811                                                                                                                 </td><td>131                                                                                                                       </td><td>87.91946                                                                                                                  </td><td>2.333284                                                                                                                  </td><td>Bacteria                                                                                                                  </td><td>Proteobacteria                                                                                                            </td><td>Gammaproteobacteria                                                                                                       </td><td>Enterobacteriales                                                                                                         </td><td>Enterobacteriaceae                                                                                                        </td><td>Escherichia/Shigella                                                                                                      </td><td>NA                                                                                                                        </td><td>GGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGG</td></tr>
	<tr><td> 3504728                                                                                                                  </td><td> 3.780642                                                                                                                 </td><td>121                                                                                                                       </td><td>81.20805                                                                                                                  </td><td>3.872899                                                                                                                  </td><td>Bacteria                                                                                                                  </td><td>Firmicutes                                                                                                                </td><td>Bacilli                                                                                                                   </td><td>Bacillales                                                                                                                </td><td>Staphylococcaceae                                                                                                         </td><td>Staphylococcus                                                                                                            </td><td>NA                                                                                                                        </td><td>GGACCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAAATCTTGACATCCTTTGACAACTCTAGAGATAGAGCCTTCCCCTTCGGGGGACAAAGTGACA</td></tr>
	<tr><td> 2967771                                                                                                                  </td><td> 3.201412                                                                                                                 </td><td>123                                                                                                                       </td><td>82.55034                                                                                                                  </td><td>2.977641                                                                                                                  </td><td>Bacteria                                                                                                                  </td><td>Firmicutes                                                                                                                </td><td>Bacilli                                                                                                                   </td><td>Lactobacillales                                                                                                           </td><td>Enterococcaceae                                                                                                           </td><td>Enterococcus                                                                                                              </td><td>NA                                                                                                                        </td><td>GGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCTTTGACCACTCTAGAGATAGAGCTTTCCCTTCGGGGACAAAGTGACAGG</td></tr>
	<tr><td> 2834047                                                                                                                  </td><td> 3.057161                                                                                                                 </td><td>102                                                                                                                       </td><td>68.45638                                                                                                                  </td><td>2.920868                                                                                                                  </td><td>Bacteria                                                                                                                  </td><td>Firmicutes                                                                                                                </td><td>Bacilli                                                                                                                   </td><td>Lactobacillales                                                                                                           </td><td>Streptococcaceae                                                                                                          </td><td>Streptococcus                                                                                                             </td><td>NA                                                                                                                        </td><td>GGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCCGATGCTATTCCTAGAGATAGGAAGTTTCTTCGGAACATCGGTGACAGG</td></tr>
	<tr><td> 2677002                                                                                                                  </td><td> 2.887752                                                                                                                 </td><td>112                                                                                                                       </td><td>75.16779                                                                                                                  </td><td>4.086987                                                                                                                  </td><td>Bacteria                                                                                                                  </td><td>Firmicutes                                                                                                                </td><td>Bacilli                                                                                                                   </td><td>Bacillales                                                                                                                </td><td>Staphylococcaceae                                                                                                         </td><td>Staphylococcus                                                                                                            </td><td>NA                                                                                                                        </td><td>GGACCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAAATCTTGACATCCTCTGACCCCTCTAGAGATAGAGTTTTCCCCTTCGGGGGACAGAGTGACA</td></tr>
</tbody>
</table>




```R
RA.stats <- as.threshold$relative.abundance.filtering.stats
CV.stats <- as.threshold$CV.filtering.stats
P.stats <- as.threshold$prevalence.filtering.stats
ASV.stats <- as.threshold$ASV.filtering.stats
ASV.stats.c <- ASV.stats[complete.cases(ASV.stats[,5]),] #remove filtered taxa


plot4 <- ggplot(data = RA.stats, aes(x = relative.abundance.filter, y = ASV.count)) + geom_line(size=2, color="steelblue2") + 
  theme(axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 10, colour = "black", angle = 315, vjust = 0.7),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"))

plot5 <- ggplot(data = CV.stats, aes(x = CV.filter, y = ASV.count)) + geom_line(size=2, color="orangered1") + 
  theme(axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 10, colour = "black", angle = 315, vjust = 0.7),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black")) 

plot6 <- ggplot(data = P.stats, aes(x = prevalence.filter, y = ASV.count)) + geom_line(size=2, color="forestgreen") + 
  theme(axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 10, colour = "black", angle = 315, vjust = 0.7),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black")) 


plot_grid(plot4, plot5, plot6, labels = "AUTO")

```




![png](output_28_1.png)


Now that we've taken a peek at the effect of multiple thresholds and filters on the dataset read count, let's see how these filters will affect the taxonomic profile of our ASVs. For each plot, I've included a solid black line at the filter threshold that I'll apply in the next step. __Unfiltered taxa exist in the space above and right of the horizontal and vertical lines, respectively.__ <br>


```R
ggplot(data = ASV.stats.c, aes(x = ASV.CV, y = ASV.prevalence.percent, color=Phylum)) + geom_point(size=3) + 
  facet_wrap("Phylum", scales = "fixed") + 
  theme(axis.text.x = element_text(size = 10, angle = 315, vjust = 0.2), legend.position = "none") + 
geom_hline(yintercept = 2)

```




![png](output_30_1.png)



```R
ggplot(data = ASV.stats.c, aes(x = ASV.CV, y = ASV.read.percent, color=Phylum)) + geom_point(size=3) + 
  facet_wrap("Phylum", scales = "fixed") + scale_y_log10() +
  theme(axis.text.x = element_text(size = 10, angle = 315, vjust = 0.2), legend.position = "none") + 
geom_hline(yintercept = 1e-5)

```




![png](output_31_1.png)



```R
ggplot(data = ASV.stats.c, aes(x = ASV.prevalence.percent, y = ASV.read.percent, color=Phylum)) + geom_point(size=3) + 
  facet_wrap("Phylum", scales = "fixed") + scale_y_log10() +
  theme(axis.text.x = element_text(size = 10, angle = 315, vjust = 0.2), legend.position = "none") +
geom_hline(yintercept = 1e-5) + geom_vline(xintercept = 2) 

```




![png](output_32_1.png)


<div class="alert alert-block alert-warning">
<b>Commentary:</b> You might notice that I only chose to apply relative abundance and prevalence filters. 
While I think that all three have their associated merits, CV calculations and filters are particularly 
susceptible to compositional effects. <b>Prior to estimation and application of CV parameters, estimate.ASthreshold() normalizes each sample to the median sequencing depth of the dataset, but this still does not properly account 
for the compositional nature of marker gene studies.</b> <i>With that said, you should feel free to apply 
whatever approaches that you feel are appropriate!</i>
</div>

Now that we've estimated our WS and AS filtering parameters, we can apply these thresholds to our dataset with the __filter.dataset()__ function. __filter.dataset()__ is a multifaceted filtering funtion that first applies our previously estimated WS filter, removes samples below a read count threshold and/or that match a metadata indentifier, and then applies one or multiple AS filters.

<div class="alert alert-block alert-info">
<b><font size="5">filter.dataset</font></b><br>
A highly flexible function that applies WS and AS (relative abundance, prevalence, CV), as well as metadata (positive and negative) and read count based filters. <br>
<br>
<b>Input:</b> <br>
<b>ps:</b> a phyloseq-class dataset.<br>
<br>
<b>Parameters (all are optional):</b><br>
<b>controlID:</b> As above, used only for reporting ASV read count, taxonomy, and sequence identity of taxa in the positive control.<br>
<b>mdCAT:</b> The sample data category used in conjunction with mdFACTOR for metadata based filtering.<br>
<b>mdFACTOR:</b> The factor within mdCAT used for metadata based filtering.<br>
<b>mdNEGATIVE:</b> A logical flag indicating whether metadata-based filtering should be positive (remove matching samples) or negative (remove non matching samples).<br>
<b>minLIB:</b> The minimum sequence depth for a sample to be retained.<br>
<b>RAT:</b> Relative Abundance Threshold; the percentage of reads across the dataset for which an ASV must be present for that ASV to be retained. <br>
<b>CVT:</b> Coefficient of Variation Threshold; the lowest CV value of an ASV for that ASV to be retained.<br>
<b>PFT:</b> Prevalence Filter Threshold; the percentage of samples across the dataset for which an ASV must be present for that ASV to be retained. <br>
<b>return.all:</b> A logical flag indicating whether to return all output or only the filtered dataset. Default is FALSE.<br>
<br><br>
<b>Output:</b> <br>
<b>control.taxa.sequences</b>
-  These are the nucleotide sequences of the taxa remaining in the positive control after applying the within-sample filter threshold. Useful for cross referenceing annotations against additional databases.<br>
<b>taxonomy.of.control.taxa</b>
-   The taxonomy of taxa in the positive control are output in the same order as the sequences, above. Generally useful for validating that the bugs detected in your mock community are as intended. <br>
<b>read.count.table</b>
-  This table lists the read count and percentage remaining after each filtering step. Samples matching given metadata identifiers are removed after application of the WS filter and prior to the AS filter(s) and, as such, will list NA values for AS filter statistics.  <br>
<b>relative abundance.filter.threshold</b>
-  The read count corresponding to the percentage given to the RAT parameter. <br>
<b>prevalence.filter.threshold</b>
-  The sample count corresponding to the percentage given to the PFT parameter. <br>
<b>filtered.phyloseq</b>
-  the shiny new (phyloseq class) filtered dataset with control samples removed. <br>

</div>



```R
NI.f <- filter.dataset(ps = NI, controlID = "pc1", mdCAT = "timepoint", 
                      mdFACTOR = "control", WST = 1.4e-4, PFT = 0.02, RAT = 1e-5, minLIB = 4000, return.all = TRUE)
```

    Removing 13 samples with read count < 4000
    Applying WS filter threshold of 0.00014
    Removing 6 samples matching metadata identifiers timepoint:control
    Applying relative abundance threshold of 1e-05
    Applying prevalence threshold of 0.02



```R
NI.f$ntaxa.in.control
```


25



```R
NI.f$control.taxa.sequences
```


<ol class=list-inline>
	<li>'GGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCCGATGCTATTTCTAGAGATAGGAAGTTACTTCGGTACATCGGTGACAGG'</li>
	<li>'GGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACAGAACTTTCCAGAGATGGATTGGTGCCTTCGGGAACTGTGAGACAGG'</li>
	<li>'GGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGG'</li>
	<li>'GGACCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAAATCTTGACATCCTTTGACAACTCTAGAGATAGAGCCTTCCCCTTCGGGGGACAAAGTGACA'</li>
	<li>'GGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCTTTGACCACTCTAGAGATAGAGCTTTCCCTTCGGGGACAAAGTGACAGG'</li>
	<li>'GGACCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAAATCTTGACATCCTCTGACCCCTCTAGAGATAGAGTTTTCCCCTTCGGGGGACAGAGTGACA'</li>
	<li>'GGACCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAAATCTTGACATCCTTTGACAACTCTAGAGATAGAGCTTTCCCCTTCGGGGGACAAAGTGACA'</li>
	<li>'GGGCCCGCACAAGCGGCGGAGCATGCGGATTAATTCGATGCAACGCGAAGAACCTTACCAAGGCTTGACATGCACGGCGGCACTGCAGAGATGTGGTGGCATTTAGTTGGTCGTGTGCAGGT'</li>
	<li>'GGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCCTCTGACCGCTCTAGAGATAGAGCTTTCCTTCGGGACAGAGGTGACAGG'</li>
	<li>'GGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCCGATGCTATTCTTAGAGATAGGAAGTTACTTCGGTACATCGGAGACAGG'</li>
	<li>'GGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGCCTTGACATCCAATGAACTTTCTAGAGATAGATTGGTGCCTTCGGGAACATTGAGACAGG'</li>
	<li>'GGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGCAGAACCTTACCAACCCTTGACATGGCGATCGCGGTTCCAGAGATGGTTCCTTCAGTTCGGCTGGATCGCACACA'</li>
	<li>'GGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCTTCTGACCGGCCTAGAGATAGGCTTTCTCTTCGGAGCAGAAGTGACAGG'</li>
	<li>'GGGCCCGCACAAGCGGAGGAACATGTGGTTTAATTCGATGATACGCGAGGAACCTTACCCGGGCTTAAATTGCAGATGAATTACGGTGAAAGCCGTAAGCCGCAAGGCATCTGTGAAGGTGC'</li>
	<li>'GGACCCGCACAAGCGGTGGATGATGTGGATTAATTCGATGCAACGCGAAAAACCTTACCTACCCTTGACATGTCTGGAATGCCGAAGAGATTTGGTAGTGCTCGCAAGAGAACCGGAACACA'</li>
	<li>'GGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCAGTGCAAACCTAAGAGATTAGGTGTTCCCTTCGGGGACGCTGAGACAGG'</li>
	<li>'GGCCCCGCACAAGCGGCGGAGCATGCGGATTAATTCGATGCAACGCGTAGAACCTTACCTGGGTTTGACATGGATCGGGAGTGCTCAGAGATGGGTGTGCCTCTTTTGGGGTCGGTTCACAG'</li>
	<li>'GGACCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAAATCTTGACATCCTCTGATCCCTCTAGAGATAGAGGTTTCCCCTTCGGGGGACAGAGTGACA'</li>
	<li>'GGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCTCTGACAACCCTAGAGATAGGGCTTCTCCTTCGGGAGCAGAGTGACAGG'</li>
	<li>'GGGCCCGCACAAGCAGCGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCTAGACTTGACATCTCCTGAATTACCCTTAATCGGGGAAGCCCTTCGGGGCAGGAAGACAGGTG'</li>
	<li>'GGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCTGGCCTTGACATGCTGAGAACTTTCCAGAGATGGATTGGTGCCTTCGGGAACTCAGACACAGG'</li>
	<li>'GGACCCGCACAAGCGGTGGATGATGTGGATTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATGTACGGAATCCTCCGGAGACGGAGGAGTGCCTTCGGGAGCCGTAACACAGG'</li>
	<li>'GGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCTTTGACCACTCTGGAGACAGAGCTTTCCCTTCGGGGACAAAGTGACAGG'</li>
	<li>'GGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGCCTTGACATACTAGAAACTTTCCAGAGATGGATTGGTGCCTTCGGGAATCTAGATACAGG'</li>
	<li>'GGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATGCTAGGAACTTTGCAGAGATGCAGAGGTGCCCTTCGGGGAACCTAGACACA'</li>
</ol>




```R
NI.f$taxonomy.of.control.taxa
```


<table>
<thead><tr><th scope=col>Kingdom</th><th scope=col>Phylum</th><th scope=col>Class</th><th scope=col>Order</th><th scope=col>Family</th><th scope=col>Genus</th><th scope=col>Species</th></tr></thead>
<tbody>
	<tr><td>Bacteria                                    </td><td>Firmicutes                                  </td><td>Bacilli                                     </td><td>Lactobacillales                             </td><td>Streptococcaceae                            </td><td>Streptococcus                               </td><td>lactarius/salivarius                        </td></tr>
	<tr><td>Bacteria                                    </td><td>Proteobacteria                              </td><td>Gammaproteobacteria                         </td><td>Enterobacteriales                           </td><td>Enterobacteriaceae                          </td><td>Salmonella                                  </td><td>NA                                          </td></tr>
	<tr><td>Bacteria                                    </td><td>Proteobacteria                              </td><td>Gammaproteobacteria                         </td><td>Enterobacteriales                           </td><td>Enterobacteriaceae                          </td><td>Escherichia/Shigella                        </td><td>NA                                          </td></tr>
	<tr><td>Bacteria                                    </td><td>Firmicutes                                  </td><td>Bacilli                                     </td><td>Bacillales                                  </td><td>Staphylococcaceae                           </td><td>Staphylococcus                              </td><td>NA                                          </td></tr>
	<tr><td>Bacteria                                    </td><td>Firmicutes                                  </td><td>Bacilli                                     </td><td>Lactobacillales                             </td><td>Enterococcaceae                             </td><td>Enterococcus                                </td><td>NA                                          </td></tr>
	<tr><td>Bacteria                                    </td><td>Firmicutes                                  </td><td>Bacilli                                     </td><td>Bacillales                                  </td><td>Staphylococcaceae                           </td><td>Staphylococcus                              </td><td>NA                                          </td></tr>
	<tr><td>Bacteria                                    </td><td>Firmicutes                                  </td><td>Bacilli                                     </td><td>Bacillales                                  </td><td>Staphylococcaceae                           </td><td>Staphylococcus                              </td><td>NA                                          </td></tr>
	<tr><td>Bacteria                                    </td><td>Actinobacteria                              </td><td>Actinobacteria                              </td><td>Actinomycetales                             </td><td>Actinomycetaceae                            </td><td>Actinomyces                                 </td><td>lingnae/odontolyticus/turicensis            </td></tr>
	<tr><td>Bacteria                                    </td><td>Firmicutes                                  </td><td>Bacilli                                     </td><td>Lactobacillales                             </td><td>Streptococcaceae                            </td><td>Streptococcus                               </td><td>NA                                          </td></tr>
	<tr><td>Bacteria                                    </td><td>Firmicutes                                  </td><td>Bacilli                                     </td><td>Lactobacillales                             </td><td>Streptococcaceae                            </td><td>Streptococcus                               </td><td>mutans                                      </td></tr>
	<tr><td>Bacteria                                    </td><td>Proteobacteria                              </td><td>Gammaproteobacteria                         </td><td>Pseudomonadales                             </td><td>Pseudomonadaceae                            </td><td>Pseudomonas                                 </td><td>NA                                          </td></tr>
	<tr><td>Bacteria                                    </td><td>Proteobacteria                              </td><td>Alphaproteobacteria                         </td><td>Rhodobacterales                             </td><td>Rhodobacteraceae                            </td><td>Rhodobacter                                 </td><td>azotoformans/johrii/megalophilus/sphaeroides</td></tr>
	<tr><td>Bacteria                                    </td><td>Firmicutes                                  </td><td>Bacilli                                     </td><td>Lactobacillales                             </td><td>NA                                          </td><td>NA                                          </td><td>NA                                          </td></tr>
	<tr><td>Bacteria                                    </td><td>Bacteroidetes                               </td><td>Bacteroidia                                 </td><td>Bacteroidales                               </td><td>Bacteroidaceae                              </td><td>Bacteroides                                 </td><td>vulgatus                                    </td></tr>
	<tr><td>Bacteria                                    </td><td>Proteobacteria                              </td><td>Betaproteobacteria                          </td><td>Burkholderiales                             </td><td>Alcaligenaceae                              </td><td>Achromobacter                               </td><td>NA                                          </td></tr>
	<tr><td>Bacteria                                    </td><td>Firmicutes                                  </td><td>Bacilli                                     </td><td>Lactobacillales                             </td><td>Lactobacillaceae                            </td><td>Lactobacillus                               </td><td>NA                                          </td></tr>
	<tr><td>Bacteria                                    </td><td>Actinobacteria                              </td><td>Actinobacteria                              </td><td>Actinomycetales                             </td><td>Propionibacteriaceae                        </td><td>Propionibacterium                           </td><td>acnes                                       </td></tr>
	<tr><td>Bacteria                                    </td><td>Firmicutes                                  </td><td>Bacilli                                     </td><td>Bacillales                                  </td><td>NA                                          </td><td>NA                                          </td><td>NA                                          </td></tr>
	<tr><td>Bacteria                                    </td><td>Firmicutes                                  </td><td>Bacilli                                     </td><td>Bacillales                                  </td><td>Bacillaceae_1                               </td><td>Bacillus                                    </td><td>NA                                          </td></tr>
	<tr><td>Bacteria                                    </td><td>Firmicutes                                  </td><td>Clostridia                                  </td><td>Clostridiales                               </td><td>Clostridiaceae_1                            </td><td>NA                                          </td><td>NA                                          </td></tr>
	<tr><td>Bacteria                                    </td><td>Proteobacteria                              </td><td>Gammaproteobacteria                         </td><td>Pseudomonadales                             </td><td>Pseudomonadaceae                            </td><td>NA                                          </td><td>NA                                          </td></tr>
	<tr><td>Bacteria                                    </td><td>Proteobacteria                              </td><td>Betaproteobacteria                          </td><td>Neisseriales                                </td><td>Neisseriaceae                               </td><td>Neisseria                                   </td><td>NA                                          </td></tr>
	<tr><td>Bacteria                                    </td><td>Firmicutes                                  </td><td>Bacilli                                     </td><td>Bacillales                                  </td><td>Listeriaceae                                </td><td>Listeria                                    </td><td>NA                                          </td></tr>
	<tr><td>Bacteria                                    </td><td>Proteobacteria                              </td><td>Gammaproteobacteria                         </td><td>Pseudomonadales                             </td><td>Moraxellaceae                               </td><td>Acinetobacter                               </td><td>NA                                          </td></tr>
	<tr><td>Bacteria                                    </td><td>Deinococcus-Thermus                         </td><td>Deinococci                                  </td><td>Deinococcales                               </td><td>Deinococcaceae                              </td><td>Deinococcus                                 </td><td>radiodurans                                 </td></tr>
</tbody>
</table>




```R
NI.f$read.count.table[135:149, ]
```


<table>
<thead><tr><th></th><th scope=col>unfiltered.read.count</th><th scope=col>WSfiltered.read.count</th><th scope=col>WSfiltered.read.percent</th><th scope=col>ASfiltered.read.count</th><th scope=col>ASfiltered.read.percent</th></tr></thead>
<tbody>
	<tr><th scope=row>N.0395..Week.36</th><td>1256538  </td><td>1252944  </td><td> 99.71398</td><td>1212704  </td><td>96.51153 </td></tr>
	<tr><th scope=row>N.0395..Week.4</th><td> 523876  </td><td> 522239  </td><td> 99.68752</td><td> 521647  </td><td>99.57452 </td></tr>
	<tr><th scope=row>N.0704..baseline</th><td>  71953  </td><td>  71268  </td><td> 99.04799</td><td>  69124  </td><td>96.06827 </td></tr>
	<tr><th scope=row>N.0704..Week.15</th><td> 782531  </td><td> 780733  </td><td> 99.77023</td><td> 779930  </td><td>99.66762 </td></tr>
	<tr><th scope=row>N.0704..Week.36</th><td> 886216  </td><td> 883769  </td><td> 99.72388</td><td> 874150  </td><td>98.63848 </td></tr>
	<tr><th scope=row>N.0704..Week.4</th><td> 177996  </td><td> 176089  </td><td> 98.92863</td><td> 173717  </td><td>97.59601 </td></tr>
	<tr><th scope=row>N.0705..baseline</th><td> 888231  </td><td> 885186  </td><td> 99.65718</td><td> 866188  </td><td>97.51833 </td></tr>
	<tr><th scope=row>N.0705..Week.15</th><td>1193231  </td><td>1188586  </td><td> 99.61072</td><td>1188139  </td><td>99.57326 </td></tr>
	<tr><th scope=row>N.0705..Week.4</th><td> 951118  </td><td> 947463  </td><td> 99.61572</td><td> 947129  </td><td>99.58060 </td></tr>
	<tr><th scope=row>nc1</th><td> 111821  </td><td> 110452  </td><td> 98.77572</td><td>     NA  </td><td>      NA </td></tr>
	<tr><th scope=row>nc2</th><td>  10050  </td><td>  10050  </td><td>100.00000</td><td>     NA  </td><td>      NA </td></tr>
	<tr><th scope=row>nc3</th><td> 178762  </td><td> 176944  </td><td> 98.98301</td><td>     NA  </td><td>      NA </td></tr>
	<tr><th scope=row>nc4</th><td>  16468  </td><td>  16458  </td><td> 99.93928</td><td>     NA  </td><td>      NA </td></tr>
	<tr><th scope=row>nc5</th><td>  38186  </td><td>  37955  </td><td> 99.39507</td><td>     NA  </td><td>      NA </td></tr>
	<tr><th scope=row>pc1</th><td> 745294  </td><td> 744935  </td><td> 99.95183</td><td>     NA  </td><td>      NA </td></tr>
</tbody>
</table>




```R
NI.f$relative.abundance.filter.read.count
```


915.10988



```R
NI.f$prevalence.filter.sample.count
```


2.86



```R
(NIF <- NI.f$filtered.phyloseq)
```


    phyloseq-class experiment-level object
    otu_table()   OTU Table:         [ 344 taxa and 143 samples ]
    sample_data() Sample Data:       [ 143 samples by 9 sample variables ]
    tax_table()   Taxonomy Table:    [ 344 taxa by 7 taxonomic ranks ]


## Processing and exporting
<br>Now that we've filtered the dataset to our liking, we'll export all of the data associated with our final, filtered dataset using the __write.dataset()__ function. Unfortunately, I've never been able to get the biomformat package (v1.2.0) to play nice with my version of phyloseq, but y'all may have more luck with newer versions. In case the kinks have been worked out, I'd recommend using the __write.dataset.biom()__ function to export your dataset directly to biom format. Until then, we'll use __write.dataset()__ to perform some final minor modifications and export our dataset. 

<div class="alert alert-block alert-info">
<b><font size="5">write.dataset, write.dataset.biom</font></b><br>
Exporting functions for writing marker gene datasets to external files.<br>
<br>
<b>Input:</b> <br>
<b>ps:</b> a phyloseq-class dataset.<br>
<br>
<b>Parameters:</b><br>
<b>Filepath:</b> The filepath where output should be written.<br>
<b>Fileprefix:</b> The prefix identifier for all output written.<br>
<b>Rename:</b> A logical flag indicating if ASV identifiers should be renamed to ASV1, ASV2, ..., ASVN. Default is FALSE.<br>
<br>
<b>Output:</b> <br>
<b>_ASVs.fasta:</b>  FASTA of all ASV sequences. Even if rename=FALSE, FASTA identifiers are labeled ASV1, ASV2,... <br>
<br>
<b>write.dataset.biom:</b><br>
<b>_ASV_table.biom:</b> A biom file with observation and sample metadata integrated. I've run into problems attempting to import the resulting biom back into phyloseq, but this might warrant a closer look with a newer version. <br>
<br>
<b>write.dataset:</b><br>
<b>_ASV_table.txt:</b>  A standard count table (eg. OTU table) of read counts of taxa (ASVs) across all samples. <br>
<b>_ASV_taxonomy.txt:</b> A taxonomy table with the sequence identifier (eg. ASV1) in the first column and semicolon-separated taxonomy in the second column.<br>
<b>_sample_data.txt:</b> A sample data table of all associated metadata. 
<br>

</div>
<br>
<div class="alert alert-block alert-success">
<b>TIP:</b> Before writing, all files written by <b>write.dataset()</b> are formatted for easy merging into a biom file using the  <A HREF = "http://biom-format.org/documentation/biom_conversion.html">standalone biom package</A> without any further manipulation. 
</div>


```R
fp <- paste(getwd(), "/", sep = "")
NIF.f <- write.dataset(ps = NIF, filepath = fp, fileprefix = "Nif.f", rename = TRUE)
```

    [1] "/Users/bpb/Desktop/marker_gene_filtering_functions/marker_gene_processing_vignette/Nif.f_ASVs.fasta"
    [1] "/Users/bpb/Desktop/marker_gene_filtering_functions/marker_gene_processing_vignette/Nif.f_ASV_table.txt"
    [1] "/Users/bpb/Desktop/marker_gene_filtering_functions/marker_gene_processing_vignette/Nif.f_ASV_taxonomy.txt"
    [1] "/Users/bpb/Desktop/marker_gene_filtering_functions/marker_gene_processing_vignette/Nif.f_sample_data.txt"


Boom! We've now written all of the metadata associated with our filtered marker gene study to various output files. Now, let's take a final peek at our ASV count table in a format that's a little more friendly on the eyes. 


```R
phyloseq::otu_table(NIF.f)[1:10, 1:10]
```


<table>
<thead><tr><th></th><th scope=col>ASV1</th><th scope=col>ASV2</th><th scope=col>ASV3</th><th scope=col>ASV4</th><th scope=col>ASV5</th><th scope=col>ASV6</th><th scope=col>ASV7</th><th scope=col>ASV8</th><th scope=col>ASV9</th><th scope=col>ASV10</th></tr></thead>
<tbody>
	<tr><th scope=row>N.0055..baseline</th><td>     0</td><td> 57092</td><td>    0 </td><td>   0  </td><td>    0 </td><td>    0 </td><td>532461</td><td>  8519</td><td>     0</td><td>55386 </td></tr>
	<tr><th scope=row>N.0055..Week.15</th><td>  3137</td><td>    50</td><td>  126 </td><td> 170  </td><td>   12 </td><td>   16 </td><td>    59</td><td>    26</td><td>    27</td><td>    0 </td></tr>
	<tr><th scope=row>N.0055..Week.36</th><td>  1745</td><td>583068</td><td>  295 </td><td>7571  </td><td>    0 </td><td>  117 </td><td>     0</td><td>   225</td><td>161094</td><td>    0 </td></tr>
	<tr><th scope=row>N.0055..Week.4</th><td>  1552</td><td>  4311</td><td>24642 </td><td>3267  </td><td>  422 </td><td>    0 </td><td>    86</td><td>104530</td><td>     0</td><td> 1295 </td></tr>
	<tr><th scope=row>N.0170..baseline</th><td>  1438</td><td>   368</td><td>  362 </td><td> 772  </td><td>  431 </td><td>  392 </td><td>  1023</td><td>   361</td><td>   265</td><td> 1382 </td></tr>
	<tr><th scope=row>N.0170.Week.15</th><td>136177</td><td>   269</td><td>99767 </td><td>   0  </td><td>  642 </td><td>  505 </td><td>  1201</td><td>    72</td><td>  3606</td><td>   91 </td></tr>
	<tr><th scope=row>N.0170.Week.36</th><td>367401</td><td>  8053</td><td>17969 </td><td>   0  </td><td>40399 </td><td>31963 </td><td>   271</td><td>     0</td><td> 50041</td><td>  752 </td></tr>
	<tr><th scope=row>N.0170.Week.4</th><td>335868</td><td>   788</td><td>77043 </td><td>   0  </td><td>11923 </td><td> 6800 </td><td>  1818</td><td> 42446</td><td>  4513</td><td>    0 </td></tr>
	<tr><th scope=row>N.0183..baseline</th><td> 12981</td><td>  5611</td><td>44012 </td><td>1452  </td><td> 1553 </td><td> 1915 </td><td>  1577</td><td>  1821</td><td>   231</td><td> 2531 </td></tr>
	<tr><th scope=row>N.0183.Week.15</th><td> 58948</td><td>    97</td><td>23709 </td><td>  24  </td><td>  111 </td><td>    0 </td><td>     0</td><td>    56</td><td>  2227</td><td>    0 </td></tr>
</tbody>
</table>



Finally, let's subset our filtered dataset to include only samples that were collected during the first week of the study. We can enable negative metadata filtering with the flag mdNEGATIVE.


```R
NIF.fw1 <- filter.dataset(ps = NIF.f, mdCAT = "timepoint", mdFACTOR = "Week 1", mdNEGATIVE = TRUE, return.all = FALSE)
```

    Removing 111 samples not matching metadata identifiers timepoint:Week 1



```R
NIF.fw1
```


    phyloseq-class experiment-level object
    otu_table()   OTU Table:         [ 344 taxa and 32 samples ]
    sample_data() Sample Data:       [ 32 samples by 9 sample variables ]
    tax_table()   Taxonomy Table:    [ 344 taxa by 7 taxonomic ranks ]

