
<font size=7>Filtering marker gene datasets</font>

Here, we'll walk through some filtering strategies applicable to compositional marker gene studies. We will employ several custom R scripts designed to manipulate phyloseq-class datasets that aid in the filtering and reproducability of marker gene studies. We apply this workflow to a subset of a larger 16S community dataset, though these scripts are scalable to any marker gene and large datasets with thousands of taxa and samples. 

<div class="alert alert-block alert-info">
<font size=3><b>A preface on filtering:</b></font> Filtering is an <u>essential</u> step in the analyis of marker gene datasets. This workflow uses a two-pronged approach to filtering contaminants from read datasets. That is, we attempt to remove two types of contamination: <br>
<br>
1. Real biological sequences that have cross contaminated unassociated samples. <br>
<br>
2. Exogenous contaminants introduced during sample preparation. <br>
<br>
To achieve this, we'll first filter out taxa below an abundance threshold individually within each sample. We estimate this threshold empirically by measuring the number of amplicon sequence variants (ASVs) within known mock communities under various thresholds and select the threshold that returns the desired results. Next, we'll remove control samples from the dataset and then apply a dataset-wide prevalence filter to remove exogenous contaminants present as low abundance ASVs. 
</div>

***

<font size=5><b>Before we begin...</b></font><br>
This workflow assumes that your sequencing project has been:
-  processed with __[dada2](https://benjjneb.github.io/dada2/index.html)__ (namely, that your ASV names are the ASV nucleotide sequences)
-  converted into a __[phyloseq](https://joey711.github.io/phyloseq/)__ object


```R
library(RCurl); packageVersion("RCurl")
```


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


    phyloseq-class experiment-level object
    otu_table()   OTU Table:         [ 2195 taxa and 149 samples ]
    sample_data() Sample Data:       [ 149 samples by 9 sample variables ]
    tax_table()   Taxonomy Table:    [ 2195 taxa by 7 taxonomic ranks ]



```R
(ws.threshold <- estimate.threshold(ps = NI, WSmin = 1e-4, WSmax = 2e-4, WSstep = 1e-5, controlID = "pc1"))
```


<table>
<thead><tr><th scope=col>control.taxa.count</th><th scope=col>read.count</th><th scope=col>read.percent</th><th scope=col>threshold.value</th></tr></thead>
<tbody>
	<tr><td>23      </td><td>89567283</td><td>99.80368</td><td>0.00010 </td></tr>
	<tr><td>23      </td><td>89553735</td><td>99.78859</td><td>0.00011 </td></tr>
	<tr><td>23      </td><td>89541319</td><td>99.77475</td><td>0.00012 </td></tr>
	<tr><td>23      </td><td>89527748</td><td>99.75963</td><td>0.00013 </td></tr>
	<tr><td>23      </td><td>89513016</td><td>99.74321</td><td>0.00014 </td></tr>
	<tr><td>23      </td><td>89500419</td><td>99.72918</td><td>0.00015 </td></tr>
	<tr><td>22      </td><td>89487204</td><td>99.71445</td><td>0.00016 </td></tr>
	<tr><td>22      </td><td>89472732</td><td>99.69833</td><td>0.00017 </td></tr>
	<tr><td>22      </td><td>89460411</td><td>99.68460</td><td>0.00018 </td></tr>
	<tr><td>22      </td><td>89449162</td><td>99.67206</td><td>0.00019 </td></tr>
	<tr><td>20      </td><td>89437647</td><td>99.65923</td><td>0.00020 </td></tr>
</tbody>
</table>



And let's take a peek at a quick visual summary


```R
ggplot(data = ws.threshold, aes(x = threshold.value, y = control.taxa.count)) + geom_line(size=2, color="steelblue2") + 
  theme(axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"))

```




![png](output_12_1.png)



```R
ggplot(data = ws.threshold, aes(x = threshold.value, y = read.percent)) + geom_line(size=2, color="orangered1") + 
  theme(axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"))

```




![png](output_13_1.png)


***
It looks like a within-sample threshold of 2e-4 is sufficient for recovering the desired taxa count for our mock communities, without sacrificing too many reads. Now that we've determined our threshold, let's apply this cutoff to our dataset with the __filter.dataset()__ function


```R
NI.f <- filter.dataset(ps = NI, controlID = "pc1", controlCAT = "timepoint", 
                      controlFACTOR = "control", WSF = 1e-3, PF = 1e-5, return.all = TRUE)
```


```R
NI.f$control.taxa.sequences
```


<ol class=list-inline>
	<li>'GGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACAGAACTTTCCAGAGATGGATTGGTGCCTTCGGGAACTGTGAGACAGG'</li>
	<li>'GGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGG'</li>
	<li>'GGACCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAAATCTTGACATCCTTTGACAACTCTAGAGATAGAGCCTTCCCCTTCGGGGGACAAAGTGACA'</li>
	<li>'GGACCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAAATCTTGACATCCTCTGACCCCTCTAGAGATAGAGTTTTCCCCTTCGGGGGACAGAGTGACA'</li>
	<li>'GGACCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAAATCTTGACATCCTTTGACAACTCTAGAGATAGAGCTTTCCCCTTCGGGGGACAAAGTGACA'</li>
	<li>'GGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCCGATGCTATTCTTAGAGATAGGAAGTTACTTCGGTACATCGGAGACAGG'</li>
	<li>'GGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGCAGAACCTTACCAACCCTTGACATGGCGATCGCGGTTCCAGAGATGGTTCCTTCAGTTCGGCTGGATCGCACACA'</li>
	<li>'GGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCAGTGCAAACCTAAGAGATTAGGTGTTCCCTTCGGGGACGCTGAGACAGG'</li>
	<li>'GGCCCCGCACAAGCGGCGGAGCATGCGGATTAATTCGATGCAACGCGTAGAACCTTACCTGGGTTTGACATGGATCGGGAGTGCTCAGAGATGGGTGTGCCTCTTTTGGGGTCGGTTCACAG'</li>
	<li>'GGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCTCTGACAACCCTAGAGATAGGGCTTCTCCTTCGGGAGCAGAGTGACAGG'</li>
	<li>'GGGCCCGCACAAGCAGCGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCTAGACTTGACATCTCCTGAATTACCCTTAATCGGGGAAGCCCTTCGGGGCAGGAAGACAGGTG'</li>
	<li>'GGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCTGGCCTTGACATGCTGAGAACTTTCCAGAGATGGATTGGTGCCTTCGGGAACTCAGACACAGG'</li>
	<li>'GGACCCGCACAAGCGGTGGATGATGTGGATTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATGTACGGAATCCTCCGGAGACGGAGGAGTGCCTTCGGGAGCCGTAACACAGG'</li>
	<li>'GGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCTTTGACCACTCTGGAGACAGAGCTTTCCCTTCGGGGACAAAGTGACAGG'</li>
</ol>




```R
NI.f$taxonomy.of.control.taxa
```


<table>
<thead><tr><th scope=col>Kingdom</th><th scope=col>Phylum</th><th scope=col>Class</th><th scope=col>Order</th><th scope=col>Family</th><th scope=col>Genus</th><th scope=col>Species</th></tr></thead>
<tbody>
	<tr><td>Bacteria                                    </td><td>Proteobacteria                              </td><td>Gammaproteobacteria                         </td><td>Enterobacteriales                           </td><td>Enterobacteriaceae                          </td><td>Salmonella                                  </td><td>NA                                          </td></tr>
	<tr><td>Bacteria                                    </td><td>Proteobacteria                              </td><td>Gammaproteobacteria                         </td><td>Enterobacteriales                           </td><td>Enterobacteriaceae                          </td><td>Escherichia/Shigella                        </td><td>NA                                          </td></tr>
	<tr><td>Bacteria                                    </td><td>Firmicutes                                  </td><td>Bacilli                                     </td><td>Bacillales                                  </td><td>Staphylococcaceae                           </td><td>Staphylococcus                              </td><td>NA                                          </td></tr>
	<tr><td>Bacteria                                    </td><td>Firmicutes                                  </td><td>Bacilli                                     </td><td>Bacillales                                  </td><td>Staphylococcaceae                           </td><td>Staphylococcus                              </td><td>NA                                          </td></tr>
	<tr><td>Bacteria                                    </td><td>Firmicutes                                  </td><td>Bacilli                                     </td><td>Bacillales                                  </td><td>Staphylococcaceae                           </td><td>Staphylococcus                              </td><td>NA                                          </td></tr>
	<tr><td>Bacteria                                    </td><td>Firmicutes                                  </td><td>Bacilli                                     </td><td>Lactobacillales                             </td><td>Streptococcaceae                            </td><td>Streptococcus                               </td><td>mutans                                      </td></tr>
	<tr><td>Bacteria                                    </td><td>Proteobacteria                              </td><td>Alphaproteobacteria                         </td><td>Rhodobacterales                             </td><td>Rhodobacteraceae                            </td><td>Rhodobacter                                 </td><td>azotoformans/johrii/megalophilus/sphaeroides</td></tr>
	<tr><td>Bacteria                                    </td><td>Firmicutes                                  </td><td>Bacilli                                     </td><td>Lactobacillales                             </td><td>Lactobacillaceae                            </td><td>Lactobacillus                               </td><td>NA                                          </td></tr>
	<tr><td>Bacteria                                    </td><td>Actinobacteria                              </td><td>Actinobacteria                              </td><td>Actinomycetales                             </td><td>Propionibacteriaceae                        </td><td>Propionibacterium                           </td><td>acnes                                       </td></tr>
	<tr><td>Bacteria                                    </td><td>Firmicutes                                  </td><td>Bacilli                                     </td><td>Bacillales                                  </td><td>Bacillaceae_1                               </td><td>Bacillus                                    </td><td>NA                                          </td></tr>
	<tr><td>Bacteria                                    </td><td>Firmicutes                                  </td><td>Clostridia                                  </td><td>Clostridiales                               </td><td>Clostridiaceae_1                            </td><td>NA                                          </td><td>NA                                          </td></tr>
	<tr><td>Bacteria                                    </td><td>Proteobacteria                              </td><td>Gammaproteobacteria                         </td><td>Pseudomonadales                             </td><td>Pseudomonadaceae                            </td><td>NA                                          </td><td>NA                                          </td></tr>
	<tr><td>Bacteria                                    </td><td>Proteobacteria                              </td><td>Betaproteobacteria                          </td><td>Neisseriales                                </td><td>Neisseriaceae                               </td><td>Neisseria                                   </td><td>NA                                          </td></tr>
	<tr><td>Bacteria                                    </td><td>Firmicutes                                  </td><td>Bacilli                                     </td><td>Bacillales                                  </td><td>Listeriaceae                                </td><td>Listeria                                    </td><td>NA                                          </td></tr>
</tbody>
</table>




```R
NI.f$read.count.table[142:149, ]
```


<table>
<thead><tr><th></th><th scope=col>unfiltered.read.count</th><th scope=col>WSfiltered.read.count</th><th scope=col>WSfiltered.read.percent</th><th scope=col>Pfiltered.read.count</th><th scope=col>Pfiltered.read.percent</th></tr></thead>
<tbody>
	<tr><th scope=row>N.0705..Week.15</th><td>1163553 </td><td>1142012 </td><td>98.14869</td><td>1142012 </td><td>98.14869</td></tr>
	<tr><th scope=row>N.0705..Week.4</th><td> 934760 </td><td> 923799 </td><td>98.82740</td><td> 923799 </td><td>98.82740</td></tr>
	<tr><th scope=row>nc1</th><td> 109438 </td><td> 105591 </td><td>96.48477</td><td>     NA </td><td>      NA</td></tr>
	<tr><th scope=row>nc2</th><td>   9943 </td><td>   9787 </td><td>98.43106</td><td>     NA </td><td>      NA</td></tr>
	<tr><th scope=row>nc3</th><td> 177463 </td><td> 172042 </td><td>96.94528</td><td>     NA </td><td>      NA</td></tr>
	<tr><th scope=row>nc4</th><td>  15991 </td><td>  15638 </td><td>97.79251</td><td>     NA </td><td>      NA</td></tr>
	<tr><th scope=row>nc5</th><td>  37114 </td><td>  35909 </td><td>96.75325</td><td>     NA </td><td>      NA</td></tr>
	<tr><th scope=row>pc1</th><td> 697410 </td><td> 694799 </td><td>99.62561</td><td>     NA </td><td>      NA</td></tr>
</tbody>
</table>




```R
NI.f$prevalence.filter.threshold
```


877.30769



```R
(NIF <- NI.f$filtered.phyloseq)
```


    phyloseq-class experiment-level object
    otu_table()   OTU Table:         [ 322 taxa and 143 samples ]
    sample_data() Sample Data:       [ 143 samples by 9 sample variables ]
    tax_table()   Taxonomy Table:    [ 322 taxa by 7 taxonomic ranks ]


Next, we'll export all of the data associated with our final, filtered dataset using the __write.dataset()__ function. Unfortunately, I've never been able to get the biomformat package (v1.2.0) to play nice with my version of phyloseq, but y'all may have more luck with newer versions. In case the kinks have been worked out, I'd recommend using the __write.dataset.biom()__ function to export your dataset directly to biom format. Until then, we'll use __write.dataset()__ to perform some final minor modifications and export our dataset.

<div class="alert alert-block alert-info">
<font size=3><b>Output:</b></font> The output of these scripts are as follows: <br>
<br>
<font size=3><b>Both:</b></font><br>
-  FASTA of all ASV sequences. dada2 does not assign names to each ASV, but rather uses the nucleotide sequence as the identifier; we can leverage this feature to write an output in FASTA format of each ASV before we assign names.  <br>
<br>
<font size=3><b>write.dataset.biom()</b></font><br>
-   A biom file with observation and sample metadata integrated. I've run into problems attempting to import the resulting biom back into phyloseq, but this might warrant a closer look with a newer version. <br>
<br>
<font size=3><b>write.dataset()</b></font><br>
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

    [1] "/Users/bpb/Desktop/marker_gene_filtering_functions/Nif.f_ASVs.fasta"
    [1] "/Users/bpb/Desktop/marker_gene_filtering_functions/Nif.f_ASV_table.txt"
    [1] "/Users/bpb/Desktop/marker_gene_filtering_functions/Nif.f_ASV_taxonomy.txt"
    [1] "/Users/bpb/Desktop/marker_gene_filtering_functions/Nif.f_sample_data.txt"



```R
phyloseq::otu_table(NIF.f)[1:10, 1:10]
```


<table>
<thead><tr><th></th><th scope=col>ASV1</th><th scope=col>ASV2</th><th scope=col>ASV3</th><th scope=col>ASV4</th><th scope=col>ASV5</th><th scope=col>ASV6</th><th scope=col>ASV7</th><th scope=col>ASV8</th><th scope=col>ASV9</th><th scope=col>ASV10</th></tr></thead>
<tbody>
	<tr><th scope=row>N.0055..baseline</th><td>     0</td><td> 57092</td><td>    0 </td><td>   0  </td><td>    0 </td><td>    0 </td><td>532461</td><td>  8519</td><td>     0</td><td>55386 </td></tr>
	<tr><th scope=row>N.0055..Week.15</th><td>  3137</td><td>    50</td><td>  126 </td><td> 170  </td><td>   12 </td><td>   16 </td><td>    59</td><td>    26</td><td>    27</td><td>    0 </td></tr>
	<tr><th scope=row>N.0055..Week.36</th><td>  1745</td><td>583068</td><td>    0 </td><td>7571  </td><td>    0 </td><td>    0 </td><td>     0</td><td>     0</td><td>161094</td><td>    0 </td></tr>
	<tr><th scope=row>N.0055..Week.4</th><td>  1552</td><td>  4311</td><td>24642 </td><td>3267  </td><td>    0 </td><td>    0 </td><td>     0</td><td>104530</td><td>     0</td><td> 1295 </td></tr>
	<tr><th scope=row>N.0170..baseline</th><td>  1438</td><td>   368</td><td>  362 </td><td> 772  </td><td>  431 </td><td>  392 </td><td>  1023</td><td>   361</td><td>   265</td><td> 1382 </td></tr>
	<tr><th scope=row>N.0170.Week.15</th><td>136177</td><td>     0</td><td>99767 </td><td>   0  </td><td>  642 </td><td>  505 </td><td>  1201</td><td>     0</td><td>  3606</td><td>    0 </td></tr>
	<tr><th scope=row>N.0170.Week.36</th><td>367401</td><td>  8053</td><td>17969 </td><td>   0  </td><td>40399 </td><td>31963 </td><td>     0</td><td>     0</td><td> 50041</td><td>  752 </td></tr>
	<tr><th scope=row>N.0170.Week.4</th><td>335868</td><td>     0</td><td>77043 </td><td>   0  </td><td>11923 </td><td> 6800 </td><td>  1818</td><td> 42446</td><td>  4513</td><td>    0 </td></tr>
	<tr><th scope=row>N.0183..baseline</th><td> 12981</td><td>  5611</td><td>44012 </td><td>1452  </td><td> 1553 </td><td> 1915 </td><td>  1577</td><td>  1821</td><td>     0</td><td> 2531 </td></tr>
	<tr><th scope=row>N.0183.Week.15</th><td> 58948</td><td>     0</td><td>23709 </td><td>   0  </td><td>    0 </td><td>    0 </td><td>     0</td><td>     0</td><td>  2227</td><td>    0 </td></tr>
</tbody>
</table>


