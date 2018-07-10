##functions for processing compositional microbiome data
#all scripts written by BPB, 070418


estimate.threshold <- function(ps, Tmin=1e-4, Tmax=2e-4, Tstep=1e-5, controlID) {
  
  #Build param lists
  l.t <- seq(from = Tmin, to = Tmax, by = Tstep)
  nt <- length(l.t)
  tvec <- c()
  svec <- c()
  
  for (i in 1:nt){
    tryCatch({
      #loop through values
      filterfx = function(x){
        x[(x / sum(x)) < (l.t[i])] <- 0
        return(x)
      }
      
      ps.if <- phyloseq::transform_sample_counts(ps, fun = filterfx)
      
      #FILTERING
      if(as.logical(class(phyloseq::otu_table(ps.if))[1] == "otu_table") && 
         as.logical(taxa_are_rows(phyloseq::otu_table(ps.if)) == TRUE)){
        otu.tab <- as.matrix(phyloseq::otu_table(ps.if))
      } else {
        otu.tab <- as.matrix(t(phyloseq::otu_table(ps.if)))
      }
      
      tvec[i] <- nrow(otu.tab[which(otu.tab[,match(controlID, colnames(otu.tab))] != 0),])
      
      svec[i] <- sum(phyloseq::sample_sums(ps.if))
                     
    },
    error=function(e){cat("Warning :",conditionMessage(e), "\n")})
  }
  names(tvec) <- c(l.t)
  names(svec) <- c(l.t)
  
  # Build return list
  l.return = list()
    l.return[['ASVs in control at IF threshold']] <- tvec
    l.return[['reads across dataset retained at IF threshold']] <- svec
   return(l.return)
  
}


filter.dataset <- function(ps, controlID=NULL, controlCAT=NULL, controlFACTOR=NULL, TF=1e-4, PF=1e-5, return.all=TRUE){
  
  #per-sample filter function
  filterfx = function(x){
    x[(x / sum(x)) < TF] <- 0
    return(x)
  }
  #create unfiltered sample sum vector
  ov <- phyloseq::sample_sums(ps)
  
  #create filtered object
  ps.if <- phyloseq::transform_sample_counts(ps, fun = filterfx)
  sampledf <- suppressWarnings(as.matrix(phyloseq::sample_data(ps.if)))
  #create per sample filtered sample sum vector
  ifv <- phyloseq::sample_sums(ps.if)
  #calculate percent filtered, individual
  p.if <- phyloseq::sample_sums(ps.if)/phyloseq::sample_sums(ps)*100
  
  #prevalence filtering
  if(as.logical(class(phyloseq::otu_table(ps.if))[1] == "otu_table") && 
     as.logical(taxa_are_rows(phyloseq::otu_table(ps.if)) == TRUE)){
    otu.tab <- as.matrix(phyloseq::otu_table(ps.if))
  } else {
    otu.tab <- as.matrix(t(phyloseq::otu_table(ps.if)))
  }
  npos <- nrow(otu.tab[which(otu.tab[,match(controlID, colnames(otu.tab))] != 0),])
  #taxonomy of taxa in positive control
  tax.tab <- phyloseq::tax_table(ps.if)
  taxanames.control <- rownames(otu.tab[which(otu.tab[,match(controlID, colnames(otu.tab))] != 0),])
  tax.tab.subset <- tax.tab[taxanames.control] #taxonomy of taxa in positive control
  ttsn <- tax.tab.subset
  rownames(ttsn) <- NULL
  
  #remove controls
  filtered.names <- rownames(sampledf[which(sampledf[,match(controlCAT, colnames(sampledf))] != controlFACTOR),])
  sampledf.s <- as.data.frame(sampledf[filtered.names,])
  phyloseq::sample_data(ps.if) <- phyloseq::sample_data(sampledf.s)
  
  #dataset-wide prevalence filter %
  #determine taxa sums for cutoff %
  prevf <- sum(phyloseq::taxa_sums(ps.if)) * PF
  ps.cf <- phyloseq::prune_taxa(taxa_sums(ps.if)>=prevf, ps.if)
  #create prevalence filter sample sum vector
  pfv <- phyloseq::sample_sums(ps.cf)
  #calculate percent filtered, prevalence
  p.pf <- suppressWarnings(phyloseq::sample_sums(ps.cf)/phyloseq::sample_sums(ps)[names(phyloseq::sample_sums(ps.cf))]*100)
  
  #order vectors
  pfv <- pfv[names(p.if)]
  p.pf <- p.pf[names(p.if)]
  
  #cbind vectors into df
  sstab <- cbind(ov, ifv, p.if, pfv, p.pf)
  colnames(sstab) <- c("unfiltered.read.count", "filterI.read.count", "filterI.read.percent", "filterP.read.count", "filterP.read.percent")
  
  # Build return list
  l.return = list()
  if (return.all==FALSE){
    return(ps.cf)
  } else {
    l.return[['filtered.phyloseq']] <- ps.cf
    l.return[['ntaxa.in.control']] <- npos
    l.return[['control.taxa.sequences']] <- rownames(tax.tab.subset)
    l.return[['taxonomy.of.control.taxa']] <- ttsn
    l.return[['read.count.table']] <- sstab
    l.return[['prevalence.filter.threshold']] <- prevf
    
  }
  
  return(l.return)
  
}

write.dataset.biom <- function(ps, filepath, fileprefix){
  
  #save ASV sequences to vector and rename for fasta format
  f.onames <- phyloseq::taxa_names(ps)
  phyloseq::taxa_names(ps) <- paste("ASV", 1:length(phyloseq::taxa_names(ps)), sep = "")
  names(f.onames) <- paste0(">", phyloseq::taxa_names(ps))
  
  #generate biom file
  suppressWarnings(ps.b <- biomformat::make_biom(
    data = as.matrix(t(phyloseq::otu_table(ps))),
    sample_metadata = as.data.frame(phyloseq::sample_data(ps)),
    observation_metadata = as.data.frame(phyloseq::tax_table(ps)), 
    matrix_element_type = "int"
  )
  )
  
  fa <- print(paste0(filepath, fileprefix, "_ASVs.fasta"))
  bo <- print(paste0(filepath, fileprefix, "_ASV_table.biom"))
  
  #write output
  write.table(x = f.onames, file = fa, quote = FALSE, sep = "\n", col.names = FALSE)
  #biom export
  biomformat::write_biom(x = ps.b, biom_file = bo)
  
  #return phyloseq object with taxa renamed to ASV1, etc.
  return(ps)
}

write.dataset <- function(ps, filepath, fileprefix){
  
  #save ASV sequences to vector and rename for fasta format
  f.onames <- phyloseq::taxa_names(ps)
  phyloseq::taxa_names(ps) <- paste("ASV", 1:length(phyloseq::taxa_names(ps)), sep = "")
  names(f.onames) <- paste0(">", phyloseq::taxa_names(ps))
  
  #generate otu table formatted for biom generation
  suppressWarnings(otu.tab <- as.matrix(t(phyloseq::otu_table(ps))))
  cb <- as.matrix(cbind(rownames(otu.tab), otu.tab))
  rcb <- as.matrix(rbind(colnames(cb), cb))
  rcb[1,1] <- "#ASVID"
  rownames(rcb) <- NULL
  colnames(rcb) <- NULL
  
  #generate tax table formatted for biom generation
  tax.tab <- as.data.frame(phyloseq::tax_table(ps))
  tax.tab$taxonomy <- paste(tax.tab$Kingdom, tax.tab$Phylum, tax.tab$Class, 
                            tax.tab$Order, tax.tab$Family, tax.tab$Genus, tax.tab$Species, sep = ";")
  cbt <- as.matrix(cbind(rownames(tax.tab), tax.tab$taxonomy))
  rcbt <- as.matrix(rbind(c("#ASVID", "taxonomy"), cbt))
  rownames(cbt) <- NULL
  colnames(cbt) <- NULL
  
  #generate sampledf table formatted for biom generation
  samdf <- suppressWarnings(as.matrix(phyloseq::sample_data(ps)))
  cbs <- as.matrix(cbind(rownames(samdf), samdf))
  rcbs <- as.matrix(rbind(colnames(cbs), cbs))
  rcbs[1,1] <- "#SampleID"
  rownames(rcbs) <- NULL
  colnames(rcbs) <- NULL
  
  #create output string
  fa <- print(paste0(filepath, fileprefix, "_ASVs.fasta"))
  otb <- print(paste0(filepath, fileprefix, "_ASV_table.txt"))
  ttb <- print(paste0(filepath, fileprefix, "_ASV_taxonomy.txt"))
  stb <- print(paste0(filepath, fileprefix, "_sample_data.txt"))
  
  
  #write output
  #ASV fasta 
  write.table(x = f.onames, file = fa, quote = FALSE, sep = "\n", col.names = FALSE)
  #otu.tab
  write.table(x = rcb, file = otb, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  #tax.tab
  write.table(x = rcbt, file = ttb, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  #sampledf
  write.table(x = rcbs, file = stb, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  
  #return phyloseq object with taxa renamed to ASV1, etc.
  return(ps)
}
