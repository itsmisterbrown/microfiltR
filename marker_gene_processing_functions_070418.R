##functions for processing compositional microbiome data
#all scripts written by BPB, 070418


filter.threshold <- function(ps, Tmin=1e-4, Tmax=2e-4, Tstep=1e-5, controlID) {
  
  #Build param lists
  l.t <- seq(from = Tmin, to = Tmax, by = Tstep)
  nt <- length(l.t)
  tvec <- c()
  
  for (i in 1:nt){
      tryCatch({
        #loop through values
        filterfx = function(x){
          x[(x / sum(x)) < (l.t[i])] <- 0
          return(x)
        }
        
        ps.if <- phyloseq::transform_sample_counts(ps, fun = filterfx)
        
        #FILTERING
        otu.tab <- as.matrix(t(phyloseq::otu_table(ps.if)))
        tvec[i] <- nrow(otu.tab[which(otu.tab[,match(controlID, colnames(otu.tab))] != 0),])

      },
      error=function(e){cat("Warning :",conditionMessage(e), "\n")})
  }
  names(tvec) <- c(l.t)
  return(tvec)
  
}


filter.object <- function(ps, controlID=NULL, controlCAT=NULL, controlFACTOR=NULL, TF=1e-4, PF=1e-5, return.all=TRUE){
  
  #per-sample filter function
  filterfx = function(x){
    x[(x / sum(x)) < TF] <- 0
    return(x)
  }
  
  #create filtered object
  ps.if <- phyloseq::transform_sample_counts(ps, fun = filterfx)
  sampledf <- suppressWarnings(as.matrix(phyloseq::sample_data(ps.if)))
  
  #FILTERING
  otu.tab <- as.matrix(t(phyloseq::otu_table(ps.if)))
  npos <- nrow(otu.tab[which(otu.tab[,match(controlID, colnames(otu.tab))] != 0),])
  #taxonomy of taxa in positive control
  tax.tab <- phyloseq::tax_table(ps.if)
  taxanames.control <- rownames(otu.tab[which(otu.tab[,match(controlID, colnames(otu.tab))] != 0),])
  tax.tab.subset <- tax.tab[taxanames.control] #taxonomy of taxa in positive control
  ttsn <- tax.tab.subset
  rownames(ttsn) <- NULL
  
  
  #calculate percentage of taxa remaining after per-samplefilter
  ds.filtered <- (sum(phyloseq::sample_sums(ps.if))/sum(phyloseq::sample_sums(ps))*100) #dataset
  pc.filtered <- (phyloseq::sample_sums(ps.if)[controlID]/phyloseq::sample_sums(ps)[controlID]*100) #pc
  
  #remove controls
  filtered.names <- rownames(sampledf[which(sampledf[,match(controlCAT, colnames(sampledf))] != controlFACTOR),])
  sampledf.s <- as.data.frame(sampledf[filtered.names,])
  phyloseq::sample_data(ps.if) <- phyloseq::sample_data(sampledf.s)
  
  #dataset-wide prevalence filter %
  #determine taxa sums for cutoff %
  prevf <- sum(phyloseq::taxa_sums(ps.if)) * PF
  ps.cf <- phyloseq::prune_taxa(taxa_sums(ps.if)>=prevf, ps.if)
  
  #calculate percentage of taxa remaining after prevalence filter
  ds.filteredp <- (sum(phyloseq::sample_sums(ps.cf))/sum(phyloseq::sample_sums(ps))*100) #dataset, prevalence
  filter.vec <- c(pc.filtered, ds.filtered, ds.filteredp)
  names(filter.vec) <- c("pc reads filterI", "dataset reads filterI", "dataset reads filterP")
  
  # Build return list
  l.return = list()
  if (return.all==FALSE){
    return(ps.cf)
    } else {
    l.return[['filtered.phyloseq']] <- ps.cf
    l.return[['ntaxa.in.control']] <- npos
    l.return[['control.taxa.sequences']] <- rownames(tax.tab.subset)
    l.return[['taxonomy.of.control.taxa']] <- ttsn
    l.return[['read.filtering.stats']] <- filter.vec
    l.return[['prevalence.filter.threshold']] <- prevf
    
  }
  
  return(l.return)
  
}

write.phyloseq.biom <- function(ps, filepath, fileprefix){
  
  #save ASV sequences to vector and rename for fasta format
  f.onames <- phyloseq::taxa_names(ps)
  phyloseq::taxa_names(ps) <- paste(">ASV", 1:length(phyloseq::taxa_names(ps)), sep = "")
  names(f.onames) <- phyloseq::taxa_names(ps)
  phyloseq::taxa_names(ps) <- paste("ASV", 1:length(phyloseq::taxa_names(ps)), sep = "")
  
  #generate biom file
  suppressWarnings(ps.b <- biomformat::make_biom(
    data = as.matrix(t(phyloseq::otu_table(ps))),
    sample_metadata = as.data.frame(phyloseq::sample_data(ps)),
    observation_metadata = as.data.frame(phyloseq::tax_table(ps)), 
    matrix_element_type = "int"
  )
  )
  
  fa <- print(paste0(filepath, fileprefix, ".fasta"))
  bo <- print(paste0(filepath, fileprefix, ".biom"))
  
  #write output
  write.table(x = f.onames, file = fa, quote = FALSE, sep = "\n", col.names = FALSE)
  #biom export
  biomformat::write_biom(x = ps.b, biom_file = bo)
  
  return(ps)
}

write.phyloseq <- function(ps, filepath, fileprefix){
  
  #save ASV sequences to vector and rename for fasta format
  f.onames <- phyloseq::taxa_names(ps)
  phyloseq::taxa_names(ps) <- paste(">ASV", 1:length(phyloseq::taxa_names(ps)), sep = "")
  names(f.onames) <- phyloseq::taxa_names(ps)
  phyloseq::taxa_names(ps) <- paste("ASV", 1:length(phyloseq::taxa_names(ps)), sep = "")
  
  #generate otu table 
  otu.tab <- as.matrix(t(phyloseq::otu_table(ps)))
  cb <- as.matrix(cbind(rownames(otu.tab), otu.tab))
  rcb <- as.matrix(rbind(colnames(cb), cb))
  rcb[1,1] <- "#OTU ID"
  rownames(rcb) <- NULL
  colnames(rcb) <- NULL
  
  fa <- print(paste0(filepath, fileprefix, ".fasta"))
  otb <- print(paste0(filepath, fileprefix, ".txt"))
  
  #write output
  write.table(x = f.onames, file = fa, quote = FALSE, sep = "\n", col.names = FALSE)
  #biom export
  write.table(x = rcb, file = otb, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  
  return(ps)
}
