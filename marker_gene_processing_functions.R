##functions for processing compositional microbiome data
#all scripts written by BPB, 070418


estimate.WSthreshold <- function(ps, WSmin=1e-4, WSmax=2e-4, WSstep=1e-5, controlID) {
  
  #Build param lists
  l.t <- seq(from = WSmin, to = WSmax, by = WSstep)
  nt <- length(l.t)
  tvec <- c()
  svec <- c()
  pvec <- c()
  
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
      pvec[i] <- sum(phyloseq::sample_sums(ps.if))/sum(phyloseq::sample_sums(ps))*100
      
    },
    error=function(e){cat("Warning :",conditionMessage(e), "\n")})
  }
  names(tvec) <- c(l.t)
  names(svec) <- c(l.t)
  names(pvec) <- c(l.t)
  
  df <- as.data.frame(cbind(tvec, svec, pvec, as.numeric(paste0(names(tvec)))))
  colnames(df) <- c("control.taxa.count", "read.count", "read.percent", "threshold.value")
  rownames(df) <- seq(1:length(l.t))
  df
  
}

estimate.ASthreshold <- function(ps, WSF, Pmin=NULL, Pmax=NULL, Pstep=NULL, CVmin=NULL, 
                                 CVmax=NULL, CVstep=NULL, RAFmin=NULL, RAFmax=NULL, RAFstep=NULL, independent=FALSE){
  
  if(is.null(WSF)){
    message('Warning: Estimation of AS filters will not be accurate without first applying WS filter')
  }
  
  #perform WS filtering
  filterfx = function(x){
    x[(x / sum(x)) < WSF] <- 0
    return(x)
  }
  
  ps.ws <- phyloseq::transform_sample_counts(ps, fun = filterfx)
  
  if(isTRUE(independent)){
    ps.wsi <- ps.ws
  }
  
  #RELATIVE ABUNDANCE
  #build param lists
  if(any(is.null(RAFmin), is.null(RAFmax), is.null(RAFstep))){
    dfr <- NULL
  } else {
    l.r <- seq(from = RAFmin, to = RAFmax, by = RAFstep)
    nr <- length(l.r)
    rvec <- c()
    
    for (i in 1:nr){
      tryCatch({
        #loop through values and filter
        raf <- sum(phyloseq::taxa_sums(ps.ws)) * l.r[i]
        ps.ws <- phyloseq::prune_taxa(taxa_sums(ps.ws)>=raf, ps.ws)
        rvec[i] <- phyloseq::ntaxa(ps.ws)
        
      },
      error=function(e){cat("Warning :r",conditionMessage(e), "\n")})
    }
    #create df
    dfr <- as.data.frame(cbind(l.r, rvec))
    colnames(dfr) <- c("relative.abundance.filter", "ASV.count")
    rownames(dfr) <- seq(1:length(l.r))
  }
  
  #revert, if desired
  if(isTRUE(independent)){
    ps.ws <- ps.wsi
  }
  
  #CV
  #build param lists
  if(any(is.null(CVmin), is.null(CVmax), is.null(CVstep))){
    dfc <- NULL
  } else {
    l.c <- seq(from = CVmin, to = CVmax, by = CVstep)
    nc <- length(l.c)
    cvec <- c()
    
    for (i in 1:nc){
      tryCatch({
        #loop through values and filter
        ps.ws <- phyloseq::filter_taxa(ps.ws, function(x) sd(x)/mean(x) > l.c[i], TRUE)
        cvec[i] <- phyloseq::ntaxa(ps.ws)
        
      },
      error=function(e){cat("Warning :c",conditionMessage(e), "\n")})
    }
    #create df
    dfc <- as.data.frame(cbind(l.c, cvec))
    colnames(dfc) <- c("CV.filter", "ASV.count")
    rownames(dfc) <- seq(1:length(l.c))
  }
  
  #revert, if desired
  if(isTRUE(independent)){
    ps.ws <- ps.wsi
  }
    
  #PREVALENCE
  #Build param lists
  if(any(is.null(Pmin), is.null(Pmax), is.null(Pstep))){
    dfp <- NULL
  } else {
  l.p <- seq(from = Pmin, to = Pmax, by = Pstep)
  np <- length(l.p)
  pvec <- c()
  taxa.cvec <- c()
  
  #ASV table error protection
  if(as.logical(class(phyloseq::otu_table(ps.ws))[1] == "otu_table") && 
     as.logical(taxa_are_rows(phyloseq::otu_table(ps.ws)) == TRUE)){
    otu.tab <- as.matrix(phyloseq::otu_table(ps.ws))
  } else {
    otu.tab <- as.matrix(t(phyloseq::otu_table(ps.ws)))
  }
  
  #get prevalence list for each taxon
  taxa.plist <- apply(X = otu.tab, MARGIN = 1, FUN = function(x){names(x)[which(x!=0)]})
  #populate vector of sample counts per ASV
  for(j in 1:length(taxa.plist)){
    taxa.cvec[j] <- length(taxa.plist[[j]])
  }
  
  for (i in 1:np){
    tryCatch({
      #apply ASV names and filter ASVs below PF
      names(taxa.cvec) <- names(taxa.plist)
      prev.count <- phyloseq::nsamples(ps.ws)* l.p[i]
      taxa.cvec.f <- taxa.cvec[which(taxa.cvec > prev.count)]
      tn.cvec.f <- names(taxa.cvec.f)
      #filter ps
      ps.ws <- phyloseq::prune_taxa(tn.cvec.f, ps.ws)
      pvec[i] <- phyloseq::ntaxa(ps.ws)
      
    },
    error=function(e){cat("Warning :p",conditionMessage(e), "\n")})
  }
  
  #create df
  dfp <- as.data.frame(cbind(l.p, pvec))
  colnames(dfp) <- c("prevalence.filter", "ASV.count")
  rownames(dfp) <- seq(1:length(l.p))
  }
  
  #create ASV df
  #reload ps.ws
  ps.ws <- phyloseq::transform_sample_counts(ps, fun = filterfx)
  #build df vectors
  ts <- taxa_sums(ps.ws)
  tsp <- taxa_sums(ps.ws)/sum(taxa_sums(ps.ws)) * 100
  namevec <- names(ts)
  #set tvec to null if no prev stats desired
  if(any(is.null(Pmin), is.null(Pmax), is.null(Pstep))){
    prev <- NULL
    prevp <- NULL
  } else {
    #revert to unfiltered taxa.cvec and otu.tab
    taxa.cvec <- c()
    #ASV table error protection
    if(as.logical(class(phyloseq::otu_table(ps.ws))[1] == "otu_table") && 
       as.logical(taxa_are_rows(phyloseq::otu_table(ps.ws)) == TRUE)){
      otu.tab <- as.matrix(phyloseq::otu_table(ps.ws))
    } else {
      otu.tab <- as.matrix(t(phyloseq::otu_table(ps.ws)))
    }
    #get prevalence list for each taxon
    taxa.plist <- apply(X = otu.tab, MARGIN = 1, FUN = function(x){names(x)[which(x!=0)]})
    #populate vector of sample counts per ASV
    for(j in 1:length(taxa.plist)){
      taxa.cvec[j] <- length(taxa.plist[[j]])
    }
    names(taxa.cvec) <- names(taxa.plist)
    prev <- as.numeric(paste0(taxa.cvec[namevec])) 
    prevp <- as.numeric(paste0(taxa.cvec[namevec]))/phyloseq::nsamples(ps.ws) * 100
  }
  
  #back to building df vectors
  cv.asv <- apply(otu.tab[namevec,], MARGIN = 1, FUN = function(x) sd(x)/mean(x))
  tax.tab <- phyloseq::tax_table(ps.ws)[namevec,]
  #create df
  df.asv <- cbind.data.frame(ts, tsp, prev, prevp, cv.asv, tax.tab)
  colnames(df.asv)[1:5] <- c("ASV.read.count", "ASV.read.percent", "ASV.prevalence", "ASV.prevalence.percent", "ASV.CV")
  rownames(df.asv) <- seq(1:nrow(df.asv))
  
  # Build return list
  l.return = list()
    l.return[['relative.abundance.filtering.stats']] <- dfr
    l.return[['CV.filtering.stats']] <- dfc
    l.return[['prevalence.filtering.stats']] <- dfp
    l.return[['ASV.filtering.stats']] <- df.asv
    
  return(l.return)
}


filter.dataset <- function(ps, controlID=NULL, controlCAT=NULL, controlFACTOR=NULL, minLIB=NULL, WSF=NULL, RAF=NULL, CVF=NULL, PF=NULL, return.all=TRUE){
  
  #remove samples < minlib
  if(is.null(minLIB)){
    ps = ps
  } else {
    ps = phyloseq::prune_samples(phyloseq::sample_sums(ps)>=minLIB, ps)
  }
  
  #create unfiltered sample sum vector
  ov <- phyloseq::sample_sums(ps)
  
  #WS filtering
  if(is.null(WSF)){
    ps.ws <- ps
  } else {
  filterfx = function(x){
    x[(x / sum(x)) < WSF] <- 0
    return(x)
  }
  
  #create WS filtered object
  ps.ws <- phyloseq::transform_sample_counts(ps, fun = filterfx)
  }
  sampledf <- suppressWarnings(as.matrix(phyloseq::sample_data(ps.ws)))
  #create per sample filtered sample sum vector
  ifv <- phyloseq::sample_sums(ps.ws)
  #calculate percent filtered, individual
  p.if <- phyloseq::sample_sums(ps.ws)/phyloseq::sample_sums(ps)*100
  
  #ASV table error protection
  if(as.logical(class(phyloseq::otu_table(ps.ws))[1] == "otu_table") && 
     as.logical(taxa_are_rows(phyloseq::otu_table(ps.ws)) == TRUE)){
    otu.tab <- as.matrix(phyloseq::otu_table(ps.ws))
  } else {
    otu.tab <- as.matrix(t(phyloseq::otu_table(ps.ws)))
  }
  
  #control sample filtering
  if(any(c(is.null(controlID), is.null(controlCAT), is.null(controlFACTOR)))){
    npos <- NULL
    tax.tab.subset <- NULL
    ttsn <- NULL
  } else {
  #control taxa count
  npos <- nrow(otu.tab[which(otu.tab[,match(controlID, colnames(otu.tab))] != 0),])
  
  #taxonomy of taxa in positive control
  tax.tab <- phyloseq::tax_table(ps.ws)
  taxanames.control <- rownames(otu.tab[which(otu.tab[,match(controlID, colnames(otu.tab))] != 0),])
  tax.tab.subset <- tax.tab[taxanames.control] #taxonomy of taxa in positive control
  ttsn <- tax.tab.subset
  rownames(ttsn) <- NULL
  
  #remove controls
  filtered.names <- rownames(sampledf[which(sampledf[,match(controlCAT, colnames(sampledf))] != controlFACTOR),])
  sampledf.s <- as.data.frame(sampledf[filtered.names,])
  phyloseq::sample_data(ps.ws) <- phyloseq::sample_data(sampledf.s)
  }
  
  #AS filtering
  #relative abundance filter
  if(is.null(RAF)){
    ps.ws <- ps.ws
  } else {
    raf <- sum(phyloseq::taxa_sums(ps.ws)) * RAF
    ps.ws <- phyloseq::prune_taxa(taxa_sums(ps.ws)>=raf, ps.ws)
  }
  
  #CV filter
  if(is.null(CVF)){
    ps.ws <- ps.ws
  } else {
    ps.ws <- phyloseq::filter_taxa(ps.ws, function(x) sd(x)/mean(x) > CVF, TRUE)
  }
  
  #prevalence filter
  if(is.null(PF)){
    ps.ws <- ps.ws
  } else {
    #more ASV table error protection
    if(as.logical(class(phyloseq::otu_table(ps.ws))[1] == "otu_table") && 
       as.logical(taxa_are_rows(phyloseq::otu_table(ps.ws)) == TRUE)){
      otu.tab <- as.matrix(phyloseq::otu_table(ps.ws))
    } else {
      otu.tab <- as.matrix(t(phyloseq::otu_table(ps.ws)))
    }
    #get prevalence list for each taxon
    taxa.plist <- apply(X = otu.tab, MARGIN = 1, FUN = function(x){names(x)[which(x!=0)]})
    #build vector of sample counts per ASV
    taxa.cvec <- c()
    for(i in 1:length(taxa.plist)){
      taxa.cvec[i] <- length(taxa.plist[[i]])
    }
    #apply ASV names and filter ASVs below PF
    names(taxa.cvec) <- names(taxa.plist)
    prev.count <- phyloseq::nsamples(ps.ws)*PF
    taxa.cvec.f <- taxa.cvec[which(taxa.cvec > prev.count)]
    tn.cvec.f <- names(taxa.cvec.f)
    #filter ps
    ps.ws <- phyloseq::prune_taxa(tn.cvec.f, ps.ws)
    
  }
  
  #create AS filter sample sum vector
  pfv <- phyloseq::sample_sums(ps.ws)
  #calculate percent filtered, prevalence
  p.pf <- suppressWarnings(phyloseq::sample_sums(ps.ws)/phyloseq::sample_sums(ps)[names(phyloseq::sample_sums(ps.ws))]*100)
  
  #order vectors
  pfv <- pfv[names(p.if)]
  p.pf <- p.pf[names(p.if)]
  
  #cbind vectors into df
  sstab <- cbind(ov, ifv, p.if, pfv, p.pf)
  colnames(sstab) <- c("unfiltered.read.count", "WSfiltered.read.count", "WSfiltered.read.percent", "ASfiltered.read.count", "ASfiltered.read.percent")
  
  # Build return list
  l.return = list()
  if (return.all==FALSE){
    return(ps.ws)
  } else {
    l.return[['filtered.phyloseq']] <- ps.ws
    l.return[['ntaxa.in.control']] <- npos
    l.return[['control.taxa.sequences']] <- rownames(tax.tab.subset)
    l.return[['taxonomy.of.control.taxa']] <- ttsn
    l.return[['read.count.table']] <- sstab
    l.return[['relative.abundance.filter.read.count']] <- raf
    l.return[['prevalence.filter.sample.count']] <- prev.count
    
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
