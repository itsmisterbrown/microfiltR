##functions for processing compositional microbiome data
#all scripts written by BPB, 071318

#ERROR PROTECTION FUNCTIONS
format.ASV.tab <- function(ps){
  if(as.logical(class(phyloseq::otu_table(ps))[1] == "otu_table") && 
     as.logical(taxa_are_rows(phyloseq::otu_table(ps)) == TRUE)){
    asv.tab <- as.matrix(phyloseq::otu_table(ps))
  } else {
    asv.tab <- as.matrix(t(phyloseq::otu_table(ps)))
  }
  asv.tab
}


#PROCESSING FUNCTIONS
#standardization 
standardize.median <- function(ps){
  median.rc <- median(phyloseq::sample_sums(ps))
  ps.t <- phyloseq::transform_sample_counts(ps, fun = function(x) round(median.rc * (x/sum(x))))
  ps.t
}

getCV <- function(ps, WSF=NULL, CVrange, CVstep){
  
  if(is.null(WSF)){
    message('Warning: Estimation of AS filters will not be accurate without first applying WS filter')
  }
  
  #perform WS filtering
  filterfx = function(x){
    x[(x / sum(x)) < WSF] <- 0
    return(x)
  }
  
  ps.ws <- phyloseq::transform_sample_counts(ps, fun = filterfx)
  #standardize to median sample depth
  ps.wsm <- standardize.median(ps.ws)
  
  #build param vectors
  if(any(is.null(CVrange), is.null(CVstep))){
    dfc <- NULL
  } else {
    l.c <- seq(from = CVrange[1], to = CVrange[2], by = CVstep)
    nc <- length(l.c)
    cvec <- c()
    
    for (i in 1:nc){
      tryCatch({
        #loop through values and filter
        ps.wsm <- phyloseq::filter_taxa(ps.wsm, function(x) sd(x)/mean(x) > l.c[i], TRUE)
        cvec[i] <- phyloseq::ntaxa(ps.wsm)
        
      },
      error=function(e){cat("Warning :c",conditionMessage(e), "\n")})
    }
    #create df
    dfc <- as.data.frame(cbind(l.c, cvec))
    colnames(dfc) <- c("CV.filter", "ASV.count")
    rownames(dfc) <- seq(1:length(l.c))
  }
  
  dfc
}

getRA <- function(ps, WSF=NULL, RAFrange, RAFstep){
  
  if(is.null(WSF)){
    message('Warning: Estimation of AS filters will not be accurate without first applying WS filter')
  }
  
  #perform WS filtering
  filterfx = function(x){
    x[(x / sum(x)) < WSF] <- 0
    return(x)
  }
  
  ps.ws <- phyloseq::transform_sample_counts(ps, fun = filterfx)
  
  #build param vectors
  if(any(is.null(RAFrange), is.null(RAFstep))){
    dfr <- NULL
  } else {
    l.r <- seq(from = RAFrange[1], to = RAFrange[2], by = RAFstep)
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
  
  dfr
}

getPrev <- function(ps, WSF=NULL, Prange, Pstep){
  
  if(is.null(WSF)){
    message('Warning: Estimation of AS filters will not be accurate without first applying WS filter')
  }
  
  #perform WS filtering
  filterfx = function(x){
    x[(x / sum(x)) < WSF] <- 0
    return(x)
  }
  
  ps.ws <- phyloseq::transform_sample_counts(ps, fun = filterfx)
  
  asv.tab <- format.ASV.tab(ps.ws)
  
  #build param vectors
  l.p <- seq(from = Prange[1], to = Prange[2], by = Pstep)
  np <- length(l.p)
  pvec <- c()
  taxa.cvec <- c()
  
  #get prevalence list for each taxon
  taxa.plist <- apply(X = asv.tab, MARGIN = 1, FUN = function(x){names(x)[which(x!=0)]})
  #populate vector of sample counts per ASV
  for(j in 1:length(taxa.plist)){
    taxa.cvec[j] <- length(taxa.plist[[j]])
  }
  
  #loop through params
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
  dfp <- cbind.data.frame(l.p, pvec)
  colnames(dfp) <- c("prevalence.filter", "ASV.count")
  rownames(dfp) <- seq(1:length(l.p))
  #name taxa prevalence vector
  names(taxa.cvec) <- names(taxa.plist)
  
  # Build return list
  l.return = list()
  l.return[['prevalence.filtering.stats']] <- dfp
  l.return[['ASV.prevalence.count']] <- taxa.cvec
  
  return(l.return)
  
}
CVfilter <- function(ps, WSF=NULL, CVF){
  
  if(is.null(WSF)){
    message('Warning: WS filtering is highly recommended to reduce cross contamination')
  }
  
  #perform WS filtering
  filterfx = function(x){
    x[(x / sum(x)) < WSF] <- 0
    return(x)
  }
  
  ps.ws <- phyloseq::transform_sample_counts(ps, fun = filterfx)

  #standardize to median sample depth
  ps.wsm <- standardize.median(ps.ws)
  
  #perform filter
  ps.ws <- phyloseq::filter_taxa(ps.ws, function(x) sd(x)/mean(x) > CVF, TRUE)
  ps.ws
  
}



RAfilter<- function(ps, WSF=NULL, RAF){
  
  if(is.null(WSF)){
    message('Warning: WS filtering is highly recommended to reduce cross contamination')
  }
  
  #perform WS filtering
  filterfx = function(x){
    x[(x / sum(x)) < WSF] <- 0
    return(x)
  }
  
  ps.ws <- phyloseq::transform_sample_counts(ps, fun = filterfx)
  
  #perform filter
  raf <- sum(phyloseq::taxa_sums(ps.ws)) * RAF
  ps.ws <- phyloseq::prune_taxa(taxa_sums(ps.ws)>=raf, ps.ws)
  ps.ws
  
}

Pfilter <- function(ps, WSF=NULL, PF){
  
  if(is.null(WSF)){
    message('Warning: WS filtering is highly recommended to reduce cross contamination')
  }
  
  #perform WS filtering
  filterfx = function(x){
    x[(x / sum(x)) < WSF] <- 0
    return(x)
  }
  
  ps.ws <- phyloseq::transform_sample_counts(ps, fun = filterfx)
  
  #format asv table
  asv.tab <- format.ASV.tab(ps.ws)
  
  #create sample count vector
  taxa.cvec <- c()
  #get prevalence list for each taxon
  taxa.plist <- apply(X = asv.tab, MARGIN = 1, FUN = function(x){names(x)[which(x!=0)]})
  #populate vector of sample counts per ASV
  for(j in 1:length(taxa.plist)){
    taxa.cvec[j] <- length(taxa.plist[[j]])
  }
  names(taxa.cvec) <- names(taxa.plist)
  prev.count <- phyloseq::nsamples(ps.ws)* PF
  taxa.cvec.f <- taxa.cvec[which(taxa.cvec > prev.count)]
  tn.cvec.f <- names(taxa.cvec.f)
  #perform filter
  ps.ws <- phyloseq::prune_taxa(tn.cvec.f, ps.ws)
  ps.ws
  
}

#WRAPPER FUNCTIONS
estimate.WSthreshold <- function(ps, WSrange, WSstep, controlID) {
  
  #Build param lists
  l.t <- seq(from = WSrange[1], to = WSrange[2], by = WSstep)
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
      
      ps.ws <- phyloseq::transform_sample_counts(ps, fun = filterfx)
      asv.tab <- format.ASV.tab(ps.ws)
      
      #FILTERING
      
      tvec[i] <- nrow(asv.tab[which(asv.tab[,match(controlID, colnames(asv.tab))] != 0),])
      svec[i] <- sum(phyloseq::sample_sums(ps.ws))
      pvec[i] <- sum(phyloseq::sample_sums(ps.ws))/sum(phyloseq::sample_sums(ps))*100
      
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

estimate.ASthreshold <- function(ps, WSF, RAF=NULL, CVF=NULL, PF=NULL, controlID=NULL, controlCAT=NULL, controlFACTOR=NULL,
                                 minLIB=NULL, Prange=NULL, Pstep=NULL, CVrange=NULL, CVstep=NULL, RAFrange=NULL, RAFstep=NULL){
  
  if(is.null(WSF)){
    message('Warning: Estimation of AS filters will not be accurate without first applying WS filter')
  }
  
  #remove samples < minlib
  if(is.null(minLIB)){
    ps = ps
  } else {
    ps = phyloseq::prune_samples(phyloseq::sample_sums(ps)>=minLIB, ps)
  }
  
  #perform WS filtering
  filterfx = function(x){
    x[(x / sum(x)) < WSF] <- 0
    return(x)
  }
  
  ps.ws <- phyloseq::transform_sample_counts(ps, fun = filterfx)
  ps.wso <- ps.ws
  
  #control sample filtering
  if(any(c(is.null(controlID), is.null(controlCAT), is.null(controlFACTOR)))){
    message('If you plan to remove controls during filtering, then it is recommended to remove during estimation')
  } else {
    #create sampledf
    sampledf <- suppressWarnings(as.matrix(phyloseq::sample_data(ps.ws)))
    #remove controls
    filtered.names <- rownames(sampledf[which(sampledf[,match(controlCAT, colnames(sampledf))] != controlFACTOR),])
    sampledf.s <- as.data.frame(sampledf[filtered.names,])
    phyloseq::sample_data(ps.ws) <- phyloseq::sample_data(sampledf.s)
  }
  
  #RELATIVE ABUNDANCE
  #build param lists
  if(any(is.null(RAFrange), is.null(RAFstep))){
    gr <- NULL
  } else {
    gr <- getRA(ps = ps.ws, WSF = WSF, RAFrange = RAFrange, RAFstep = RAFstep)
  }
  
  #incorporate fixed threshold
  if(!is.null(RAF)){
    ps.ws <- RAfilter(ps = ps.ws, WSF = WSF, RAF = RAF)
  }
  
  #CV
  #build param lists
  if(any(is.null(CVrange), is.null(CVstep))){
    gc <- NULL
  } else {
    gc <- getCV(ps = ps.ws, WSF = WSF, CVrange = CVrange, CVstep = CVstep)
  }
  
  #incorporate fixed threshold
  if(!is.null(CVF)){
    ps.ws <- CVfilter(ps = ps.ws, WSF = WSF, CVF = CVF)
  }
  
  #PREVALENCE
  #Build param lists
  if(any(is.null(Prange), is.null(Pstep))){
    gp$prevalence.filtering.stats <- NULL
  } else {
    gp <- getPrev(ps = ps.ws, WSF = WSF, Prange = Prange, Pstep = Pstep)
  }
  
  #CREATE ASV DF
  #build df vectors
  ts <- taxa_sums(ps.wso)
  tsp <- taxa_sums(ps.wso)/sum(taxa_sums(ps.wso)) * 100
  namevec <- names(ts)
  #standardize to median sample depth for CV calculation
  ps.ws <- standardize.median(ps.wso)
  asv.tab <- format.ASV.tab(ps.ws)
  cv.asv <- apply(asv.tab[namevec,], MARGIN = 1, FUN = function(x) sd(x)/mean(x))
  tax.tab <- phyloseq::tax_table(ps.wso)[namevec,]
  
  #set prev vectors to null if no prev stats desired
  if(any(is.null(Prange), is.null(Pstep))){
    prev <- rep(NA, length(ts))
    prevp <- rep(NA, length(ts))
  } else {
    gp.reload <- getPrev(ps = ps.ws, WSF = WSF, Prange = Prange, Pstep = Pstep)
    taxa.cvec <- gp.reload$ASV.prevalence.count
    prev <- taxa.cvec[namevec]
    prevp <- prev/phyloseq::nsamples(ps.ws) * 100
  }
  
  #build df and rename
  df.asv <- cbind.data.frame(ts, tsp, prev, prevp, cv.asv, tax.tab)
  colnames(df.asv)[1:5] <- c("ASV.read.count", "ASV.read.percent", "ASV.prevalence", "ASV.prevalence.percent", "ASV.CV")
  rownames(df.asv) <- seq(1:nrow(df.asv))
  
  # Build return list
  l.return = list()
  l.return[['relative.abundance.filtering.stats']] <- gr
  l.return[['CV.filtering.stats']] <- gc
  l.return[['prevalence.filtering.stats']] <- gp$prevalence.filtering.stats
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
    message('Warning: WS filtering is highly recommended to reduce cross contamination')
    ps.ws <- ps
  } else {
    filterfx = function(x){
      x[(x / sum(x)) < WSF] <- 0
      return(x)
    }
    
    #create WS filtered object
    ps.ws <- phyloseq::transform_sample_counts(ps, fun = filterfx)
  }
  
  #create WS filtered sample sum vector
  ifv <- phyloseq::sample_sums(ps.ws)
  #calculate percent filtered, individual
  p.if <- phyloseq::sample_sums(ps.ws)/phyloseq::sample_sums(ps)*100
  
  asv.tab <- format.ASV.tab(ps.ws)
  
  #CONTROL (METADATA-BASED) SAMPLE REMOVAL
  if(any(c(is.null(controlID), is.null(controlCAT), is.null(controlFACTOR)))){
    npos <- NULL
    tax.tab.subset <- NULL
    ttsn <- NULL
  } else {
    #control taxa count
    npos <- nrow(asv.tab[which(asv.tab[,match(controlID, colnames(asv.tab))] != 0),])
    
    #taxonomy of taxa in positive control
    tax.tab <- phyloseq::tax_table(ps.ws)
    taxanames.control <- rownames(asv.tab[which(asv.tab[,match(controlID, colnames(asv.tab))] != 0),])
    tax.tab.subset <- tax.tab[taxanames.control] #taxonomy of taxa in positive control
    ttsn <- tax.tab.subset
    rownames(ttsn) <- NULL
    
    #remove controls
    sampledf <- suppressWarnings(as.matrix(phyloseq::sample_data(ps.ws)))
    filtered.names <- rownames(sampledf[which(sampledf[,match(controlCAT, colnames(sampledf))] != controlFACTOR),])
    sampledf.s <- as.data.frame(sampledf[filtered.names,])
    phyloseq::sample_data(ps.ws) <- phyloseq::sample_data(sampledf.s)
  }
  
  #AS filtering
  #relative abundance filter
  if(is.null(RAF)){
    ps.ws <- ps.ws
    raf <- NULL
  } else {
    ps.ws <- RAfilter(ps = ps.ws, WSF = WSF, RAF = RAF)
    raf <- RAF * sum(phyloseq::taxa_sums(ps.ws))
  }
  
  #CV filter
  if(is.null(CVF)){
    ps.ws <- ps.ws
  } else {
    ps.ws <- CVfilter(ps = ps.ws, WSF = WSF, CVF = CVF)
  }
  
  #prevalence filter
  if(is.null(PF)){
    ps.ws <- ps.ws
    prev.count <- NULL
  } else {
    ps.ws <- Pfilter(ps = ps.ws, WSF = WSF, PF = PF)
    prev.count <- phyloseq::nsamples(ps.ws) * PF
  }
  #create AS filter sample sum vector
  pfv <- phyloseq::sample_sums(ps.ws)
  #calculate percent filtered, AS
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
  
  #generate asv table formatted for biom generation
  asv.tab <- format.ASV.tab(ps)
  suppressWarnings(asv.tab <- as.matrix(asv.tab))
  cb <- as.matrix(cbind(rownames(asv.tab), asv.tab))
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
  #asv.tab
  write.table(x = rcb, file = otb, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  #tax.tab
  write.table(x = rcbt, file = ttb, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  #sampledf
  write.table(x = rcbs, file = stb, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  
  #return phyloseq object with taxa renamed to ASV1, etc.
  return(ps)
}


