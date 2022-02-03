packs <- c("qqman","optparse","data.table","R.utils")


for (p in packs) {
  if( !require(p, character.only = T)) {
    print(p)
    install.packages( p,  repos = c(CRAN = "http://cran.r-project.org") )
  }
}

rm(list = ls())

setwd('/home/ivm/drugs/results/R8_20220112/sumstats_stop/')

stats <- system('ls *rsid.gz', intern = T)

for (s in stats){
  
  file <- s
  print(paste("reading file:", file))
  
  data <- fread(file, header=T)
  
  bp_col <- "pos"
  chr_col <- "#chrom"
  info_col <- "info"
  af_col <- "af_alt"
  pcol <- "pval"
  info_threshold <- 0.6
  af_threshold <- 0.01
  
  
  print(summary(data))
  print( summary( data[[chr_col]] ) )
  
  output_prefix = paste0('/home/ivm/drugs/results/R8_20220112/plots_stop/',file)

  if( any( ! c(bp_col,chr_col, pcol) %in% colnames(data)   )) {
    stop( paste0("All required columns do not exist in the data: ", paste(pcols,sep=",", collapse=""),",", bp_col, ",",chr_col,  collapse="" ))
  }
  
  print(summary(as.factor(data[[chr_col]])))
  
  data[[chr_col]] <- gsub("chr","",data[[chr_col]])
  data[[chr_col]] <- gsub("X|chrX","23",data[[chr_col]])
  data[[chr_col]] <- gsub("Y|chrY","24",data[[chr_col]])
  data[[chr_col]] <- gsub("MT|chrMT|M|chrM","25",data[[chr_col]])
  
  data[[chr_col]] <- as.numeric(data[[chr_col]])
  data <- data[ !is.na(data[[chr_col]]) ]
  
  quants <- c(0.7,0.5,0.1,0.01, 0.001)
  
  subdata <- data[ !is.na(data[[pcol]]) & is.numeric( data[[pcol]]  ) ]
  
  # print("filtering INFO and AF")
  # subdata <- subdata[ subdata[[info_col]]>info_threshold & subdata[[af_col]]>=af_threshold  & subdata[[af_col]]<=1-af_threshold]
  # 
  # fwrite(subdata, paste0(output_prefix, ".txt.gz"), sep = "\t", quote = F, compress = "gzip")

  lambda  <- round(  quantile(  (qchisq(1-subdata[[pcol]], 1) ), probs=quants ) / qchisq(quants,1), 3)
  
  # png( paste(output_prefix,"_", pcol ,"_qqplot.png", sep="" ))
  # qq(subdata[[pcol]])
  # title(c(file, paste("\n", "\nlambda ", quants, ": ", lambda, sep="" )), line=-0.2)
  # dev.off()
  
  print("subsetting p-vals < 0.01 for manhattan...")
  subdata <- subdata[ subdata[[pcol]]<0.01 & subdata[[pcol]]>0 ]

  print( paste0("Plotting manhattan with ", nrow(subdata), " variants") )
  print( summary(subdata[[pcol]] ))
  
  png( paste(output_prefix,"_",pcol,"_manhattan_title.png", sep=""), width=1000, height=400)
  logs <- -log10(subdata[[pcol]])
  manhattan( data.table(subdata[,c(bp_col,pcol,chr_col), with=F]) , chr=chr_col, bp=bp_col, p=pcol, ylim=c( 2,max(logs)+1), main=sub(".af_0.01.info_0.6.gz.rsid.gz","",file))
  dev.off()
  
}
