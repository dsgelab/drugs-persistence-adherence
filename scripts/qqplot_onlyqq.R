packs <- c("qqman","optparse","data.table","R.utils")


for (p in packs) {
  if( !require(p, character.only = T)) {
    print(p)
    install.packages( p,  repos = c(CRAN = "http://cran.r-project.org") )
  }
}


option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-c","--chrcol"), type="character", default="CHR",
              help="chromosome column [default= %default]", metavar="character"),
  make_option(c("-p","--pval_col"), type="character", default="p.value",
              help="pvalue column [default= %default]. This can be a comma separated list and plots will be generated for each of these", metavar="character"),
  make_option(c("-b","--bp_col"), type="character", default="POS",
              help="bp column [default= %default]", metavar="character"),
  make_option(c("-l","--loglog_pval"), type="integer", default=10,
              help="-log10 p-val threshold for using log-log scale in manhattan plot [default= %default]", metavar="integer"),
  make_option(c("-y","--loglog_ylim"), type="integer", default=324,
              help="-log10 p-val limit for y-axis of log-log manhattan [default= %default]", metavar="integer"),
  make_option(c("--mlog10p"), type="logical", default=FALSE,
              help="whether the p-values are -log10 or not [default= %default]", metavar="logical"),
  make_option(c("-m","--minrep_col"), type="character",
              help="if given then chr:bp:ref:alt identifier assumed and chr and bp are read from there [default= %default]", metavar="character"),
  make_option(c("--info_col"), type="character", default="imputationInfo",
              help="info column [default= %default]", metavar="character"),
  make_option(c("--af_col"), type="character", default="AF_Allele2",
              help="AF column [default= %default]", metavar="character"),
  make_option(c("--info_threshold"), type="numeric", default=0.6,
              help="info threshold", metavar="numeric"),
  make_option(c("--af_threshold"), type="numeric", default=0.01,
              help="AF threshold", metavar="numeric")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, positional_arguments=0);

file <- opt$options$file
print(paste("reading file:", file))

data <- fread(file, header=T)

options(bitmapType='cairo')

print(str(opt))
bp_col <- opt$options$bp_col
chr_col <- opt$options$chrcol
info_col <- opt$options$info_col
af_col <- opt$options$af_col
info_threshold <- opt$options$info_threshold
af_threshold <- opt$options$af_threshold


print(summary(data))
print( summary( data[[chr_col]] ) )
#colnames(data) <- toupper( colnames(data) )

pcols <- unlist(strsplit(opt$options$pval_col,","))

output_prefix=file

if( !is.null(opt$options$out)) {
  output_prefix=opt$options$out
}


if(! is.null(opt$options$minrep_col ) ) {
  print("getting BP and CHR from minrepid")
  split <- strsplit(as.character(data[[opt$options$minrep_col]]), ":")
  data[[bp_col]] <- unlist( lapply( split, function(x) as.numeric(x[2]) ))
  data[[chr_col]] <- unlist( lapply( split, function(x) x[1] ))
}

print(append(pcols,c(bp_col,chr_col)))

if( any( ! append(pcols,c(bp_col,chr_col)) %in% colnames(data)   )) {
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


for( pcol in pcols) {
  subdata <- data[ !is.na(data[[pcol]]) & is.numeric( data[[pcol]]  ) ]
  subdata <- subdata[ subdata[[info_col]]>info_threshold & subdata[[af_col]]>af_threshold ]
  if (opt$options$mlog10p) {
    data[[pcol]] <- ifelse(10^-data[[pcol]] < 5e-324, 5e-324, 10^-data[[pcol]])
  }
  lambda  <- round(  quantile(  (qchisq(1-subdata[[pcol]], 1) ), probs=quants ) / qchisq(quants,1), 3)
  png( paste(output_prefix,"_", pcol ,"_qqplot.png", sep="" ))
  qq(subdata[[pcol]], col = "dodgerblue3")
  #title(c(file, paste("\n", "\nlambda ", quants, ": ", lambda, sep="" )), line=-0.2)
  dev.off()
  

  
}