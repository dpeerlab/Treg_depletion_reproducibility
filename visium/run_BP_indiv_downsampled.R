## run_BP_indiv_downsampled.R

# script to launch from command line for batch processing of BayesPrism 
# run the deconvolution for each section individually and downsample reads
# within each spot to 90% of total to assess the robustness of cell fraction 
# estimates to a partial loss of reads

# load libraries
require(BayesPrism)
require(Seurat)
require(dplyr)

# load processed data selected for only marker genes
load("/data/peer/sam/treg_dep/mouse/data/scrna.herman.allCT.deconvRiter1.pUp0.05.rdata")

# meta object that is loaded has cell sub-type and cell type labels 
cell.type.labels <- as.character(meta$cell.annot.labels)
cell.state.labels <- as.character(meta$cell.subtype.labels)
ref.dat.sig <- ref.dat.filtered

## change the gene names of ref.dat.sig to match visium
# convert to Seurat object
uniq_index = make.unique(rownames(ref.dat.sig), sep = "-")

rownames(ref.dat.sig) = uniq_index
rownames(meta) = uniq_index



# identify sample name and parse command line arguments

args = commandArgs(trailingOnly = T)

if(length(args) != 3){
  stop("enter the name of section and the out directory and the run number")
}

sample_names = args[1]
#sample_names = c('Ctrl1_A1')

outdir = args[2]

run_number = args[3]

outdir = file.path(outdir, paste0("run_", run_number))

# make out dir if it doesn't exit
if(!(dir.exists(outdir))){
  dir.create(outdir, recursive = T)
}

# load visium data
visium_dat_dir = '/data/peer/sam/treg_dep/visium/data/requant/results/'

visium_samples = list.dirs(visium_dat_dir, recursive = F) %>% basename()
visium_samples = visium_samples[grepl("IGO", visium_samples)]
visium_samples_names = gsub("_IGO_.+", "",  visium_samples)

# get visium samples to load
visium_samples_run = visium_samples[grepl(sample_names, visium_samples)] %>% sort()


visium_samples_files = lapply(visium_samples_run, function(x){
  file.path(visium_dat_dir, x, 'outs')
})


print(paste("Running for samples", visium_samples_run))

# load using seurat
d = lapply(visium_samples_files, Load10X_Spatial)
names(d) = sample_names

# set the section name
d = lapply(names(d), function(x){
  r = d[[x]]
  r@meta.data$section = x
  return(r)
})
names(d) = sample_names


# change d.m to just be the first section
# for compatibility with rest of script
d.m = d[[1]]

lib_size_cutoff = 1000
# calculate library size
library_size = apply(d.m@assays$Spatial@counts, 2, sum)
low_lib_size_spot = library_size < lib_size_cutoff

print(paste(sum(low_lib_size_spot), "spots removed due to lib size less than", lib_size_cutoff))
d.m = d.m[,!(low_lib_size_spot)]

# index of samples
d_idx = d.m$section

#### gene manipulation ####
# handling gene ID issues


# remove genes that could be an issue

## collapse transcripts with .1 and .2 values to individual entries for matching
# get out the gene expression 
X.filtered = d.m@assays$Spatial@counts %>% as.matrix() %>% t() 

# get duplicated genes
dup_gene_names = colnames(X.filtered)[grepl("\\.[0-9]+$", colnames(X.filtered))] %>% sort()

if(length(dup_gene_names) > 0){
  dup_gene_basename = gsub("\\.[0-9]$", "",dup_gene_names) %>% unique()
  
  # get out all of the genes that I don't need to sum
  X.nondup = X.filtered[,!(colnames(X.filtered) %in% dup_gene_names)]
  
  # loop through all the dup genes and sum and change to name without the .1 etc
  X.dup_sum = do.call(cbind, lapply(dup_gene_basename, function(dg){
    dup_gene_names_check = dup_gene_names[grepl(dg, dup_gene_names)]
    
    if(length(dup_gene_names_check) == 1) return(X.filtered[,dup_gene_names_check])
    else(apply(X.filtered[,colnames(X.filtered) %in% dup_gene_names_check], 1, sum))
    
  }))
  colnames(X.dup_sum) = dup_gene_basename
  
  X.filtered = cbind(X.nondup, X.dup_sum)
}

# sum the GFP and Foxp3 reads for Foxp3
X.filtered[,'Foxp3'] = X.filtered[,'Foxp3'] + X.filtered[,'sDTR-eGFP']
X.filtered = X.filtered[,colnames(X.filtered) != 'sDTR-eGFP']

# change to title case for now to help matching with the reference I'm using
colnames(X.filtered) = stringr::str_to_title(colnames(X.filtered))

#### downsampling ####
## downsampling to 90% of reads and then running deconvolution
# create new matrix that will be the downsampled one
# remove any genes that have 0 counts
gene_counts = apply(X.filtered, 2, sum)
gene_rm = gene_counts == 0
X.downsample = X.filtered[,!(gene_rm)]

# function to downsample
downsampled <- function(data,samplerate=0.8) {
  data.test <- apply(data,1,function(q) {
    names(q) <- colnames(data)
    samplepool <- character()
    for (i in names(q)) {
      samplepool <- append(samplepool,rep(i,times=q[i]))  
    }
    sampled <- sample(samplepool,size=samplerate*length(samplepool),replace = F)
    tab <- table(sampled)
    mat <- match(names(tab),names(q))
    toret=numeric(length <- length(q))
    names(toret) <- names(q)
    toret[mat] <- tab
    return(toret)
  })
  data.test = t(data.test)
  return(data.test)
}


# for each row in the matrix, downsample 90% of reads and then reassign
# I can cut down on the time on this by just trimming to only the genes that are being
# used to deconvolve
X.downsample = downsampled(X.downsample, samplerate = .9)

# remove unnessecary files and garbage collect
rm(X.filtered)
gc()

# create BP object
myPrism <- new.prism(
  reference=ref.dat.sig, 
  mixture=X.downsample,
  input.type="count.matrix", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = cell.state.labels,
  key=NULL,
  outlier.cut=0.01,
  outlier.fraction=0.1,
)

# run BP
bp.res <- run.prism(prism = myPrism, n.cores=40)

print(bp.res)

# save this file
saveRDS(bp.res, file = file.path(outdir, paste0(sample_names,  "_bpRes.rds")))

# save cell fractions
theta = get.fraction (bp=bp.res,
                      which.theta="final",
                      state.or.type="type")
write.csv(theta, file =  file.path(outdir, paste(sample_names,  "theta.csv", sep = "_")), 
          quote = F, row.names = T)

# save posterior CV
theta.cv <- bp.res@posterior.theta_f@theta.cv
write.csv(theta.cv, file =  file.path(outdir, paste(sample_names,  "thetaCV.csv", sep = "_")), 
          quote = F, row.names = T)


