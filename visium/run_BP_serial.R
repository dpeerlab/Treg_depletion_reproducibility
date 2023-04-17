## run_BP_serial.R

# script to launch from command line for batch processing of BayesPrism 

#### load libraries ####
require(BayesPrism)
require(Seurat)
require(dplyr)


#### load scRNA-seq data ####

## this load function contains objects:
#' ref.dat.filtered = cell x gene expression matrix
#' meta = metadata for cells in the expression matrix including cell type labels and sub-state labels
#' cell state labels were assigned such that each cell population had at least 30 marker genes
load("/data/peer/sam/treg_dep/mouse/data/scrna.herman.allCT.deconvRiter1.pUp0.05.rdata")

# assign to appropriate objects
cell.type.labels <- as.character(meta$cell.annot.labels)
cell.state.labels <- as.character(meta$cell.subtype.labels)
ref.dat.sig <- ref.dat.filtered

# assign unique barcode index to each cell
uniq_index = make.unique(rownames(ref.dat.sig), sep = "-")
rownames(ref.dat.sig) = uniq_index
rownames(meta) = uniq_index

#### get command line arguments ####
args = commandArgs(trailingOnly = T)

if(length(args) != 3){
  stop("enter the names of serial sections and the out directory")
}

sample_names = args[c(1,2)] %>% sort()
#sample_names = c('Ctrl1_A1', 'Ctrl1_B1')

outdir = args[3]

# make out dir if it doesn't exit
if(!(dir.exists(outdir))){
  dir.create(outdir, recursive = T)
}

#### load visium data ####
## load 2 technical replicates at the same time for deconvolution
# below folder has subfolders containing all samples
visium_dat_dir = '/data/peer/sam/treg_dep/visium/data/requant/results/'

visium_samples = list.dirs(visium_dat_dir, recursive = F) %>% basename()
visium_samples = visium_samples[grepl("IGO", visium_samples)]
visium_samples_names = gsub("_IGO_.+", "",  visium_samples)

# get visium samples to load
visium_samples_run = visium_samples[grepl(paste(sample_names, collapse = "|"), visium_samples)] %>% sort()

# get path for SpaceRanger output for each sample
visium_samples_files = lapply(visium_samples_run, function(x){
  file.path(visium_dat_dir, x, 'outs')
})



print(paste("Running for samples", visium_samples_run))

# read in visium data using Seurat
d = lapply(visium_samples_files, Load10X_Spatial)
names(d) = sample_names

# set the section name
d = lapply(names(d), function(x){
  r = d[[x]]
  r@meta.data$section = x
  return(r)
})
names(d) = sample_names


# merge sections together for deconvolution
d.m = merge(d[[1]], d[[2]], add.cell.ids = sample_names)

lib_size_cutoff = 1000

# calculate library size
library_size = apply(d.m@assays$Spatial@counts, 2, sum)
low_lib_size_spot = library_size < lib_size_cutoff

# remove spots with library size below threshold
d.m = d.m[,!(low_lib_size_spot)]
print(paste(sum(low_lib_size_spot), "spots removed due to lib size less than", lib_size_cutoff))

# index of samples
d_idx = d.m$section

#### gene manipulation ####
# this step specific for this dataset
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
# GFP is expressed under the Foxp3 promoter in this system
X.filtered[,'Foxp3'] = X.filtered[,'Foxp3'] + X.filtered[,'sDTR-eGFP']
X.filtered = X.filtered[,colnames(X.filtered) != 'sDTR-eGFP']

# change to title case for downstream matching
colnames(X.filtered) = stringr::str_to_title(colnames(X.filtered))

# create BP object
myPrism <- new.prism(
  reference=ref.dat.sig, 
  mixture=X.filtered,
  input.type="count.matrix", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = cell.state.labels,
  key=NULL,
  outlier.cut=0.01,
  outlier.fraction=0.1,
)

# run BP
bp.res <- run.prism(prism = myPrism, n.cores=30)

print(bp.res)

# save this file
saveRDS(bp.res, file = file.path(outdir, paste0(paste(sample_names, collapse = "_"),  "_bpRes.rds")))

# save cell fractions
theta = get.fraction (bp=bp.res,
                      which.theta="final",
                      state.or.type="type")
write.csv(theta, file =  file.path(outdir, paste(paste(sample_names, collapse = '_'),  "theta.csv", sep = "_")), 
          quote = F, row.names = T)

# save posterior CV
theta.cv <- bp.res@posterior.theta_f@theta.cv
write.csv(theta.cv, file =  file.path(outdir, paste(paste(sample_names, collapse = '_'),  "thetaCV.csv", sep = "_")), 
          quote = F, row.names = T)

# save index of which section is which
write.csv(data.frame(section = d_idx), file = file.path(outdir, paste(paste(sample_names, collapse = '_'),  "deconvSectionIndex.csv", sep = "_")), 
          row.names = F, quote = F)


#' Example command from the shell to run this script:
#'  Rscript --vanilla run_BP_serial.R Ctrl1_A1 Ctrl1_B1 $OUT_DIR

