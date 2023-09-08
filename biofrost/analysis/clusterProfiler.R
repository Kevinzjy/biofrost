# "Enrichment analysis of differentially expression genes"
if (FALSE) {
    library("renv")
    renv::init("/data/public/software/sleuth")
    install.packages("optparse")
    install.packages("BiocManager")
    BiocManager::install('clusterProfiler')
    BiocManager::install("org.Mm.eg.db")
    BiocManager::install("org.Hs.eg.db")
}

rm(list=ls())
suppressPackageStartupMessages(library(renv))
renv::activate("/data/public/software/clusterProfiler")
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(org.Hs.eg.db))

help_info <- c("ENSEMBL", "ENTREZID", "GENENAME", "REFSEQ", "SYMBOL", "UNIPROT")

# Parse parameters
parser <- OptionParser()
parser <- add_option(parser, c("--input"), action="store", default=NA, type='character',
                     help="input gene list")
parser <- add_option(parser, c("--output"), action="store", default=NA, type='character',
                     help="output csv file")
parser <- add_option(parser, c("--db"), action="store", default='Hg', type='character',
                     help="reference database")
parser <- add_option(parser, c("--type"), action="store", default="SYMBOL", type='character',
                     help=paste(help_info, collapse = '; '))
parser <- add_option(parser, c("--analysis"), action="store", default='GO', type='character',
                     help="GO (Gene Ontology); KEGG")
opt <- parse_args(parser)

# main point of program is here, do this whether or not "verbose" is set
if (is.na(opt$input) || is.na(opt$output)) {
    cat("Please specify --in/--out, refer to the manual for detailed instruction!\n",
        file=stderr())
    quit()
}
if (!(opt$db %in% c("Hs", "Mm"))) {
    cat("Please specify --in/--out, refer to the manual for detailed instruction!\n",
        file=stderr())
    quit()
}

# Load gene list
gene_list <- as.character(read.csv(opt$input, header = FALSE)$V1)
db = paste(c("org", opt$db, "eg", "db"), collapse=".")
gene_id <- bitr(gene_list, fromType=opt$type, toType="ENTREZID", OrgDb=db)

if (opt$analysis == "GO") {
    # GO analysis
    bp <- enrichGO(gene_id$ENTREZID, OrgDb = db, ont="BP", keyType = 'ENTREZID')
    mf <- enrichGO(gene_id$ENTREZID, OrgDb = db, ont="MF", keyType = 'ENTREZID')
    cc <- enrichGO(gene_id$ENTREZID, OrgDb = db, ont="CC", keyType = 'ENTREZID')

    # Merge dataframe
    bp <- data.frame(simplify(bp))
    if (dim(bp)[1] > 0) {
        bp$ONTOLOGY = "BP"
    }
    mf <- data.frame(simplify(mf))
    if (dim(mf)[1] > 0) {
        mf$ONTOLOGY = "MF"
    }
    cc <- data.frame(simplify(cc))
    if (dim(cc)[1] > 0) {
        cc$ONTOLOGY = "CC"
    }
    res <- rbind(bp, mf, cc)

} else if (opt$analysis == "KEGG") {
    if (opt$db == "Hs") {
        organism <- "hsa"
    } else if (opt$db == "Mm") {
        organism <- "mmu"
    }
    kk <- enrichKEGG(gene = gene_id$ENTREZID, organism=organism)
    res <- data.frame(kk)
}

write.csv(res, file = opt$output)