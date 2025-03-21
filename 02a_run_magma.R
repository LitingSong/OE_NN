# conda activate Py39_R43_Ju10
# run step1-2 using Rscript

.libPaths(c('/sc/arion/projects/CommonMind/liting/conda/envs/Py39_R43_Ju10/lib/R/library',.libPaths()))
#.libPaths(c('/sc/arion/projects/roussp01a/liting/software/conda/Monocle3/lib/R/library',.libPaths()))


.libPaths(c("/sc/arion/projects/roussp01a/jaro/programs/R_libs_4_2", .libPaths()))
source('/sc/arion/projects/roussp01a/jaro/atacseq_ad/NYGC_AD_R01/downstream_magma.R')
#source("bb/atacseq_ad/NYGC_AD_R01/downstream_magma.R")
#source("/sc/arion/projects/roussp01a/jaro/atacseq_ad/NYGC_AD_R01/downstream_greatR.R")
#load("/sc/arion/projects/roussp01b/resources/databases/gene-info/gene-sets-misc/gseafunctions_noChrY_noChrM_v3.Rdata");  #This is an output from offlineFisherGsea in downstream_greatR.R and has var standardGeneSets

library(ggplot2)
library(RColorBrewer)
library(dplyr)

setwd('/sc/arion/projects/roussp01a/liting/Olf')

#load('/sc/arion/projects/roussp01a/liting/Olf/data/DEG_trend_Olf_k6.RData')
load('/sc/arion/projects/roussp01a/liting/Olf/data/DEG_trend_Olf_k7_graph.RData')
load('/sc/arion/projects/roussp01a/liting/Olf/data/sub_tf_targets_rep_s.RData')
load('/sc/arion/projects/roussp01a/liting/Olf/data/deggenes.RData')



g_trend_list <- split(DEG_trend_Olf$gene, DEG_trend_Olf$trend_class)
g_cluster_list <-  split(DEG_trend_Olf$gene, DEG_trend_Olf$g_cluster)
g_tf_list <-  split(sub_tf_targets_rep_s$Gene, sub_tf_targets_rep_s$tf)

deggenes1 <- subset(deggenes, cat%in%names(table(deggenes$cat))[which(table(deggenes$cat)> 30)])
g_de_list <- split(deggenes1$gene,deggenes1$cat)
g_demerge_list <- split(deggenes$gene,deggenes$catmerge)    

# step1. Prepare genesets
ENSEMBL_INFO_hg38 = "/sc/arion/projects/roussp01b/resources/databases/gene-info/ensembl/my-downloads/muchEnsemblInfo_hg38.tsv.gz"

#make a table of gene-level ensembl info
ensemblInfo=read.delim(ENSEMBL_INFO_hg38, stringsAsFactors=F)
ensemblGeneInfo=unique(ensemblInfo[,c("Ensembl.Gene.ID", "Strand", "Gene.End..bp.", "Gene.Start..bp.", "Associated.Gene.Name", "Gene.type", "Status..gene.", "Source..gene.", "Transcript.count", "Description", "Chromosome.Name","HGNC.symbol")])

g_cluster_list_ensg <- lapply(g_cluster_list, function(x){unique(subset(ensemblGeneInfo,HGNC.symbol %in% x)[,"Ensembl.Gene.ID"])})
g_trend_list_ensg <- lapply(g_trend_list, function(x){unique(subset(ensemblGeneInfo,HGNC.symbol %in% x)[,"Ensembl.Gene.ID"])})
g_tf_list_ensg <- lapply(g_tf_list, function(x){unique(subset(ensemblGeneInfo,HGNC.symbol %in% x)[,"Ensembl.Gene.ID"])})
g_de_list_ensg <- lapply(g_de_list, function(x){unique(subset(ensemblGeneInfo,HGNC.symbol %in% x)[,"Ensembl.Gene.ID"])})
g_demerge_ensg <- lapply(g_demerge_list, function(x){unique(subset(ensemblGeneInfo,HGNC.symbol %in% x)[,"Ensembl.Gene.ID"])}) # not by trend, pnly oevrlap or nonoverlap


# step2. Run MAGMA
# myAnalysis <- calcMagmaGeneSetPvals(outDir="/sc/arion/projects/CommonMind/liting/ENT/magma_graph_cluster",
#                                     geneMetaSets=g_cluster_list_ensg,
#                                     pvalInfo=grabGenePvalInfo(selectedGeneTypeName="ensemblProtCodGenes35kb10kbAutosomesNoBmhc"),
#                                     checkIfOutDirExists=F)

# myAnalysis <- calcMagmaGeneSetPvals(outDir="/sc/arion/projects/CommonMind/liting/ENT/magma_graph_trend",
#                                     geneMetaSets=g_trend_list_ensg,
#                                     pvalInfo=grabGenePvalInfo(selectedGeneTypeName="ensemblProtCodGenes35kb10kbAutosomesNoBmhc"),
#                                     checkIfOutDirExists=F)

myAnalysis <- calcMagmaGeneSetPvals(outDir="/sc/arion/projects/CommonMind/liting/ENT/magma_tf_targets",
                                    geneMetaSets=g_tf_list_ensg,
                                    pvalInfo=grabGenePvalInfo(selectedGeneTypeName="ensemblProtCodGenes35kb10kbAutosomesNoBmhc"),
                                    checkIfOutDirExists=F)

# myAnalysis <- calcMagmaGeneSetPvals(outDir="/sc/arion/projects/CommonMind/liting/ENT/magma_de_overlap",
#                                     geneMetaSets=g_de_list_ensg,
#                                     pvalInfo=grabGenePvalInfo(selectedGeneTypeName="ensemblProtCodGenes35kb10kbAutosomesNoBmhc"),
#                                     checkIfOutDirExists=F)

# myAnalysis <- calcMagmaGeneSetPvals(outDir="/sc/arion/projects/CommonMind/liting/ENT/magma_demerge_overlap",
#                                     geneMetaSets=g_demerge_ensg,
#                                     pvalInfo=grabGenePvalInfo(selectedGeneTypeName="ensemblProtCodGenes35kb10kbAutosomesNoBmhc"),
#                                     checkIfOutDirExists=F)

