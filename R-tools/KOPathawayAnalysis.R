# Differential Expression and Pathway Analysis for Seq2Fun using fgsea for pathway analysis
# August 25, 2021
# Pubudu Nawarathna


# Inputs: universal pathway list, edger_results.csv, fdr cutoff
# Usage : Rscript hicrep.R -s1 sample1_name -s2 sample2_name -f1 file_path1 -f2 file_path2 -c chromosome -o output_file_path -sm h/smooth_value -r resolution/bin size 
# -b bundary -w save_weights -cor save_correlation_matirx

# Usage

usage = function(errM) {
  cat("\nUsage : Rscript KOPathawayAnalysis.R [option] <Value>\n")
  cat("       -i         : edgeR report from genpipes\n")
  cat("       -map       : KO pathway list for the organism\n")
  cat("       -fdr       : FDR cutoff\n")
  cat("       -p         : output file prefix\n")
  cat("       -o         : output dir\n")
  #cat("       -html      : seq2fun output html file path\n")
  cat("       -rds       : KEGG library file path\n")
  cat("       -kegg      : path of file with all KEGG pathways\n")
  cat("       -h         : this help\n")
  
  stop(errM)
}


set.seed(123456789)

library(tidyr)
library(fgsea)
library(reshape2)
library(ggplot2)
library(ggridges)
library(R2HTML)
library(data.table)
library(dplyr)
#library(KEGGREST)

#download all pathway list from KEGG 
#http://rest.kegg.jp/list/pathway

# get differential expression analysis results

pathway_analysis <- function(degs, userpathways, fdr=0.05 , output_file_prefix , kegg_rds, output_dir,  all_kegg_pathways){
#degs <- fread("/lustre03/project/6007512/pubudu/seq2fun/testing_august_2021/differential_expression/seq2fun/H1ESC_GM12787/edger_results.csv")
#log2fc = 1
  degs <- fread(degs)
degs$log2FC <- abs(degs$log_FC)
degs_temp <- degs[,c(1,3)]
colnames(degs_temp) <- c("KO", "log2FC")

colnames(degs)[c(5,6)] <- c("edger.p_value", "edger.adj.p_value")

sig.genes <- subset(degs, edger.adj.p_value < fdr )

#pathway_list <- fread("KEGG_all_pathways.txt", header =F)
pathway_list <- fread(all_kegg_pathways, header =F)
colnames(pathway_list) <- c("mapID","Name")

# get pathways - uncomment this if you need to recreate the kegg.rds file with updated KEGG database

# kegg_ko <- lapply(unique(pathway_list$mapID), fun)
# names(kegg_ko) <- (pathway_list$mapID)
# kegg <- lapply(kegg_ko, substring, 4)
# 
# saveRDS(kegg, "kegg.rds")

#store this file in cvmfs

#kegg_rds="kegg.rds"
kegg <- readRDS(kegg_rds)
#userpathways = "user_pathways.txt"
user_pathways <- fread(userpathways, header=F)



kegg <- kegg[user_pathways$V1]
#gs.list <- list(kegg = kegg)


# do gene set analysis
degs.log2FC <- degs$log2FC
names(degs.log2FC) <- degs$id
degs.log2FC <- degs.log2FC[order(degs.log2FC)]

# calcualte gsea results
res.gsea <- fgsea(pathways = kegg, 
                  stats    = degs.log2FC,
                  minSize  = 5,
                  scoreType = "pos",
                  eps=0)

res.gsea <- merge(pathway_list, res.gsea, by.x = "mapID", by.y = "pathway", all = FALSE)
res.gsea <- res.gsea[order(res.gsea$pval), ]



# calculate ora results
res.gsoa <- fora(pathways = kegg,
                 genes = sig.genes,
                 universe = names(degs.log2FC),
                 minSize = 10,
                 maxSize = 500)
res.gsoa <- merge(pathway_list, res.gsoa, by.y = "pathway", by.x = "mapID", all = FALSE)
res.gsoa <- res.gsoa[order(res.gsoa$pval), ]

# gsoa results
res.gsoa <- res.gsoa[,c(1:4)]
colnames(res.gsoa) <- c("map", "name", "pval", "adj.pval")


# gsea results
res.gsea <- res.gsea[,c(1:4)]
colnames(res.gsea) <- c("map", "name", "pval", "adj.pval")

gsea.table <- res.gsea
# get sig gsea results
gsea.sig <- res.gsea[res.gsea$adj.pval < fdr, ]



# prepare data for plotting
degs$KO <- degs$id

kegg.plot <- reshape2::melt(kegg)
colnames(kegg.plot) <- c("KO", "map")

df <- merge(gsea.sig, kegg.plot, by = "map", all.x = TRUE, all.y = FALSE)
df <- merge(df, degs_temp, by = "KO", all.x = TRUE, all.y = FALSE)
df <- na.omit(df)

gsea.df <- merge(gsea.table, kegg.plot, by = "map", all.x = TRUE, all.y = FALSE)
gsea.df <- merge(gsea.df, degs_temp, by = "KO", all.x = TRUE, all.y = FALSE)
gsea.df <- na.omit(gsea.df)


write.table(gsea.df, file=paste0(output_dir,"/" ,output_file_prefix, ".txt"), sep="\t", col.names=T, row.names=F, quote=F)

#order by fold change
log2c_temp <- df[order(df$log2FC, decreasing = FALSE), ]

# make the plot
png(file=paste0(output_dir,"/" ,output_file_prefix, ".png"), width = 500, height = 350)
figure <- ggplot(df, aes(x = log2FC, y = name, fill = adj.pval)) +
  geom_density_ridges(
    jittered_points = TRUE, point_shape = "|", point_size = 2, point_color = "#898A89",
    color = "white",
    scale = 2, rel_min_height = .01, size = 0.25,
    position = position_points_jitter(height = 0)) +
  scale_y_discrete(expand = c(0, 0), name = "KEGG pathway") +
  scale_x_continuous(expand = c(0, 0), name = "KO log2FC") +
  scale_fill_gradient("Adj. p-value",
                      low = '#ffa34d',
                      high = "#994a00")+
  coord_cartesian(clip = "off") +
  theme_ridges(center = TRUE) +
  theme(legend.position = "none",
        text = element_text(size=12, color = "black"),
        axis.title = element_text(size=12, face = "bold"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))
if(dim(df)[1]>0){
print(figure)
}
dev.off()

#HTMLInsertGraph(file=paste(html_path), GraphFileName = "seq2fun_ko_pathway.png", Caption="KO pathway analysis")
}


fun <- function(x) {
  as.character(keggLink("ko",x))
  
}

ARG = commandArgs(trailingOnly=T)

## default arg values
#no default args

## get arg variables
for (i in 1:length(ARG)) {
  if (ARG[i] == "-map") {
    map_list = ARG[i+1]
  } else if (ARG[i] == "-i") {
    degs = ARG[i+1]
  } else if (ARG[i] == "-fdr") {
    fdr = ARG[i+1]
  } else if (ARG[i] == "-p") {
    output_file_prefix = ARG[i+1]
  } else if (ARG[i] == "-rds") {
    kegg_rds = ARG[i+1]
  #} else if (ARG[i] == "-html") {
    #html_path = ARG[i+1]
  } else if (ARG[i] == "-o") {
    output_dir = ARG[i+1]
  } else if (ARG[i] == "-kegg") {
    kegg_all = ARG[i+1]
  }
}

## check arg consitency
if (!(file.exists(map_list))) {
  usage("Error : pathway list is not found")
}
if (!(file.exists(degs))) {
  usage("Error : edgeR differential expression file is not found")
}
if (fdr == "") {
  usage("Error : fdr is not specified")
}
if (output_file_prefix == "") {
  usage("Error : outputfile path is not specified")
}
if (!(file.exists(kegg_rds))) {
  usage("Error : KEGG library file is not found")
}
#if (!(file.exists(html_path))) {
 # usage("Error : seq2fun All_samples.html output file is not found")
#}
if (output_file_prefix == "") {
  usage("Error : outputfile prefix is not specified")
}
if (output_dir == "") {
  usage("Error : output directory is not specified")
}
if (!(file.exists(kegg_all))) {
  usage("Error : file with all KEGG pathways is not found")
}

print(paste0("pathway list=",map_list))
print(paste0("edgeR file=",degs))
print(paste0("FDR=",fdr))
#print(paste0("html file=",html_path))
print(paste0("output file prefix=",output_file_prefix))
print(paste0("output dir=",output_dir))
print(paste0("kegg all pathways=",kegg_all))
print(paste0("kegg library rds path=",kegg_rds))


pathway_analysis(degs, map_list, fdr, output_file_prefix , kegg_rds, output_dir, kegg_all)
