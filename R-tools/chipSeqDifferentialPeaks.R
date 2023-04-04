# Generates 
# Paul Stretenowich - Apr 2021 - From Alain Pacis script

mLoad <- function(...) { sapply(sapply(match.call(), as.character)[-1], require, character.only = TRUE) }
mLoad(limma, edgeR, statmod, sva, qvalue, RColorBrewer, plyr, ggplot2, reshape2, gplots, Vennerable, ComplexHeatmap, circlize, ggrepel, enrichR, biomaRt, tibble, dplyr, stringr, magrittr)

## input and remove lowly-expressed genes
counts = read.table("/Users/pstretenowich/Mount_points/abacus/home/scratch/Kazak_ATAC-seq/allSamples_ATACseq.paired.noDup.counts", header = T, sep="\t", row.names=1)
dim(counts)

counts <- counts[, -grep("b", colnames(counts))]

## annotation & design
## normalization
anno <- data.frame("Condition" = c(rep("30C_Veh",3),rep("30C_PBZ",3),rep("6C_Veh",3),rep("6C_PBZ",3)), "Batch" = c(rep(c("B1","B2","B2"),4)))
rownames(anno) <- colnames(counts)
anno$Condition <- factor(anno$Condition)
anno$Condition <- relevel(anno$Condition, "30C_Veh")
anno$Batch <- factor(anno$Batch)

dge <- DGEList(counts=counts)
dge <- calcNormFactors(dge)
libsizes=dge$samples$lib.size
dge$samples

design=model.matrix(~ 0 + Batch + Condition, data=anno)
is.fullrank(design)
colnames(design) <- make.names(colnames(design))
v <- voom(dge,design,plot=FALSE)
expr = v$E

## PCA
x = prcomp(as.matrix(t(expr)))
perc_pc <- x$sdev^2/sum(x$sdev^2)*100

pca = x$x
pca = data.frame(pca)
pca[,"Sample"] <- rownames(pca)
#pca[,"RepID"] <- anno$RepID
pca[,"Condition"] <- anno$Condition
pca$Condition <- factor(pca$Condition, levels = c("30C_Veh","30C_PBZ","6C_Veh","6C_PBZ"))

i=1; j=2

pdf("PCA.pdf")
p1 <- ggplot(pca, aes(x = eval(parse(text=paste0("PC",i))), y = eval(parse(text=paste0("PC",j))), colour = Condition, label = Sample)) + 
geom_point(stroke = 1, size=2, alpha=1) + 
scale_color_manual(values = c("dodgerblue3","firebrick3","green4","darkgoldenrod1")) + 
geom_text_repel(size  = 2, point.padding = 0.25, segment.size  = 0.25, segment.color = "gray", direction = "both") +
xlab(paste0("PC",i,": ", round(perc_pc[[i]],2), "%")) +
ylab(paste0("PC",j,": ", round(perc_pc[[j]],2), "%")) +
theme_bw() +
theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
theme(aspect.ratio=1)
print(p1)

dev.off()

## ------------------------------------------------------
design=model.matrix(~ 0 + Condition + Batch, data=anno)
is.fullrank(design)
colnames(design) <- make.names(colnames(design))
v <- voom(dge,design,plot=FALSE)
expr = v$E

fit <-lmFit(v,design)
fit <- eBayes(fit)

contr.matrix <- makeContrasts(
   Comp_30C_PBZ_vs_30C_Veh = Condition30C_PBZ - Condition30C_Veh, 
   Comp_6C_Veh_vs_30C_Veh = Condition6C_Veh - Condition30C_Veh, 
   Comp_6C_PBZ_vs_30C_Veh = Condition6C_PBZ - Condition30C_Veh, 
   Comp_6C_PBZ_vs_6C_Veh = Condition6C_PBZ - Condition6C_Veh, 
   levels = design)
contr.matrix

fit1 <- contrasts.fit(fit, contrasts=contr.matrix)
fit1 <- eBayes(fit1)
data.frame(colnames(fit1))
t(summary(decideTests(fit1, p.value = 0.05, lfc = 0)))

res.30C_PBZ <- data.frame(topTable(fit1, coef="Comp_30C_PBZ_vs_30C_Veh", adjust="BH", n=nrow(counts), sort.by="none"))[,c("logFC", "P.Value", "adj.P.Val")]
res.6C_Veh <- data.frame(topTable(fit1, coef="Comp_6C_Veh_vs_30C_Veh", adjust="BH", n=nrow(counts), sort.by="none"))[,c("logFC", "P.Value", "adj.P.Val")]
res.6C_PBZ <- data.frame(topTable(fit1, coef="Comp_6C_PBZ_vs_30C_Veh", adjust="BH", n=nrow(counts), sort.by="none"))[,c("logFC", "P.Value", "adj.P.Val")]
res.6C_PBZ_vs_6C_Veh <- data.frame(topTable(fit1, coef="Comp_6C_PBZ_vs_6C_Veh", adjust="BH", n=nrow(counts), sort.by="none"))[,c("logFC", "P.Value", "adj.P.Val")]

res <- cbind(res.30C_PBZ, res.6C_Veh, res.6C_PBZ, res.6C_PBZ_vs_6C_Veh)
colnames(res) <- paste(rep(c("log2FC","pval","fdr"),3), c(rep("30C_PBZ",3),rep("6C_Veh",3),rep("6C_PBZ",3),rep("6C_PBZ_vs_6C_Veh",3)), sep=".")

## pval distribution
pdf("pval.pdf")
hist(res$pval.30C_PBZ, breaks = 100, xlab = "P-value", main = "30C_PBZ")
hist(res$pval.6C_Veh, breaks = 100, xlab = "P-value", main = "6C_Veh")
hist(res$pval.6C_PBZ, breaks = 100, xlab = "P-value", main = "6C_PBZ")
hist(res$pval.6C_PBZ_vs_6C_Veh, breaks = 100, xlab = "P-value", main = "6C_PBZ_vs_6C_Veh")
dev.off()

## output DEG names
lfc.cutoff = 1
fdr.cutoff = 0.01

deg.30C_PBZ <- rownames(res[which(res$fdr.30C_PBZ < fdr.cutoff & abs(res$log2FC.30C_PBZ) > lfc.cutoff),])
deg.6C_Veh <- rownames(res[which(res$fdr.6C_Veh < fdr.cutoff & abs(res$log2FC.6C_Veh) > lfc.cutoff),])
deg.6C_PBZ <- rownames(res[which(res$fdr.6C_PBZ < fdr.cutoff & abs(res$log2FC.6C_PBZ) > lfc.cutoff),])

deg.6C_PBZ_vs_6C_Veh <- rownames(res[which(((res$fdr.6C_Veh < fdr.cutoff & abs(res$log2FC.6C_Veh) > lfc.cutoff) | (res$fdr.6C_PBZ < fdr.cutoff & abs(res$log2FC.6C_PBZ) > lfc.cutoff))),])

length(deg.30C_PBZ)
length(deg.6C_Veh)
length(deg.6C_PBZ)
length(deg.6C_PBZ_vs_6C_Veh)
