library(recount3)
library(edgeR)
library(ggplot2)
library(recount)
library(ggpubr)
library(rstatix)

rse_brain <- readRDS("rse_brain.RDS")
rse_brain

rse_lung <- readRDS("rse_lung.RDS")
rse_lung

rse_liver <- readRDS("rse_liver.RDS")
rse_liver

assays(rse_brain)$counts <- transform_counts(rse_brain)
assays(rse_lung)$counts <- transform_counts(rse_lung)
assays(rse_liver)$counts <- transform_counts(rse_liver)

sum(assays(rse_brain)$counts)
sum(assays(rse_lung)$counts)
sum(assays(rse_liver)$counts)

table(rowData(rse_brain)$gbkey)

table(rowData(rse_lung)$gbkey)

table(rowData(rse_liver)$gbkey)


##RIN
colData(rse_brain)$gtex.smrin[46]
colData(rse_brain)$gtex.smrin[47]
colData(rse_brain)$gtex.smrin[52]

colData(rse_lung)$gtex.smrin[46]
colData(rse_lung)$gtex.smrin[47]
colData(rse_lung)$gtex.smrin[52]

colData(rse_liver)$gtex.smrin[46]
colData(rse_liver)$gtex.smrin[47]
colData(rse_liver)$gtex.smrin[52]


##% of rRNA

colData(rse_brain)$gtex.smrrnart[46]
colData(rse_brain)$gtex.smrrnart[47]
colData(rse_brain)$gtex.smrrnart[52]

colData(rse_lung)$gtex.smrrnart[46]
colData(rse_lung)$gtex.smrrnart[47]
colData(rse_lung)$gtex.smrrnart[52]

colData(rse_liver)$gtex.smrrnart[46]
colData(rse_liver)$gtex.smrrnart[47]
colData(rse_liver)$gtex.smrrnart[52]


##% of mapped reads
```{r}
colData(rse_brain)$"recount_qc.star.uniquely_mapped_reads_%_both"[46]
colData(rse_brain)$"recount_qc.star.uniquely_mapped_reads_%_both"[47]
colData(rse_brain)$"recount_qc.star.uniquely_mapped_reads_%_both"[52]

colData(rse_lung)$"recount_qc.star.uniquely_mapped_reads_%_both"[46]
colData(rse_lung)$"recount_qc.star.uniquely_mapped_reads_%_both"[47]
colData(rse_lung)$"recount_qc.star.uniquely_mapped_reads_%_both"[52]

colData(rse_liver)$"recount_qc.star.uniquely_mapped_reads_%_both"[46]
colData(rse_liver)$"recount_qc.star.uniquely_mapped_reads_%_both"[47]
colData(rse_liver)$"recount_qc.star.uniquely_mapped_reads_%_both"[52]


##replicates 46,47,52
rse_brain_selected <- rse_brain[,c(46,47,52)]

rse_lung_selected <- rse_lung[,c(46,47,52)]

rse_liver_selected <- rse_liver[,c(46,47,52)]

counts_brain_selected <- assays(rse_brain_selected)$counts

counts_lung_selected <- assays(rse_lung_selected)$counts

counts_liver_selected <- assays(rse_liver_selected)$counts

x <- cbind(counts_brain_selected,counts_lung_selected,counts_liver_selected)

colnames(x) <- c("Brain46", "Brain47","Brain52","Lung46","Lung47","Lung52","Liver46","Liver47","Liver52")

rownames(x) <- rowData(rse_brain_selected)$gene_name

y <- DGEList(counts=x)

group <- as.factor(c("Brain","Brain","Brain","Lung","Lung","Lung","Liver","Liver","Liver"))

y$samples$group <- group

y$samples$rin <- as.factor(c(colData(rse_brain_selected)$gtex.smrin,colData(rse_lung_selected)$gtex.smrin,colData(rse_liver_selected)$gtex.smrin))

y$samples$slice <- as.factor(c(colData(rse_brain_selected)$gtex.smtsd,colData(rse_lung_selected)$gtex.smtsd,colData(rse_liver_selected)$gtex.smtsd))

y$samples$sex <- as.factor(c(colData(rse_brain_selected)$gtex.sex,colData(rse_lung_selected)$gtex.sex,colData(rse_liver_selected)$gtex.sex))

y$samples$age <- as.factor(c(colData(rse_brain_selected)$gtex.age,colData(rse_lung_selected)$gtex.age,colData(rse_liver_selected)$gtex.age))

y$samples$rRNA <- as.factor(c(colData(rse_brain_selected)$gtex.smrrnart,colData(rse_lung_selected)$gtex.smrrnart,colData(rse_liver_selected)$gtex.smrrnart))

y$samples$mapped <- as.factor(c(colData(rse_brain_selected)$"recount_qc.star.uniquely_mapped_reads_%_both",colData(rse_lung_selected)$"recount_qc.star.uniquely_mapped_reads_%_both",colData(rse_liver_selected)$"recount_qc.star.uniquely_mapped_reads_%_both"))

y$samples$chrm <- as.factor(c(colData(rse_brain_selected)$"recount_qc.aligned_reads%.chrm",colData(rse_lung_selected)$"recount_qc.aligned_reads%.chrm",colData(rse_liver_selected)$"recount_qc.aligned_reads%.chrm"))

y

table(rowSums(y$counts==0)==9)

keep.exprs <- filterByExpr(y, group=group)

y <- y[keep.exprs,, keep.lib.sizes=FALSE]
dim(y)

logcpm_before <- cpm(y, log=TRUE)

y <- calcNormFactors(y, method = "TMM")

y

boxplot(logcpm_before)

logcpm_after <- cpm(y, log=TRUE)

boxplot(logcpm_after)

design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(y$samples$group)
design

logcpm <- cpm(y, log=TRUE)
plotMDS(logcpm, labels=group)

plotMDS(logcpm, labels=y$samples$rRNA)

plotMDS(logcpm, labels=y$samples$chrm)

plotMDS(logcpm, labels=y$samples$mapped)

plotMDS(logcpm, labels=y$samples$rin)

plotMDS(logcpm, labels=y$samples$slice)

plotMDS(logcpm, labels=y$samples$sex)

plotMDS(logcpm, labels=y$samples$age)


y <- estimateDisp(y, design)
plotBCV(y)


fit <- glmQLFit(y, design)
fit


#lung (top) vs brain (bottom)
qlfLungB <- glmQLFTest(fit, contrast=c(-1,0,1))
#liver (top) vs brain (bottom)
qlfLivB <- glmQLFTest(fit, contrast=c(-1,1,0))
#liver (top) vs lung (bottom)
qlfLivLung <- glmQLFTest(fit, contrast=c(0,1,-1))

topTags(qlfLungB, n=10,adjust.method = "BH", sort.by = "PValue")

resultsLungB <- topTags(qlfLungB, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
write.table(resultsLungB, "resultsLungB.txt")

resultsLivB <- topTags(qlfLivB, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
write.table(resultsLivB, "resultsLivB.txt")

resultsLivLung <- topTags(qlfLivLung, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
write.table(resultsLivLung, "resultsLivLung.txt")


summary(decideTests(qlfLungB, adjust.method="BH", p.value=0.01,lfc=0))
summary(decideTests(qlfLivB, adjust.method="BH", p.value=0.01, lfc=0))
summary(decideTests(qlfLivLung, adjust.method="BH", p.value=0.01, lfc=0))


##Genes over-represented in brain with respect to lung and liver
brain_BrainvsLung <- rownames(as.data.frame(resultsLungB)[as.data.frame(resultsLungB)$logFC < -1 & as.data.frame(resultsLungB)$FDR < 0.01 & as.data.frame(resultsLungB)$logCPM > 0 ,])
brain_BrainvsLung

brain_BrainvsLiver <- rownames(as.data.frame(resultsLivB)[as.data.frame(resultsLivB)$logFC < -1 & as.data.frame(resultsLivB)$FDR < 0.01 & as.data.frame(resultsLivB)$logCPM > 0 ,])

brain_vs_LungLiver <- intersect(brain_BrainvsLung, brain_BrainvsLiver)
head(brain_vs_LungLiver)




starts = c('LOC', 'LINC', 'MIR', 'SNORD', 'RPL')

for(s in starts) {
  brain_vs_LungLiver <- brain_vs_LungLiver[which(!startsWith(brain_vs_LungLiver, s))]
}

write.table(brain_vs_LungLiver, "Brain_vs_LungLiver.txt")



##Kruskal test and boxplot
which(rowData(rse_brain)$gene_name == "FABP7")
assays(rse_brain)$TPM <- recount::getTPM(rse_brain)
assays(rse_lung)$TPM <- recount::getTPM(rse_lung)
assays(rse_liver)$TPM <- recount::getTPM(rse_liver)

df_lung=data.frame(TPM=assays(rse_lung)$TPM[44843,],group="Lung")
df_br=data.frame(TPM=assays(rse_brain)$TPM[44843,],group="Brain")
df_liv=data.frame(TPM=assays(rse_liver)$TPM[44843,],group="Liver")
data_mrg=rbind(df_br,df_lung,df_liv)

res_kruskal=data_mrg %>% kruskal_test(TPM ~ group)
res_kruskal

pwc2=data_mrg %>% wilcox_test(TPM ~ group, p.adjust.method = "BH")
pwc2

pwc = pwc2 %>% add_xy_position(x = "group")

ggboxplot(data_mrg, x = "group", y = "TPM",outlier.shape = NA,width = 0.5,title="FABP7 expression across organs")


##Over-represented in lung
Lung_BrainvsLung <- rownames(as.data.frame(resultsLungB)[as.data.frame(resultsLungB)$logFC > 1 & as.data.frame(resultsLungB)$FDR < 0.01 & as.data.frame(resultsLungB)$logCPM > 0 ,])
Lung_LivervsLung <- rownames(as.data.frame(resultsLivLung)[as.data.frame(resultsLivLung)$logFC < -1 & as.data.frame(resultsLivLung)$FDR < 0.01 & as.data.frame(resultsLivLung)$logCPM > 0 ,])

Lung_vs_brainliver <- intersect(Lung_BrainvsLung, Lung_LivervsLung)
Lung_vs_brainliver 

for(s in starts) {
  Lung_vs_brainliver <- Lung_vs_brainliver[which(!startsWith(Lung_vs_brainliver, s))]
}

write.table(Lung_vs_brainliver, "Lung.txt")


##Over-represented in liver
Liver_BrainvsLiver <- rownames(as.data.frame(resultsLivB)[as.data.frame(resultsLivB)$logFC > 1 & as.data.frame(resultsLivB)$FDR < 0.01 & as.data.frame(resultsLivB)$logCPM > 0 ,])
Liver_LivervsLung <- rownames(as.data.frame(resultsLivLung)[as.data.frame(resultsLivLung)$logFC > 1 & as.data.frame(resultsLivLung)$FDR < 0.01 & as.data.frame(resultsLivLung)$logCPM > 0 ,])

Liver_vs_brainlung <- intersect(Liver_BrainvsLiver, Liver_LivervsLung)
Liver_vs_brainlung 

for(s in starts) {
  Liver_vs_brainlung <- Liver_vs_brainlung[which(!startsWith(Liver_vs_brainlung, s))]
}

write.table(Liver_vs_brainlung, "Liver.txt")