#sauger DAPc
library(rhdf5)
library(abind)
library(RColorBrewer)
library(vcfR)
library(dartR)
library(tidyverse)
library(ggpubr)
library(viridis)
library(matrixStats)

randomize.dapc <- function(gl.obj,pop,npca,niter=10,verbose=TRUE){
  gl.obj$pop <- pop
  orig.mat <- as.matrix(gl.obj)
  rand.loadings <- c()
  
  for (i in 1:niter){
    if(verbose==TRUE){
      print(paste("running iteration ",i,sep=""))
    }
    rand.mat <- apply(orig.mat, 2, function(x) sample(x, replace = T))
    rand.gl <- as.genlight(rand.mat)
    rm(rand.mat)
    indNames(rand.gl) <- indNames(gl.obj)
    ploidy(rand.gl) <- 2
    rand.gl$pop <- gl.obj$pop
    
    # remove any NA loci
    toRemove <- is.na(glMean(rand.gl, alleleAsUnit = FALSE)) # TRUE where NA
    #which(toRemove) # position of entirely non-typed loci
    glNoNA <- rand.gl[, !toRemove]
    
    dapc <- dapc(glNoNA,pop=pop(glNoNA),
                 n.pca=npca, n.da=3)
    loadings <- dapc$var.contr[,1]
    rm(rand.gl)
    #hist(loadings)
    #print(summary(loadings))
    
    rand.loadings <- append(rand.loadings,loadings)
  }
  #rand.hist <- hist(rand.loadings)
  #quantile(rand.loadings,c(0.95,0.975,0.99,1),na.rm=TRUE)
  return(quantile(rand.loadings,c(0.99),na.rm=TRUE))
}



'''
vcf <- read.vcfR("./../pop_struct/all_saug/vcf_file/sauger_miss0.5_maf0.03.recode.vcf")

geno <- read.geno("./../pop_struct/all_saug/sauger_miss0.5_maf0.03.recode.geno")
missing <- apply(geno,1,function(x) length(which(x == 9))/length(x))
keep <- which(missing < 0.7)
keep_cols <- c(1:9,keep+9)
vcf <- read.vcfR("./../pop_struct/all_saug/vcf_file/sauger_miss0.5_maf0.03.recode.vcf",cols=keep_cols)
write.vcf(vcf,file="./../pop_struct/all_saug/vcf_file/sauger_0.5_0.03_filt.vcf.gz")
'''

#start here
vcf <- read.vcfR("./../pop_struct/all_saug/vcf_file/sauger_miss0.5_maf0.03.recode.vcf")
gl <- vcfR2genlight(vcf)

meta <- read.csv("./../entropy/results_meta_sarwae_0.6_0.05.csv")
meta <- meta[which(meta$q < 0.1),]
vcf_names <- colnames(vcf@gt)[-1]
vcf_names <- gsub(".*aln_(EGM18_[[:digit:]]+).sorted.bam","\\1",vcf_names)
meta <- meta[which(meta$EGM_ID %in% vcf_names),]
meta <- meta[match(vcf_names,meta$EGM_ID),]

pop(gl) <- meta$Sex
gl@other$loc.metrics.flags$monomorphs <- TRUE
gl <- gl.drop.pop(gl,pop.list=c(NA,"I"))
gl <- gl.drop.pop(gl,pop.list=c("I"))
dapc <- dapc(gl,n.pca=150,n.da=2)
dapc_threshold <- randomize.dapc(gl,npca=150,pop=pop(gl))

levels(dapc$grp) <- c("Female","Male")
colors.vir <- c("#5f0f40","#fb8b24")
layout(matrix(c(1,2,2,2), nrow = 1, ncol = 4, byrow = TRUE))
scatter(dapc, scree.da=TRUE, 
        bg="white", pch=20, cell=0, 
        cstar=0, solid=0.4, cex=2, clab=0, 
        leg=TRUE, col=colors.vir)

loadingplot(dapc$var.contr[,1],
            cex.lab=0.75,lab="",srt=90,byfac=FALSE,
            xlab="SNP location", main="",threshold = 1)
chr_start <- min(grep("041335",gl$loc.names))
chr_end <- max(grep("041335",gl$loc.names))
abline(v=chr_start,lty=5,lwd=2,col="gray40")
abline(v=chr_end,lty=5,lwd=2,col="gray40")
text(chr_end+300,0.017,labels="Chr 7",cex=1.5)
abline(h=dapc_threshold, lty=2,lwd=2)


#cutoff 99% = 0.001389319 
dapc_threshold <- 0.001389319 

sex_loci <- data.frame(gl@loc.names[dapc$var.contr[,1] > dapc_threshold],which(dapc$var.contr[,1] > dapc_threshold))
colnames(sex_loci) <- c("CHROM","num")
sex_loci$chrom_mod <- gsub("(NC_[[:digit:]]+.1)_[[:digit:]]+","\\1",sex_loci$CHROM)
sex_loci$pos <- as.numeric(gsub("NC_[[:digit:]]+.1_([[:digit:]]+)","\\1",sex_loci$CHROM))
sex_loci$loading <- dapc$var.contr[which(dapc$var.contr[,1] > dapc_threshold),1]
write.csv(sex_loci,file="sex_loci_0.5_0.03.csv",row.names=F,quote=F)

all_locs <- data.frame(gl@loc.names)
colnames(all_locs)[1] <- "chrom"
all_locs$chrom_mod <- gsub("(NC_[[:digit:]]+.1)_[[:digit:]]+","\\1",all_locs$chrom)
all_locs$pos <- as.numeric(gsub("NC_[[:digit:]]+.1_([[:digit:]]+)","\\1",all_locs$chrom))


gl_rm <-  gl.drop.loc(gl,gl@loc.names[dapc$var.contr[,1] > dapc_threshold])


test_vcf <- vcf@gt[-as.vector(which(dapc$var.contr[,1] > dapc_threshold)),]
test_fix <- vcf@fix[-as.vector(which(dapc$var.contr[,1] > dapc_threshold)),]
vcf_filt <- vcf
vcf_filt@gt <- test_vcf
vcf_filt@fix <- test_fix

write.vcf(vcf_filt,file="./sauger_0.5_0.03_nosex.vcf.gz")

#try optimizing 
gl_mat <- as.matrix(gl)
gl_mat_noNA <- gtools::na.replace(gl_mat, mean, na.rm=T)
dapc_xval <- xvalDapc(gl_mat_noNA,grp=pop(gl),n.da=2,xval.plot=TRUE,center=TRUE,training.set=0.9, parallel="snow",ncpus=3)
#150 PCs





#try to figure out weird PC3 structure in all sauger

source("C:/users/wrose/Documents/Masters/scripts/pca_funcs.R")

file <- read.table("./../DAPC/pntest_mean_saug_0.5_0.03_nosex.txt",header=F)

vcf <- vcfR::read.vcfR("./../DAPC/sauger_0.5_0.03_nosex.vcf",nrows=1)

meta <- read.csv("./../entropy/results_meta_sarwae_0.6_0.05.csv",header=T,stringsAsFactors = F)

vcf_names <- colnames(vcf@gt)[-1]
vcf_names <- gsub(".*aln_(EGM18_[[:digit:]]+).sorted.bam","\\1",vcf_names)

meta_filt <- meta[which(meta$EGM_ID %in% vcf_names),]
meta_filt <- meta_filt[match(vcf_names,meta_filt$EGM_ID),]
fishinfo <- meta_filt
colnames(fishinfo)[4] <- "Tributary"

#LEA::vcf2geno("./../DAPC/sauger_0.5_0.03_nosex.vcf")
gt <- LEA::read.geno("./../DAPC/sauger_0.5_0.03_nosex.geno")
missing <- apply(gt,1,function(x) length(which(x == 9))/length(x))
file <- file[,which(missing < 0.7)]
fishinfo <- fishinfo[which(missing < 0.7),]
meta_filt <- meta_filt[which(missing < 0.7),]
names <- vcf_names[which(missing < 0.7)]


fishinfo$Tributary[grep("Bighorn",fishinfo$Tributary)] <- "Bighorn River"
fishinfo$Tributary[grep("Wind",fishinfo$Tributary)] <- "Wind River Basin"
fishinfo$Tributary[grep("Popo",fishinfo$Tributary)] <- "Wind River Basin"
fishinfo$Tributary <- gsub("Boysen","Boysen Reservoir",fishinfo$Tributary)
pca <- doPCA(file,fishinfo)

group1 <- names[which(pca[[1]]$PC3 < -0.01)]
group2 <- names[which(pca[[1]]$PC3 > 0)]
group_assignment <- rep(NA,length(names))
group_assignment[which(names %in% group1)] <- "1"
group_assignment[which(names %in% group2)] <- "2"


vcf <- read.vcfR("./sauger_0.5_0.03_nosex_filt.vcf")
meta <- read.csv("./../entropy/results_meta_sarwae_0.6_0.05.csv",header=T,stringsAsFactors = F)
vcf_names <- colnames(vcf@gt)[-1]
vcf_names <- gsub(".*aln_(EGM18_[[:digit:]]+).sorted.bam","\\1",vcf_names)

meta_filt <- meta[which(meta$EGM_ID %in% vcf_names),]
meta_filt <- meta_filt[match(vcf_names,meta_filt$EGM_ID),]
fishinfo <- meta_filt
colnames(fishinfo)[4] <- "Tributary"

gl <- vcfR2genlight(vcf)
gl@other$loc.metrics.flags$monomorphs <- TRUE
pop(gl) <- group_assignment
gl <- gl.drop.pop(gl,NA)
#gl_mat <- as.matrix(gl)
#gl_mat_noNA <- gtools::na.replace(gl_mat, mean, na.rm=T)
#dapc_xval <- xvalDapc(gl_mat_noNA,grp=pop(gl),n.da=2,xval.plot=TRUE,center=TRUE,training.set=0.9, parallel="snow",ncpus=3)
#350 PCs


dapc <- dapc(gl,n.pca=350,n.da=2)
dapc_df <- data.frame(gl$loc.names,dapc$var.contr[,1])
write.csv(dapc_df,file="./other_struct_dapc_all.csv",row.names=F,quote=F)
#dapc_threshold <- randomize.dapc(gl,npca=150,pop=pop(gl))
#threshold = 0.001455475
dapc_threshold <- 0.001455475

colors.vir <- c("#f16d22","#1d4e89")
layout(matrix(c(1,2,2,2), nrow = 1, ncol = 4, byrow = TRUE))
scatter(dapc, scree.da=TRUE, 
        bg="white", pch=20, cell=0, 
        cstar=0, solid=0.4, cex=2, clab=0, 
        leg=TRUE, col=colors.vir)

loadingplot(dapc$var.contr[,1],
            cex.lab=0.8,srt=90,byfac=FALSE,
            xlab="SNP location", main="", threshold=1)
abline(h=dapc_threshold, lty=2)

other_loci <- data.frame(gl@loc.names[dapc$var.contr[,1] > dapc_threshold],which(dapc$var.contr[,1] > dapc_threshold))
colnames(other_loci) <- c("CHROM","num")
other_loci$chrom_mod <- gsub("(NC_[[:digit:]]+.1)_[[:digit:]]+","\\1",other_loci$CHROM)
other_loci$pos <- as.numeric(gsub("NC_[[:digit:]]+.1_([[:digit:]]+)","\\1",other_loci$CHROM))
other_loci$loading <- dapc$var.contr[which(dapc$var.contr[,1] > dapc_threshold),1]
write.csv(other_loci,file="./other_struct_loci_0.5_0.03.csv",row.names = F,quote=F)

#Fst between mystery groups
group_assignment[which(is.na(group_assignment))] <- "hyb"
pop(gl) <- group_assignment
mystery_fst <- gl.fst.pop(gl)
#get 95% confidence interval
upper_confint <- mystery_fst$Bootstraps$`Upper bound CI limit`[1]
lower_confint <- mystery_fst$Bootstraps$`Lower bound CI limit`[1]

#look at heterozygosity for sex loci and other loci

sex_loci <- read.csv("./sex_loci_0.5_0.03.csv")
vcf <- vcfR::read.vcfR("./../DAPC/sauger_0.5_0.03_nosex.vcf",nrows=1)
meta <- read.csv("./../entropy/results_meta_sarwae_0.6_0.05.csv",header=T,stringsAsFactors = F)
vcf_names <- colnames(vcf@gt)[-1]
vcf_names <- gsub(".*aln_(EGM18_[[:digit:]]+).sorted.bam","\\1",vcf_names)

meta_filt <- meta[which(meta$EGM_ID %in% vcf_names),]
meta_filt <- meta_filt[match(vcf_names,meta_filt$EGM_ID),]
fishinfo <- meta_filt
colnames(fishinfo)[4] <- "Tributary"

gt <- LEA::read.geno("./sauger_0.5_0.03_nosex.geno")
missing <- apply(gt,1,function(x) length(which(x == 9))/length(x))
fishinfo <- fishinfo[which(missing < 0.7),]
meta_filt <- meta_filt[which(missing < 0.7 & meta_filt$Sex %in% c("M","F")),]
names <- vcf_names[which(missing < 0.7)]

gt_sex <- gt[which(meta_filt$Sex %in% c("M","F")),sex_loci$num]
sex_het <- apply(gt_sex,1,function(x) length(which(x == 1))/length(which(x != 9)))
sex_het.loc <- apply(gt_sex,2,function(x) length(which(x == 1))/length(which(x != 9)))
het.loc <- apply(gt,2,function(x) length(which(x == 1))/length(which(x != 9)))

gt_sex_m <- gt_sex[which(meta_filt$Sex == "M"),]
gt_sex_f <- gt_sex[which(meta_filt$Sex == "F"),]
sex_het_f <- apply(gt_sex_f,2,function(x) length(which(x == 1))/length(which(x != 9)))
sex_het_m <- apply(gt_sex_m,2,function(x) length(which(x == 1))/length(which(x != 9)))

f_lm <- lm(sex_het_f ~ sex_loci$loading)
m_lm <- lm(sex_het_m ~ sex_loci$loading)
summary(f_lm)
summary(m_lm)

gt_m <- gt[which(meta_filt$Sex == "M"),]
gt_f <- gt[which(meta_filt$Sex == "F"),]
het_f <- apply(gt_f,2,function(x) length(which(x == 1))/length(which(x != 9)))
het_m <- apply(gt_m,2,function(x) length(which(x == 1))/length(which(x != 9)))

t.test(het_f,het.loc)
t.test(het_m,het.loc) 
t.test(het_m,het_f)

#other loci
other_loci <- read.csv("./other_struct_loci_0.5_0.03.csv")
vcf <- vcfR::read.vcfR("./../DAPC/sauger_0.5_0.03_nosex.vcf",nrows=1)
meta <- read.csv("./../entropy/results_meta_sarwae_0.6_0.05.csv",header=T,stringsAsFactors = F)
vcf_names <- colnames(vcf@gt)[-1]
vcf_names <- gsub(".*aln_(EGM18_[[:digit:]]+).sorted.bam","\\1",vcf_names)

meta_filt <- meta[which(meta$EGM_ID %in% vcf_names),]
meta_filt <- meta_filt[match(vcf_names,meta_filt$EGM_ID),]
fishinfo <- meta_filt
colnames(fishinfo)[4] <- "Tributary"

gt <- LEA::read.geno("./sauger_0.5_0.03_nosex.geno")
missing <- apply(gt,1,function(x) length(which(x == 9))/length(x))
fishinfo <- fishinfo[which(missing < 0.7),]
meta_filt <- meta_filt[which(missing < 0.7),]
names <- vcf_names[which(missing < 0.7)]

group_assignment <- read.csv("./group_assignment_other.csv")
meta_filt <- cbind(meta_filt,group_assignment)
meta_other <- meta_filt[which(is.na(group_assignment$group_assignment) == FALSE),]

gt_other <- gt[which(is.na(meta_other$group_assignment) == FALSE),other_loci$num]
other_het <- apply(gt_other,2,function(x) length(which(x == 1))/length(x))

gt_other1 <- gt_other[which(meta_other$group_assignment == "1"),]
gt_other2 <- gt_other[which(meta_other$group_assignment == "2"),]

het_other1 <- apply(gt_other1,2,function(x) length(which(x == 1))/length(which(x != 9)))
het_other2 <- apply(gt_other2,2,function(x) length(which(x == 1))/length(which(x != 9)))

t.test(het_other1,het_other2)
t.test(het_other1,het.loc)
t.test(het_other2,het.loc)

lm1 <- lm(het_other1~other_loci$loading)
lm2 <- lm(het_other2~other_loci$loading)




#see if similar loci are seen in walleye
library(vcfR)
library(dartR)

meta <- read.csv("./../entropy/results_meta_sarwae_0.6_0.05.csv")
vcf <- read.vcfR("./sarwae_miss0.5_maf0.03_filt.vcf")
colnames(vcf@gt) <- gsub(".*aln_(EGM18_[[:digit:]]+).sorted.bam","\\1",colnames(vcf@gt))

meta <- meta[which(meta$EGM_ID %in% colnames(vcf@gt) ==TRUE),]
meta <- meta[match(colnames(vcf@gt)[-1], meta$EGM_ID),]

gl <- vcfR2genlight(vcf)
sp <- rep("sar",nrow(meta))
sp[which(meta$q > 0.9)] <- "wae"
pop(gl) <- sp
gl_wae <- gl.drop.pop(gl,"sar")

meta_wae <- meta[which(sp == "wae"),]
pop(gl_wae) <- meta_wae$Sex

gl_wae <- gl.drop.pop(gl_wae,NA)
gl_mat <- as.matrix(gl_wae)
gl_mat_noNA <- gtools::na.replace(gl_mat, mean, na.rm=T)
dapc_xval <- xvalDapc(gl_mat_noNA,grp=pop(gl_wae),n.da=2,xval.plot=TRUE,center=TRUE,training.set=0.9, parallel="snow",ncpus=3)
#20 PCs

dapc <- dapc(gl_wae,n.pca=20,n.da=2)
dapc_threshold <- randomize.dapc(gl_wae,npca=20,pop=pop(gl_wae))
dapc_threshold <- 0.001083235 

colors.vir <- plasma(6)
layout(matrix(c(1,2,2,2), nrow = 1, ncol = 4, byrow = TRUE))
scatter(dapc, scree.da=TRUE, 
        bg="white", pch=20, cell=0, 
        cstar=0, solid=0.4, cex=2, clab=0, 
        leg=TRUE, col=colors.vir)

loadingplot(dapc$var.contr[,1],
            cex.lab=0.5,srt=90,byfac=FALSE,
            xlab="SNP location", main="", threshold=dapc_threshold)
abline(h=dapc_threshold, lty=2)
#pretty shoddy


sex_loci <- read.csv("./sex_loci_0.5_0.03.csv",header=T)
all_loci <- paste(vcf@fix[,1],vcf@fix[,2],sep="_")
remove <- which(all_loci %in% sex_loci$CHROM)

vcf@fix <- vcf@fix[-remove,]
vcf@gt <- vcf@gt[-remove,]
write.vcf(vcf,file="./sarwae_0.5_0.03_nosex_filt.vcf.gz")



 