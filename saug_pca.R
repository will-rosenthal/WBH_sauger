
#Sauger PCA for detecting population structure
#re-do from last time -- cleaner script

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

plotPCApng(pca,"saug_more_pcs",fishinfo$Tributary)
#plotPCApng(pca,"report_pca",fishinfo$Tributary)

pal <- RColorBrewer::brewer.pal(n=3,name="Dark2")
point_pal <- c(21:23)
{png(filename="report_pca.png",width=6,height=6,units="in",res=300)
plot(pca[[1]]$PC1,pca[[1]]$PC2,pch=point_pal[as.factor(fishinfo$Tributary)],col=pal[as.factor(fishinfo$Tributary)],bg=adjustcolor(pal[as.factor(fishinfo$Tributary)],alpha.f = 0.6),xlab=paste0("PC1 (",round(pca[[2]]$importance[2,2]*100,digits=2),"%)"),ylab=paste0("PC2 (",round(pca[[2]]$importance[2,3]*100,digits=2),"%)"),sub="6,536 loci")
legend("bottomright", legend = unique(fishinfo$Tributary), pch=point_pal,col=pal[as.factor(unique(fishinfo$Tributary))], pt.bg=adjustcolor(pal[as.factor(unique(fishinfo$Tributary))],alpha.f = 0.6))
dev.off()
}


pal <- RColorBrewer::brewer.pal(n=3,name="Dark2")
{png(filename="report_pca3.png",width=6,height=6,units="in",res=300)
  plot(pca[[1]]$PC2,pca[[1]]$PC3,pch=21,col=pal[as.factor(fishinfo$Tributary)],bg=adjustcolor(pal[as.factor(fishinfo$Tributary)],alpha.f = 0.6),xlab=paste0("PC2 (",round(pca[[2]]$importance[2,1]*100,digits=2),"%)"),ylab=paste0("PC3 (",round(pca[[2]]$importance[2,2]*100,digits=2),"%)"),sub="6,536 loci")
  legend("bottomright", legend = unique(fishinfo$Tributary), pch=21,col=pal[as.factor(unique(fishinfo$Tributary))], pt.bg=adjustcolor(pal[as.factor(unique(fishinfo$Tributary))],alpha.f = 0.6))
  dev.off()
}


order <- findInterval(fishinfo$Length,sort(fishinfo$Length))
cols <- RColorBrewer::brewer.pal(n=10,"RdYlGn")
pal <- colorRampPalette(cols)


# pal <- RColorBrewer::brewer.pal(n=5,name="Set1")
{png(filename="saug_pca3_size.png",width=6,height=6,units="in",res=300)
  plot(bh_pca[[1]]$PC2,bh_pca[[1]]$PC3,pch=21,col=pal(nrow(fishinfo))[order],bg=adjustcolor(pal(nrow(fishinfo))[order],alpha.f = 0.6),xlab=paste0("PC2 (",round(bh_pca[[2]]$importance[2,2]*100,digits=2),"%)"),ylab=paste0("PC3 (",round(bh_pca[[2]]$importance[2,3]*100,digits=2),"%)"))
  legend("bottomleft", legend = c("8 cm","30 cm"), pch=21,col=pal(2), pt.bg=adjustcolor(pal(2),alpha.f = 0.6))
  dev.off()
}



{png(filename="report_pca_sex.png",width=6,height=6,units="in",res=300)
  plot(pca[[1]]$PC1,pca[[1]]$PC2,pch=21,col=pal[as.factor(fishinfo$Sex)],bg=adjustcolor(pal[as.factor(fishinfo$Sex)],alpha.f = 0.6),xlab=paste0("PC1 (",round(pca[[2]]$importance[2,1]*100,digits=2),"%)"),ylab=paste0("PC2 (",round(pca[[2]]$importance[2,2]*100,digits=2),"%)"))
  legend("bottomright", legend = unique(fishinfo$Sex), pch=21,col=pal[as.factor(unique(fishinfo$Sex))], pt.bg=adjustcolor(pal[as.factor(unique(fishinfo$Sex))],alpha.f = 0.6))
  dev.off()
}


pal <- RColorBrewer::brewer.pal(n=5,name="Set1")
{png(filename="report_pca3_sex.png",width=6,height=6,units="in",res=300)
  plot(pca[[1]]$PC2,pca[[1]]$PC3,pch=21,col=pal[as.factor(fishinfo$Sex)],bg=adjustcolor(pal[as.factor(fishinfo$Sex)],alpha.f = 0.6),xlab=paste0("PC1 (",round(pca[[2]]$importance[2,2]*100,digits=2),"%)"),ylab=paste0("PC2 (",round(pca[[2]]$importance[2,3]*100,digits=2),"%)"))
  legend("bottomright", legend = unique(fishinfo$Sex)[-3], pch=21,col=pal[as.factor(unique(fishinfo$Sex)[-3])], pt.bg=adjustcolor(pal[as.factor(unique(fishinfo$Sex)[-3])],alpha.f = 0.6))
  dev.off()
}


#pca of spawning bighorn sauger
meta_filt2 <- meta_filt[which(fishinfo$Temporal.Event == "Spawn"),]
file2 <- file[,which(fishinfo$Temporal.Event == "Spawn")]
names2 <- names[which(fishinfo$Temporal.Event == "Spawn")]
colnames(meta_filt2)[4] <- "Tributary"
spawn_pca <- doPCA(file2,meta_filt2)
plotPCApng(spawn_pca,"bighorn_spawn_extra_pcs_sex",meta_filt2$Sex)

pal <- RColorBrewer::brewer.pal(n=5,name="Set1")
{png(filename="bighorn_spawn_pca.png",width=6,height=6,units="in",res=300)
plot(spawn_pca[[1]]$PC1,spawn_pca[[1]]$PC2,pch=21,col=pal[as.factor(meta_filt2$Tributary)],bg=adjustcolor(pal[as.factor(meta_filt2$Tributary)],alpha.f = 0.6),xlab=paste0("PC1 (",round(spawn_pca[[2]]$importance[2,1]*100,digits=2),"%)"),ylab=paste0("PC2 (",round(spawn_pca[[2]]$importance[2,2]*100,digits=2),"%)"))
legend("bottomleft", legend = c("Middle Bighorn","Upper Bighorn","Lower Bighorn"), pch=21,col=pal[as.factor(unique(meta_filt2$Tributary))], pt.bg=adjustcolor(pal[as.factor(unique(meta_filt2$Tributary))],alpha.f = 0.6))
dev.off()
}

#pca of all bighorn sauger
meta_filt3 <- meta_filt[which(fishinfo$Tributary == "Bighorn River"),]
file3 <- file[,which(fishinfo$Tributary == "Bighorn River")]
colnames(meta_filt3)[4] <- "Tributary"
names3 <- names[which(fishinfo$Tributary == "Bighorn River")]
bh_pca <- doPCA(file3,meta_filt3)
plotPCApng(bh_pca,"more_pcs_bh_sex",meta_filt3$Sex)

pal <- RColorBrewer::brewer.pal(n=5,name="Set1")
{png(filename="bighorn_pca.png",width=6,height=6,units="in",res=300)
  plot(bh_pca[[1]]$PC1,bh_pca[[1]]$PC2,pch=21,col=pal[as.factor(meta_filt3$Tributary)],bg=adjustcolor(pal[as.factor(meta_filt3$Tributary)],alpha.f = 0.6),xlab=paste0("PC1 (",round(bh_pca[[2]]$importance[2,1]*100,digits=2),"%)"),ylab=paste0("PC2 (",round(bh_pca[[2]]$importance[2,2]*100,digits=2),"%)"))
  legend("topleft", legend = c("Middle Bighorn","Upper Bighorn","Lower Bighorn","Bighorn Lake"), pch=21,col=pal[as.factor(unique(meta_filt3$Tributary))], pt.bg=adjustcolor(pal[as.factor(unique(meta_filt3$Tributary))],alpha.f = 0.6))
  dev.off()
}

pal <- RColorBrewer::brewer.pal(n=5,name="Set1")
{png(filename="bighorn_pca_sex.png",width=6,height=6,units="in",res=300)
  plot(bh_pca[[1]]$PC1,bh_pca[[1]]$PC2,pch=21,col=pal[as.factor(meta_filt3$Sex)],bg=adjustcolor(pal[as.factor(meta_filt3$Sex)],alpha.f = 0.6),xlab=paste0("PC1 (",round(bh_pca[[2]]$importance[2,1]*100,digits=2),"%)"),ylab=paste0("PC2 (",round(bh_pca[[2]]$importance[2,2]*100,digits=2),"%)"))
  legend("topleft", legend = unique(meta_filt3$Sex)[-3], pch=21,col=pal[as.factor(unique(meta_filt3$Sex)[-3])], pt.bg=adjustcolor(pal[as.factor(unique(meta_filt3$Sex)[-3])],alpha.f = 0.6))
  dev.off()
}


library(RColorBrewer)
points <- c(21,22,23,24)
points_vec <- as.integer(as.factor(meta_filt3$Tributary))+20
points_vec[which(points_vec == 21)] <- 24
meta_filt3$Temporal.Event <- gsub("Home","Non-spawning",meta_filt3$Temporal.Event)
#pal <- list(brewer.pal(n=11,"RdYlGn")[c(1,2,3,4,9,10,11)],brewer.pal(n=9,"YlGnBu")[c(4,6,7:9)],brewer.pal(n=9,"RdPu")[c(6,8)])
pal <- brewer.pal(n=11,"Spectral")[c(2,3,9,10)]
date_fac <- as.numeric(as.factor(meta_filt3$Temporal.Event))

{png(filename="bighorn_pca_TemporalEvent.png",width=6,height=6,units="in",res=300)
  par(mfrow=c(2,2),mar=c(2,2,3,2),oma=c(4,4,0,0))
  for (i in 1:length(unique(meta_filt3$Tributary)[-4])){
    plot(bh_pca[[1]]$PC1,bh_pca[[1]]$PC2,type="n",xlab="",ylab="")
    loc <- unique(meta_filt3$Tributary)[i]
    ind <- which(meta_filt3$Tributary == loc)
    points(bh_pca[[1]]$PC1[which(meta_filt3$Tributary != loc)],bh_pca[[1]]$PC2[which(meta_filt3$Tributary != loc)],pch=points_vec[which(meta_filt3$Tributary != loc)],col="gray80",bg=adjustcolor("gray80",alpha.f = 0.4))
    points(bh_pca[[1]]$PC1[which(meta_filt3$Tributary == loc)],bh_pca[[1]]$PC2[which(meta_filt3$Tributary == loc)],pch=points_vec[which(meta_filt3$Tributary == loc)],col=pal[date_fac[ind]],bg=adjustcolor(pal[date_fac[ind]],alpha.f = 0.8))
  }
  plot(1,1,type="n",axes=F,xlab="",ylab="")
  legend("center", legend = c(unique(meta_filt3$Temporal.Event)[-3],"Middle Bighorn","Upper Bighorn","Lower Bighorn"), pch=c(rep(21,3),23,24,22),col=c(pal[unique(date_fac)[-3]],rep("gray40",3)), pt.bg=adjustcolor(c(pal[unique(date_fac)[-3]],rep("gray40",3)),alpha.f = 0.8),bty="n",ncol=2,cex=0.9,pt.cex = 1,text.width = 0.32,x.intersp = 0.8)
  mtext(paste0("PC1 (",round(bh_pca[[2]]$importance[2,1]*100,digits=2),"%)"), side = 1, outer=TRUE, line=2,adj=0.5,cex = 1)
  mtext(paste0("PC2 (",round(bh_pca[[2]]$importance[2,2]*100,digits=2),"%)"),side=2,outer=TRUE,line=2,adj=0.5,cex=1)
  dev.off()
}





#pca of wr basin sauger
meta_filt4 <- meta_filt[which(fishinfo$Tributary != "Bighorn River"),]
file4 <- file[,which(fishinfo$Tributary != "Bighorn River")]
colnames(meta_filt4)[4] <- "Tributary"
names4 <- names[which(fishinfo$Tributary != "Bighorn River")]
wr_pca <- doPCA(file4,meta_filt4)
plotPCApng(bh_pca,"more_wr_basin",meta_filt4$Tributary)

meta_filt4$Tributary[which(meta_filt4$Tributary == "Boysen")] <- "Boysen Reservoir"

pal <- RColorBrewer::brewer.pal(n=5,name="Set1")
{png(filename="windriver_pca.png",width=6,height=6,units="in",res=300)
  #par(mfrow(c(1,2)))
  plot(wr_pca[[1]]$PC1,wr_pca[[1]]$PC2,pch=21,col=pal[as.factor(meta_filt4$Tributary)],bg=adjustcolor(pal[as.factor(meta_filt4$Tributary)],alpha.f = 0.6),xlab=paste0("PC1 (",round(wr_pca[[2]]$importance[2,1]*100,digits=2),"%)"),ylab=paste0("PC2 (",round(wr_pca[[2]]$importance[2,2]*100,digits=2),"%)"))
  legend("bottomleft", legend = as.factor(unique(meta_filt4$Tributary)), pch=21,col=pal[as.factor(unique(meta_filt4$Tributary))], pt.bg=adjustcolor(pal[as.factor(unique(meta_filt4$Tributary))],alpha.f = 0.6))
  dev.off()
}


#pca of all bighorn sauger, colored by size
meta_filt3 <- meta_filt[which(fishinfo$Tributary == "Bighorn River"),]
file3 <- file[,which(fishinfo$Tributary == "Bighorn River")]
colnames(meta_filt3)[4] <- "Tributary"
names3 <- names[which(fishinfo$Tributary == "Bighorn River")]
bh_pca <- doPCA(file3,meta_filt3)
plotPCApng(bh_pca,"more_pcs_bh_size",meta_filt3$Length)


order <- findInterval(meta_filt3$Length,sort(meta_filt3$Length))
cols <- RColorBrewer::brewer.pal(n=10,"RdYlGn")
pal <- colorRampPalette(cols)


# pal <- RColorBrewer::brewer.pal(n=5,name="Set1")
{png(filename="bighorn_pca_size.png",width=6,height=6,units="in",res=300)
  plot(bh_pca[[1]]$PC1,bh_pca[[1]]$PC2,pch=21,col=pal(nrow(meta_filt3))[order],bg=adjustcolor(pal(nrow(meta_filt3))[order],alpha.f = 0.6),xlab=paste0("PC1 (",round(bh_pca[[2]]$importance[2,1]*100,digits=2),"%)"),ylab=paste0("PC2 (",round(bh_pca[[2]]$importance[2,2]*100,digits=2),"%)"))
  legend("bottomleft", legend = c("8 cm","30 cm"), pch=21,col=pal(2), pt.bg=adjustcolor(pal(2),alpha.f = 0.6))
  dev.off()
}


#PCA of bighorn sauger by location and date
meta_filt3 <- meta_filt[which(fishinfo$Tributary == "Bighorn River"),]
file3 <- file[,which(fishinfo$Tributary == "Bighorn River")]
colnames(meta_filt3)[4] <- "Tributary"
names3 <- names[which(fishinfo$Tributary == "Bighorn River")]
bh_pca <- doPCA(file3,meta_filt3)
#plotPCApng(bh_pca,"more_pcs_bh_sex",meta_filt3$Sex)

library(RColorBrewer)
points <- c(21,22,23,24)
points_vec <- as.integer(as.factor(meta_filt3$Tributary))+20
pal <- list(brewer.pal(n=9,"BuGn"),brewer.pal(n=9,"Blues")[c(6:11)],brewer.pal(n=9,"RdPu")[c(6:8)],brewer.pal(n=9,"Reds")[8])


{png(filename="bighorn_pca_date.png",width=6,height=6,units="in",res=300)
  par(mfrow=c(2,2))
  for (i in 1:length(unique(meta_filt3$Tributary))){
    plot(bh_pca[[1]]$PC1,bh_pca[[1]]$PC2,xlab=paste0("PC1 (",round(bh_pca[[2]]$importance[2,1]*100,digits=2),"%)"),ylab=paste0("PC2 (",round(bh_pca[[2]]$importance[2,2]*100,digits=2),"%)"),type="n")
    loc <- unique(meta_filt3$Tributary)[i]
    ind <- which(meta_filt3$Tributary == loc)
    points(bh_pca[[1]]$PC1[which(meta_filt3$Tributary != loc)],bh_pca[[1]]$PC2[which(meta_filt3$Tributary != loc)],pch=points_vec[which(meta_filt3$Tributary != loc)],col="gray90",bg=adjustcolor("gray90",alpha.f = 0.4))
    points(bh_pca[[1]]$PC1[which(meta_filt3$Tributary == loc)],bh_pca[[1]]$PC2[which(meta_filt3$Tributary == loc)],pch=points_vec[which(meta_filt3$Tributary == loc)],col=pal[[i]][as.factor(meta_filt3$Date[ind])],bg=adjustcolor(pal[[i]][as.factor(meta_filt3$Date[ind])],alpha.f = 0.6))
  }
  #legend("bottomleft", legend = c("Middle Bighorn","Upper Bighorn","Lower Bighorn"), pch=21,col=pal[as.factor(unique(meta_filt2$Tributary))], pt.bg=adjustcolor(pal[as.factor(unique(meta_filt2$Tributary))],alpha.f = 0.6))
  dev.off()
}


######################################################################
#PCA for each location, with temporal event colored
meta_filt5 <- meta_filt3[which(meta_filt3$Tributary == "Bighorn_Middle"),]
file5 <- file3[,which(meta_filt3$Tributary == "Bighorn_Middle")]
mbh_pca <- doPCA(file5,meta_filt5)

pal <- RColorBrewer::brewer.pal(n=11,name="Spectral")[c(3,9,10)]
{png(filename="midbighorn_temporalevent_pca.png",width=6,height=6,units="in",res=300)
  plot(mbh_pca[[1]]$PC1,mbh_pca[[1]]$PC2,pch=21,col=pal[as.factor(meta_filt5$Temporal.Event)],bg=adjustcolor(pal[as.factor(meta_filt5$Temporal.Event)],alpha.f = 0.6),xlab=paste0("PC1 (",round(mbh_pca[[2]]$importance[2,1]*100,digits=2),"%)"),ylab=paste0("PC2 (",round(mbh_pca[[2]]$importance[2,2]*100,digits=2),"%)"))
  legend("topleft", legend = unique(meta_filt5$Temporal.Event), pch=21,col=pal[as.factor(unique(meta_filt5$Temporal.Event))], pt.bg=adjustcolor(pal[as.factor(unique(meta_filt5$Temporal.Event))],alpha.f = 0.6))
  dev.off()
}

meta_filt6 <- meta_filt3[which(meta_filt3$Tributary == "Bighorn_Upper"),]
file6 <- file3[,which(meta_filt3$Tributary == "Bighorn_Upper")]
ubh_pca <- doPCA(file6,meta_filt6)

pal <- RColorBrewer::brewer.pal(n=11,name="Spectral")[c(3,9,10)]
{png(filename="upbighorn_temporalevent_pca.png",width=6,height=6,units="in",res=300)
  plot(ubh_pca[[1]]$PC1,ubh_pca[[1]]$PC2,pch=21,col=pal[as.factor(meta_filt6$Temporal.Event)],bg=adjustcolor(pal[as.factor(meta_filt6$Temporal.Event)],alpha.f = 0.6),xlab=paste0("PC1 (",round(ubh_pca[[2]]$importance[2,1]*100,digits=2),"%)"),ylab=paste0("PC2 (",round(ubh_pca[[2]]$importance[2,2]*100,digits=2),"%)"))
  legend("topleft", legend = unique(meta_filt6$Temporal.Event), pch=21,col=pal[as.factor(unique(meta_filt6$Temporal.Event))], pt.bg=adjustcolor(pal[as.factor(unique(meta_filt6$Temporal.Event))],alpha.f = 0.6))
  dev.off()
}

meta_filt7 <- meta_filt3[which(meta_filt3$Tributary == "Bighorn_Lower"),]
file7 <- file3[,which(meta_filt3$Tributary == "Bighorn_Lower")]
lbh_pca <- doPCA(file7,meta_filt7)

pal <- RColorBrewer::brewer.pal(n=11,name="Spectral")[c(3,9,10)]
{png(filename="lowbighorn_temporalevent_pca.png",width=6,height=6,units="in",res=300)
  plot(lbh_pca[[1]]$PC1,lbh_pca[[1]]$PC2,pch=21,col=pal[as.factor(meta_filt7$Temporal.Event)],bg=adjustcolor(pal[as.factor(meta_filt7$Temporal.Event)],alpha.f = 0.6),xlab=paste0("PC1 (",round(lbh_pca[[2]]$importance[2,1]*100,digits=2),"%)"),ylab=paste0("PC2 (",round(lbh_pca[[2]]$importance[2,2]*100,digits=2),"%)"))
  legend("topleft", legend = unique(meta_filt7$Temporal.Event), pch=21,col=pal[as.factor(unique(meta_filt7$Temporal.Event))], pt.bg=adjustcolor(pal[as.factor(unique(meta_filt7$Temporal.Event))],alpha.f = 0.6))
  dev.off()
}

#meta_filt5$Temporal.Event <- factor(meta_filt5$Temporal.Event,levels=c("Pre-Spawn ","Spawning","Non-spawning"))
#meta_filt6$Temporal.Event <- factor(meta_filt6$Temporal.Event,levels=c("Pre-Spawn ","Spawning","Non-spawning"))
#meta_filt7$Temporal.Event <- factor(meta_filt7$Temporal.Event,levels=c("Pre-Spawn ","Spawning","Non-spawning"))

#turn into one big plot
legend_names <- c("Pre-Spawn","Spawning","Non-spawning")
pal <- RColorBrewer::brewer.pal(n=9,name="Set1")[c(4,5,3)]
{png(filename="all_location_tempevent_pca.png",width=6,height=6,units="in",res=300)
  par(mfrow=c(2,2))
  plot(mbh_pca[[1]]$PC1,mbh_pca[[1]]$PC2,pch=21,col=pal[as.factor(meta_filt5$Temporal.Event)],main="Middle Bighorn sauger",bg=adjustcolor(pal[as.factor(meta_filt5$Temporal.Event)],alpha.f = 0.6),xlab=paste0("PC1 (",round(mbh_pca[[2]]$importance[2,1]*100,digits=2),"%)"),ylab=paste0("PC2 (",round(mbh_pca[[2]]$importance[2,2]*100,digits=2),"%)"))
  plot(ubh_pca[[1]]$PC1,ubh_pca[[1]]$PC2,pch=21,col=pal[as.factor(meta_filt6$Temporal.Event)],main="Upper Bighorn sauger",bg=adjustcolor(pal[as.factor(meta_filt6$Temporal.Event)],alpha.f = 0.6),xlab=paste0("PC1 (",round(ubh_pca[[2]]$importance[2,1]*100,digits=2),"%)"),ylab=paste0("PC2 (",round(ubh_pca[[2]]$importance[2,2]*100,digits=2),"%)"))
  plot(lbh_pca[[1]]$PC1,lbh_pca[[1]]$PC2,pch=21,col=pal[as.factor(meta_filt7$Temporal.Event)],main="Lower Bighorn sauger",bg=adjustcolor(pal[as.factor(meta_filt7$Temporal.Event)],alpha.f = 0.6),xlab=paste0("PC1 (",round(lbh_pca[[2]]$importance[2,1]*100,digits=2),"%)"),ylab=paste0("PC2 (",round(lbh_pca[[2]]$importance[2,2]*100,digits=2),"%)"))
  plot(1,1,xlab="",ylab="",type="n",axes=F)
  legend("center", bty="n",legend = legend_names, pch=21,col=pal, pt.bg=adjustcolor(pal,alpha.f = 0.6))
  dev.off()
}


#now with PC2 & 3
{png(filename="all_location_tempevent_2vs3_pca.png",width=6,height=6,units="in",res=300)
  par(mfrow=c(2,2))
  plot(mbh_pca[[1]]$PC2,mbh_pca[[1]]$PC3,pch=21,col=pal[as.factor(meta_filt5$Temporal.Event)],main="Middle Bighorn sauger",bg=adjustcolor(pal[as.factor(meta_filt5$Temporal.Event)],alpha.f = 0.6),xlab=paste0("PC2 (",round(mbh_pca[[2]]$importance[2,2]*100,digits=2),"%)"),ylab=paste0("PC3 (",round(mbh_pca[[2]]$importance[2,3]*100,digits=2),"%)"))
  plot(ubh_pca[[1]]$PC2,ubh_pca[[1]]$PC3,pch=21,col=pal[as.factor(meta_filt6$Temporal.Event)],main="Upper Bighorn sauger",bg=adjustcolor(pal[as.factor(meta_filt6$Temporal.Event)],alpha.f = 0.6),xlab=paste0("PC2 (",round(ubh_pca[[2]]$importance[2,2]*100,digits=2),"%)"),ylab=paste0("PC3 (",round(ubh_pca[[2]]$importance[2,3]*100,digits=2),"%)"))
  plot(lbh_pca[[1]]$PC2,lbh_pca[[1]]$PC3,pch=21,col=pal[as.factor(meta_filt7$Temporal.Event)],main="Lower Bighorn sauger",bg=adjustcolor(pal[as.factor(meta_filt7$Temporal.Event)],alpha.f = 0.6),xlab=paste0("PC2 (",round(lbh_pca[[2]]$importance[2,2]*100,digits=2),"%)"),ylab=paste0("PC3 (",round(lbh_pca[[2]]$importance[2,3]*100,digits=2),"%)"))
  plot(1,1,xlab="",ylab="",type="n",axes=F)
  legend("center", bty="n",legend = legend_names, pch=21,col=pal, pt.bg=adjustcolor(pal,alpha.f = 0.6))
  dev.off()
}


######################################################################
#PCA for each temporal event with location colored
meta_filt8 <- meta_filt3[which(meta_filt3$Temporal.Event == "Pre-Spawn "),]
file8 <- file3[,which(meta_filt3$Temporal.Event == "Pre-Spawn ")]
prespawn_pca <- doPCA(file8,meta_filt8)

pal <- RColorBrewer::brewer.pal(n=5,name="Set1")[-1]
{png(filename="prespawn_pca.png",width=6,height=6,units="in",res=300)
  plot(prespawn_pca[[1]]$PC1,prespawn_pca[[1]]$PC2,pch=21,col=pal[as.factor(meta_filt8$Tributary)],main="Pre-spawn sauger",bg=adjustcolor(pal[as.factor(meta_filt8$Tributary)],alpha.f = 0.6),xlab=paste0("PC1 (",round(prespawn_pca[[2]]$importance[2,1]*100,digits=2),"%)"),ylab=paste0("PC2 (",round(prespawn_pca[[2]]$importance[2,2]*100,digits=2),"%)"))
  legend("topright", legend = c("Middle Bighorn","Lower Bighorn","Upper Bighorn"), pch=21,col=pal[as.factor(unique(meta_filt8$Tributary))], pt.bg=adjustcolor(pal[as.factor(unique(meta_filt8$Tributary))],alpha.f = 0.6))
  dev.off()
}


meta_filt9 <- meta_filt3[which(meta_filt3$Temporal.Event == "Spawn"),]
file9 <- file3[,which(meta_filt3$Temporal.Event == "Spawn")]
spawn_pca <- doPCA(file9,meta_filt9)

pal <- RColorBrewer::brewer.pal(n=5,name="Set1")[-1]
{png(filename="spawn_pca.png",width=6,height=6,units="in",res=300)
  plot(spawn_pca[[1]]$PC1,spawn_pca[[1]]$PC2,pch=21,col=pal[as.factor(meta_filt9$Tributary)],main="Spawning sauger",bg=adjustcolor(pal[as.factor(meta_filt9$Tributary)],alpha.f = 0.6),xlab=paste0("PC1 (",round(spawn_pca[[2]]$importance[2,1]*100,digits=2),"%)"),ylab=paste0("PC2 (",round(spawn_pca[[2]]$importance[2,2]*100,digits=2),"%)"))
  legend("bottomleft", legend = c("Middle Bighorn","Lower Bighorn","Upper Bighorn"), pch=21,col=pal[as.factor(unique(meta_filt9$Tributary))], pt.bg=adjustcolor(pal[as.factor(unique(meta_filt9$Tributary))],alpha.f = 0.6))
  dev.off()
}

meta_filt10 <- meta_filt3[which(meta_filt3$Temporal.Event == "Non-spawning"),]
file10 <- file3[,which(meta_filt3$Temporal.Event == "Non-spawning")]
home_pca <- doPCA(file10,meta_filt10)

pal <- RColorBrewer::brewer.pal(n=5,name="Set1")[-1]
{png(filename="home_pca.png",width=6,height=6,units="in",res=300)
  plot(home_pca[[1]]$PC1,home_pca[[1]]$PC2,pch=21,col=pal[as.factor(meta_filt10$Tributary)],main="Home-range sauger",bg=adjustcolor(pal[as.factor(meta_filt10$Tributary)],alpha.f = 0.6),xlab=paste0("PC1 (",round(home_pca[[2]]$importance[2,1]*100,digits=2),"%)"),ylab=paste0("PC2 (",round(home_pca[[2]]$importance[2,2]*100,digits=2),"%)"))
  legend("bottomleft", legend = c("Middle Bighorn","Lower Bighorn","Upper Bighorn"), pch=21,col=pal[as.factor(unique(meta_filt10$Tributary))], pt.bg=adjustcolor(pal[as.factor(unique(meta_filt10$Tributary))],alpha.f = 0.6))
  dev.off()
}

#turn into one big plot

pal <- RColorBrewer::brewer.pal(n=5,name="Set1")
{png(filename="all_temporalevent_pca.png",width=6,height=6,units="in",res=300)
  par(mfrow=c(2,2))
  plot(prespawn_pca[[1]]$PC1,prespawn_pca[[1]]$PC2,pch=21,col=pal[as.factor(meta_filt8$Tributary)],main="Pre-spawn sauger",bg=adjustcolor(pal[as.factor(meta_filt8$Tributary)],alpha.f = 0.6),xlab=paste0("PC1 (",round(prespawn_pca[[2]]$importance[2,1]*100,digits=2),"%)"),ylab=paste0("PC2 (",round(prespawn_pca[[2]]$importance[2,2]*100,digits=2),"%)"))
  plot(spawn_pca[[1]]$PC1,spawn_pca[[1]]$PC2,pch=21,col=pal[as.factor(meta_filt9$Tributary)],main="Spawning sauger",bg=adjustcolor(pal[as.factor(meta_filt9$Tributary)],alpha.f = 0.6),xlab=paste0("PC1 (",round(spawn_pca[[2]]$importance[2,1]*100,digits=2),"%)"),ylab=paste0("PC2 (",round(spawn_pca[[2]]$importance[2,2]*100,digits=2),"%)"))
  plot(home_pca[[1]]$PC1,home_pca[[1]]$PC2,pch=21,col=pal[as.factor(meta_filt10$Tributary)],main="Non-spawning sauger",bg=adjustcolor(pal[as.factor(meta_filt10$Tributary)],alpha.f = 0.6),xlab=paste0("PC1 (",round(home_pca[[2]]$importance[2,1]*100,digits=2),"%)"),ylab=paste0("PC2 (",round(home_pca[[2]]$importance[2,2]*100,digits=2),"%)"))
  plot(1,1,xlab="",ylab="",type="n",axes=F)
  legend("center", bty="n",legend = c("Middle Bighorn","Lower Bighorn","Upper Bighorn"), pch=21,col=pal[as.factor(unique(meta_filt8$Tributary))], pt.bg=adjustcolor(pal[as.factor(unique(meta_filt8$Tributary))],alpha.f = 0.6))
  dev.off()
}


#####################################################################
#PCA of spawning bighorn sauger by location and date
meta_filt2 <- meta_filt[which(fishinfo$Temporal.Event == "Spawn"),]
file2 <- file[,which(fishinfo$Temporal.Event == "Spawn")]
names2 <- names[which(fishinfo$Temporal.Event == "Spawn")]
colnames(meta_filt2)[4] <- "Tributary"
spawn_pca <- doPCA(file2,meta_filt2)
#plotPCApng(bh_pca,"more_pcs_bh_sex",meta_filt3$Sex)

library(RColorBrewer)
points <- c(21,22,23,24)
points_vec <- as.integer(as.factor(meta_filt2$Tributary))+20
points_vec[which(points_vec == 21)] <- 24
#pal <- list(brewer.pal(n=11,"RdYlGn")[c(1,2,3,4,9,10,11)],brewer.pal(n=9,"YlGnBu")[c(4,6,7:9)],brewer.pal(n=9,"RdPu")[c(6,8)])
pal <- brewer.pal(n=11,"Spectral")[c(1:4,9:11)]
date_fac <- as.numeric(as.factor(meta_filt2$Date))

{png(filename="bighorn_spawn_pca_date.png",width=6,height=6,units="in",res=300)
  par(mfrow=c(2,2),mar=c(2,2,3,2),oma=c(4,4,0,0))
  for (i in 1:length(unique(meta_filt2$Tributary))){
    plot(spawn_pca[[1]]$PC1,spawn_pca[[1]]$PC2,type="n",xlab="",ylab="")
    loc <- unique(meta_filt2$Tributary)[i]
    ind <- which(meta_filt2$Tributary == loc)
    points(spawn_pca[[1]]$PC1[which(meta_filt2$Tributary != loc)],spawn_pca[[1]]$PC2[which(meta_filt2$Tributary != loc)],pch=points_vec[which(meta_filt2$Tributary != loc)],col="gray80",bg=adjustcolor("gray80",alpha.f = 0.4))
    points(spawn_pca[[1]]$PC1[which(meta_filt2$Tributary == loc)],spawn_pca[[1]]$PC2[which(meta_filt2$Tributary == loc)],pch=points_vec[which(meta_filt2$Tributary == loc)],col=pal[date_fac[ind]],bg=adjustcolor(pal[date_fac[ind]],alpha.f = 0.8))
  }
  plot(1,1,type="n",axes=F,xlab="",ylab="")
  legend("center", legend = c(unique(meta_filt2$Date),"Middle Bighorn","Upper Bighorn","Lower Bighorn",rep("",4)), pch=c(rep(21,7),22,23,24,rep(NA,4)),col=c(pal[unique(date_fac)],rep("gray40",7)), pt.bg=adjustcolor(c(pal[unique(date_fac)],rep("gray40",7)),alpha.f = 0.8),bty="n",ncol=2,cex=0.9,pt.cex = 1,text.width = 0.25,x.intersp = 0.8)
  mtext(paste0("PC1 (",round(spawn_pca[[2]]$importance[2,1]*100,digits=2),"%)"), side = 1, outer=TRUE, line=2,adj=0.5,cex = 1)
  mtext(paste0("PC2 (",round(spawn_pca[[2]]$importance[2,2]*100,digits=2),"%)"),side=2,outer=TRUE,line=2,adj=0.5,cex=1)
  dev.off()
}


#####################################################################
#PCA of bighorn sauger by date and maturity
meta_filt3 <- fishinfo[which(fishinfo$Tributary == "Bighorn River"),]
file3 <- file[,which(fishinfo$Tributary == "Bighorn River")]
colnames(meta_filt3)[4] <- "Tributary"
names3 <- names[which(fishinfo$Tributary == "Bighorn River")]
bh_pca <- doPCA(file3,meta_filt3)

meta_filt3$Maturity[grep("Inter. Firm",meta_filt3$Maturity)] <- "Firm"
meta_filt3$Maturity[grep("Inter. Soft",meta_filt3$Maturity)] <- "Soft"
meta_filt3$Maturity[grep("Intermediate",meta_filt3$Maturity)] <- "Soft"
meta_filt3$Maturity[grep("Ripe",meta_filt3$Maturity)] <- "Soft"
meta_filt3$Maturity[grep("Inter. Firm",meta_filt3$Maturity)] <- "Firm"
meta_filt3$Maturity[which(meta_filt3$Maturity == "")] <- NA

meta_filt3 <- meta_filt3[grep("2018",meta_filt3$Date),]

library(RColorBrewer)
points <- c(22,23,24,25)
points_vec <- as.integer(as.factor(meta_filt3$Maturity))+21
#points_vec[which(points_vec == 21)] <- 24
pal <- brewer.pal(n=11,"Spectral")[c(1:4,8:11)]
meta_filt3$Date <- as.Date(meta_filt3$Date,format="%m/%d/%Y")
meta_filt3$Date <- factor(meta_filt3$Date)
date_fac <- as.numeric(meta_filt3$Date)

{png(filename="bighorn_maturity_pca_date.png",width=6,height=9,units="in",res=300)
  par(mfcol=c(3,2),mar=c(2,2,3,2),oma=c(4,4,0,0))
  for (i in 2:length(unique(meta_filt3$Maturity))){
    plot(bh_pca[[1]]$PC1,bh_pca[[1]]$PC2,type="n",xlab="",ylab="")
    loc <- unique(meta_filt3$Maturity)[i]
    ind <- which(meta_filt3$Maturity == loc)
    points(bh_pca[[1]]$PC1[which(meta_filt3$Maturity != loc)],bh_pca[[1]]$PC2[which(meta_filt3$Maturity != loc)],pch=points_vec[which(meta_filt3$Maturity != loc)],col="gray80",bg=adjustcolor("gray80",alpha.f = 0.4),cex=1.3)
    points(bh_pca[[1]]$PC1[which(meta_filt3$Maturity == loc)],bh_pca[[1]]$PC2[which(meta_filt3$Maturity == loc)],pch=points_vec[which(meta_filt3$Maturity == loc)],col=pal[date_fac[ind]],bg=adjustcolor(pal[date_fac[ind]],alpha.f = 0.8),cex=1.3)
  }
  plot(1,1,type="n",axes=F,xlab="",ylab="")
  legend("center", legend = c(levels(meta_filt3$Date),"Green","Soft","Spent","Firm",rep("",3)), 
         pch=c(rep(21,8),22,23,24,25,rep(NA,4)),col=c(pal[unique(sort(date_fac))],rep("gray40",14)), 
         pt.bg=adjustcolor(c(pal[unique(sort(date_fac))],rep("gray40",7)),alpha.f = 0.8),bty="n",
         ncol=2,cex=0.9,pt.cex = 1,text.width = 0.25,x.intersp = 0.8)
  mtext(paste0("PC1 (",round(bh_pca[[2]]$importance[2,1]*100,digits=2),"%)"), side = 1, outer=TRUE, line=2,adj=0.5,cex = 1)
  mtext(paste0("PC2 (",round(bh_pca[[2]]$importance[2,2]*100,digits=2),"%)"),side=2,outer=TRUE,line=2,adj=0.5,cex=1)
  dev.off()
}




#####################################################################
#PCA of spawning bighorn sauger by location river mile
meta_filt3 <- meta_filt[which(fishinfo$Tributary == "Bighorn River"),]
file3 <- file[,which(fishinfo$Tributary == "Bighorn River")]
colnames(meta_filt3)[4] <- "Tributary"
names3 <- names[which(fishinfo$Tributary == "Bighorn River")]
bh_pca <- doPCA(file3,meta_filt3)
#plotPCApng(bh_pca,"more_pcs_bh_sex",meta_filt3$Sex)

#meta_filt3 <- meta_filt3[-which(meta_filt3$Tributary == "Bighorn_Lake"),]
#meta_filt_bhl <- meta_filt3[which(meta_filt3$Tributary == "Bighorn_Lake"),]

library(RColorBrewer)
library(wesanderson)
points <- c(21,22,23,24)
points_vec <- as.integer(as.factor(meta_filt3$Tributary))+20
points_vec[which(points_vec == 21)] <- 24


pal <- list(list(brewer.pal(n=11,"RdYlGn")[c(1:3,9:11)],brewer.pal(n=11,"PuOr")[c(3,10)],brewer.pal(n=11,"RdBu")[2]),
            list(brewer.pal(n=11,"RdYlGn")[-6],brewer.pal(n=11,"PuOr")[c(1:3,8:11)],brewer.pal(n=11,"RdBu")[c(1:3,8:11)]))

years <- c("2016","2018")

{for (y in 1:length(years)){
  
  meta_year <- meta_filt3[grep(years[y],meta_filt3$Date),]
  year_ind <- grep(years[y],meta_filt3$Date)
  filename <- paste0("./bighorn_pca_rivermile_",years[y],".png")
  png(filename=filename,width=6,height=9.5,units="in",res=300)
  par(mfrow=c(3,1),mar=c(2,2,3,2),oma=c(4,4,0,0))
  for (i in 1:length(unique(meta_year$Tributary)[-4])){
    plot(bh_pca[[1]]$PC1,bh_pca[[1]]$PC2,type="n",xlab="",ylab="")
    pca_year <- bh_pca[[1]][year_ind,3:4]
    loc <- unique(meta_year$Tributary)[i]
    ind <- which(meta_year$Tributary == loc)
    ind2 <- intersect(which(meta_filt3$Tributary == loc),year_ind)
    met <- meta_year[ind,]
    
    river_mile_clean <- met$River_Mile
    river_mile_clean <- as.vector(sapply(river_mile_clean,function(x) mean(as.numeric(strsplit(x,"-")[[1]]))))
    order <- order(unique(river_mile_clean))
    met$River_Mile <- factor(met$River_Mile,levels = unique(met$River_Mile)[order])
    
    
    points(bh_pca[[1]]$PC1[which(meta_filt3$EGM_ID %in% met$EGM_ID == FALSE)],bh_pca[[1]]$PC2[which( meta_filt3$EGM_ID %in% met$EGM_ID == FALSE)],cex=1.3,pch=points_vec[which(meta_filt3$Tributary != loc| meta_filt3$EGM_ID %in% met$EGM_ID == FALSE)],col="gray80",bg=adjustcolor("gray80",alpha.f = 0.4))
    points(pca_year$PC1[ind],pca_year$PC2[ind],pch=points_vec[meta_filt3$EGM_ID %in% met$EGM_ID],cex=1.3,col=pal[[y]][[i]][as.integer(met$River_Mile)],bg=adjustcolor(pal[[y]][[i]][as.integer(met$River_Mile)],alpha.f = 0.8))
    legend("topleft", legend = c(unique(as.vector(met$River_Mile)),"Middle Bighorn","Upper Bighorn","Lower Bighorn"), 
           pch=c(rep(21,length(unique(met$River_Mile))),23,24,22),
           col=c(pal[[y]][[i]][as.integer(unique(met$River_Mile))],rep("gray40",3)), 
           pt.bg=adjustcolor(c(pal[[y]][[i]][as.integer(unique(met$River_Mile))],rep("gray40",3)),alpha.f = 0.8),
           cex=0.9,pt.cex = 1,x.intersp = 0.8,bg=adjustcolor("white",alpha.f = 0.8))
  }
  mtext(paste0("PC1 (",round(bh_pca[[2]]$importance[2,1]*100,digits=2),"%)"), side = 1, outer=TRUE, line=2,adj=0.5,cex = 1)
  mtext(paste0("PC2 (",round(bh_pca[[2]]$importance[2,2]*100,digits=2),"%)"),side=2,outer=TRUE,line=2,adj=0.5,cex=1)
  dev.off()
}
}



####################################################################
#see if divergence on saug PCA 3 is due to library effects

lib1 <- read.csv(file="D:/Masters/sauger/misc/SARWAE1_barcode_key.csv",header=T)
lib2 <- read.csv(file="D:/Masters/sauger/misc/SARWAE2_barcode_key.csv",header=T)
lib3 <- read.csv(file="D:/Masters/sauger/misc/SARWAE3_barcode_key.csv",header=T)
lib4 <- read.csv(file="D:/Masters/sauger/misc/SARWAE4_barcode_key.csv",header=T)

lib_data <- rbind(lib1,lib2,lib3,lib4)
lib_data$lib <- c(rep("1",207),rep("2",192),rep("3",399),rep("4",399))
dups <- as.vector(lib_data$Idcode[which(duplicated(lib_data$Idcode))])
lib_data$lib[which(lib_data$Idcode %in% dups)] <- "3+4"
lib_data <- unique(lib_data)


source("D:/Masters/scripts/pca_funcs.R")

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
fishinfo$lib <- lib_data$lib[match(fishinfo$EGM_ID,lib_data$Idcode)]


pca <- doPCA(file,fishinfo)

plotPCApng(pca,"saug_more_pcs_lib",fishinfo$lib)

pal <- RColorBrewer::brewer.pal(n=5,name="Dark2")
{tiff(filename="report_pca3_lib.tiff",width=6,height=6,units="in",res=300)
  plot(pca[[1]]$PC2,pca[[1]]$PC3,pch=21,col=pal[as.factor(fishinfo$lib)],bg=adjustcolor(pal[as.factor(fishinfo$lib)],alpha.f = 0.6),xlab=paste0("PC2 (",round(pca[[2]]$importance[2,1]*100,digits=2),"%)"),ylab=paste0("PC3 (",round(pca[[2]]$importance[2,2]*100,digits=2),"%)"),sub="6,536 loci")
  legend("bottomright", legend = unique(fishinfo$lib), pch=21,col=pal[as.factor(unique(fishinfo$lib))], pt.bg=adjustcolor(pal[as.factor(unique(fishinfo$lib))],alpha.f = 0.6))
  dev.off()
}






#bighorn
meta_filt3 <- fishinfo[which(fishinfo$Tributary == "Bighorn River"),]
file3 <- file[,which(fishinfo$Tributary == "Bighorn River")]
colnames(meta_filt3)[4] <- "Tributary"
names3 <- names[which(fishinfo$Tributary == "Bighorn River")]
bh_pca <- doPCA(file3,meta_filt3)

pal <- RColorBrewer::brewer.pal(n=5,name="Dark2")
{png(filename="bighorn_pca_lib.png",width=6,height=6,units="in",res=300)
  plot(bh_pca[[1]]$PC1,bh_pca[[1]]$PC2,pch=21,col=pal[as.factor(meta_filt3$lib)],bg=adjustcolor(pal[as.factor(meta_filt3$lib)],alpha.f = 0.6),xlab=paste0("PC1 (",round(bh_pca[[2]]$importance[2,1]*100,digits=2),"%)"),ylab=paste0("PC2 (",round(bh_pca[[2]]$importance[2,2]*100,digits=2),"%)"))
  legend("bottomright", legend = unique(meta_filt3$lib), pch=21,col=pal[as.factor(unique(meta_filt3$lib))], pt.bg=adjustcolor(pal[as.factor(unique(meta_filt3$lib))],alpha.f = 0.6))
  dev.off()
}




#no effect
#double check to make sure library matched correctly
lib_data_filt <- lib_data[which(lib_data$Idcode %in% fishinfo$EGM_ID),]
lib_data_filt <- lib_data_filt[match(fishinfo$EGM_ID,lib_data_filt$Idcode),]
lib_data_filt$Idcode == fishinfo$EGM_ID
lib_data_filt$lib == fishinfo$lib

#maybe there are spawning cohorts?

meta_filt2 <- meta_filt[which(fishinfo$Temporal.Event == "Spawn"),]
meta_filt2$Maturity[grep("Inter. Firm",meta_filt2$Maturity)] <- "Firm"
meta_filt2$Maturity[grep("Inter. Soft",meta_filt2$Maturity)] <- "Soft"
meta_filt2$Maturity[grep("Intermediate",meta_filt2$Maturity)] <- "Soft"
meta_filt2$Maturity[grep("Ripe",meta_filt2$Maturity)] <- "Soft"
meta_filt2$Maturity[grep("Inter. Firm",meta_filt2$Maturity)] <- "Firm"


file2 <- file[,which(fishinfo$Temporal.Event == "Spawn")]
names2 <- names[which(fishinfo$Temporal.Event == "Spawn")]
colnames(meta_filt2)[4] <- "Tributary"
spawn_pca <- doPCA(file2,meta_filt2)
plotPCApng(spawn_pca,"bighorn_spawn_extra_pcs_maturity",meta_filt2$Maturity)

pal <- RColorBrewer::brewer.pal(n=5,name="Set1")
{png(filename="bighorn_spawn_pca.png",width=6,height=6,units="in",res=300)
  plot(spawn_pca[[1]]$PC1,spawn_pca[[1]]$PC2,pch=21,col=pal[as.factor(meta_filt2$Tributary)],bg=adjustcolor(pal[as.factor(meta_filt2$Tributary)],alpha.f = 0.6),xlab=paste0("PC1 (",round(spawn_pca[[2]]$importance[2,1]*100,digits=2),"%)"),ylab=paste0("PC2 (",round(spawn_pca[[2]]$importance[2,2]*100,digits=2),"%)"))
  legend("bottomleft", legend = c("Middle Bighorn","Upper Bighorn","Lower Bighorn"), pch=21,col=pal[as.factor(unique(meta_filt2$Tributary))], pt.bg=adjustcolor(pal[as.factor(unique(meta_filt2$Tributary))],alpha.f = 0.6))
  dev.off()
}


#maybe there are spawning cohorts?

fishinfo$Maturity[grep("Inter. Firm",fishinfo$Maturity)] <- "Firm"
fishinfo$Maturity[grep("Inter. Soft",fishinfo$Maturity)] <- "Soft"
fishinfo$Maturity[grep("Intermediate",fishinfo$Maturity)] <- "Soft"
fishinfo$Maturity[grep("Ripe",fishinfo$Maturity)] <- "Soft"
fishinfo$Maturity[grep("Inter. Firm",fishinfo$Maturity)] <- "Firm"
fishinfo$Maturity[which(fishinfo$Maturity == "")] <- NA

pal <- RColorBrewer::brewer.pal(n=5,name="Dark2")
{png(filename="maturity_pca.png",width=6,height=6,units="in",res=300)
  plot(pca[[1]]$PC2,pca[[1]]$PC3,pch=21,col=pal[as.factor(fishinfo$Maturity)],bg=adjustcolor(pal[as.factor(fishinfo$Maturity)],alpha.f = 0.6),xlab=paste0("PC2 (",round(pca[[2]]$importance[2,2]*100,digits=2),"%)"),ylab=paste0("PC3 (",round(pca[[2]]$importance[2,3]*100,digits=2),"%)"))
  legend("bottomleft", legend = unique(fishinfo$Maturity)[-1], pch=21,col=pal[as.factor(unique(fishinfo$Maturity)[-1])], pt.bg=adjustcolor(pal[as.factor(unique(fishinfo$Maturity)[-1])],alpha.f = 0.6))
  dev.off()
}


meta_filt3 <- meta_filt[which(fishinfo$Tributary == "Bighorn River"),]
file3 <- file[,which(fishinfo$Tributary == "Bighorn River")]
colnames(meta_filt3)[4] <- "Tributary"
names3 <- names[which(fishinfo$Tributary == "Bighorn River")]
bh_pca <- doPCA(file3,meta_filt3)

meta_filt3$Maturity[grep("Inter. Firm",meta_filt3$Maturity)] <- "Firm"
meta_filt3$Maturity[grep("Inter. Soft",meta_filt3$Maturity)] <- "Soft"
meta_filt3$Maturity[grep("Intermediate",meta_filt3$Maturity)] <- "Soft"
meta_filt3$Maturity[grep("Ripe",meta_filt3$Maturity)] <- "Soft"
meta_filt3$Maturity[grep("Inter. Firm",meta_filt3$Maturity)] <- "Firm"
meta_filt3$Maturity[which(meta_filt3$Maturity == "")] <- NA


pal <- RColorBrewer::brewer.pal(n=5,name="Set1")
{png(filename="bighorn_maturity_pca.png",width=6,height=6,units="in",res=300)
  plot(bh_pca[[1]]$PC1,bh_pca[[1]]$PC2,pch=21,col=pal[as.factor(meta_filt3$Maturity)],bg=adjustcolor(pal[as.factor(meta_filt3$Maturity)],alpha.f = 0.6),xlab=paste0("PC1 (",round(bh_pca[[2]]$importance[2,1]*100,digits=2),"%)"),ylab=paste0("PC2 (",round(bh_pca[[2]]$importance[2,2]*100,digits=2),"%)"))
  legend("bottomleft", legend = unique(meta_filt3$Maturity)[-1], pch=21,col=pal[as.factor(unique(meta_filt3$Maturity)[-1])], pt.bg=adjustcolor(pal[as.factor(unique(meta_filt3$Maturity)[-1])],alpha.f = 0.6))
  dev.off()
}





###################################################################
#correlation between bh_pca1 and river mile
#spawning sauger
river_mile_clean <- meta_filt2$River_Mile
river_mile_clean <- as.vector(sapply(river_mile_clean,function(x) mean(as.numeric(strsplit(x,"-")[[1]]))))

spawn_df <- data.frame("rm"=river_mile_clean,"PC1"=spawn_pca[[1]]$PC1,"PC2"=spawn_pca[[1]]$PC2)

#all bighorn sauger
river_mile_clean <- meta_filt3$River_Mile
river_mile_clean <- as.vector(sapply(river_mile_clean,function(x) mean(as.numeric(strsplit(x,"-")[[1]]))))

bh_df <- data.frame("rm"=as.numeric(river_mile_clean),"PC1"=bh_pca[[1]]$PC1,"PC2"=bh_pca[[1]]$PC2)

bh_df$rm <- gsub("NaN",1,bh_df$rm)
bh_df$rm <- as.numeric(bh_df$rm)

#plot results
library(ggpubr)
pal <- RColorBrewer::brewer.pal(n=8,"Dark2")
ggbh <- ggscatter(bh_df,x="rm",y="PC1",color=pal[6],shape=21,fill=adjustcolor(pal[6],alpha.f = 0.6),add="reg.line",conf.int = TRUE,add.params = list(color="red",fill="gray60"),size=4) +
  stat_regline_equation(label.x = 16,label.y=-0.16,aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),size=5,fontface="bold")
ggbh1 <- ggpar(ggbh,xlab="River mile",ylab="PC1 loading")

ggbh2 <- ggscatter(bh_df,x="rm",y="PC2",color=pal[6],shape=21,fill=adjustcolor(pal[6],alpha.f = 0.6),add="reg.line",conf.int = TRUE,add.params = list(color="red",fill="gray60"),size=4) +
  stat_regline_equation(label.x = 16,label.y=-0.12,aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),size=5,fontface="bold")
ggbh3 <- ggpar(ggbh2,xlab="River mile",ylab="PC2 loading")

ggs <- ggscatter(bh_df,x="rm",y="PC1",color=pal[1],shape=21,fill=adjustcolor(pal[1],alpha.f = 0.6),add="reg.line",conf.int = TRUE,add.params = list(color="red",fill="gray60"),size=4) +
  stat_regline_equation(label.x = 14,label.y=-0.12,aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),size=5,fontface="bold")
ggs1 <- ggpar(ggs,xlab="River mile",ylab="PC1 loading")

ggs2 <- ggscatter(bh_df,x="rm",y="PC2",color=pal[1],shape=21,fill=adjustcolor(pal[1],alpha.f = 0.6),add="reg.line",conf.int = TRUE,add.params = list(color="red",fill="gray60"),size=4) +
  stat_regline_equation(label.x = 13,label.y=-0.10,aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),size=5,fontface="bold")
ggs3 <- ggpar(ggs2,xlab="River mile",ylab="PC2 loading")



#################################################################
#filter vcf
keep <- which(missing < 0.7)
keep_cols <- c(1:9,keep+9)
vcf2 <- vcfR::read.vcfR("./../DAPC/sauger_0.5_0.03_nosex.vcf",cols=keep_cols)
vcfR::write.vcf(vcf2,file="./../DAPC/sauger_0.5_0.03_nosex_filt.vcf.gz")



##################################################################
#sex PCA for 0.5 0.03 dataset
vcf <- read.vcfR("./../pop_struct/all_saug/vcf_file/sauger_miss0.5_maf0.03.recode.vcf")
#vcf2geno("./../pop_struct/all_saug/vcf_file/sauger_miss0.5_maf0.03.recode.vcf")

source("D:/Masters/scripts/pca_funcs.R")

file <- read.table("./../DAPC/pntest_mean_sauger_0.5_0.03.txt",header=F)

meta <- read.csv("./../entropy/results_meta_sarwae_0.6_0.05.csv",header=T,stringsAsFactors = F)

vcf_names <- colnames(vcf@gt)[-1]
vcf_names <- gsub(".*aln_(EGM18_[[:digit:]]+).sorted.bam","\\1",vcf_names)

meta_filt <- meta[which(meta$EGM_ID %in% vcf_names),]
meta_filt <- meta_filt[match(vcf_names,meta_filt$EGM_ID),]
fishinfo <- meta_filt
colnames(fishinfo)[4] <- "Tributary"

#LEA::vcf2geno("./../DAPC/sauger_0.5_0.03_nosex.vcf")
gt <- LEA::read.geno("./../pop_struct/all_saug/vcf_file/sauger_miss0.5_maf0.03.recode.geno")
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

pal <- wesanderson::wes_palette("Zissou1",5)[c(3,1,5)]
{png(filename="report_pca_sex.png",width=6,height=6,units="in",res=300)
  plot(pca[[1]]$PC1,pca[[1]]$PC2,pch=21,col=pal[as.factor(fishinfo$Sex)],bg=adjustcolor(pal[as.factor(fishinfo$Sex)],alpha.f = 0.6),xlab=paste0("PC1 (",round(pca[[2]]$importance[2,1]*100,digits=2),"%)"),ylab=paste0("PC2 (",round(pca[[2]]$importance[2,2]*100,digits=2),"%)"))
  legend("bottomright", legend = c("Female","Male","Immature"), pch=21,col=pal[as.factor(unique(fishinfo$Sex)[-3])], pt.bg=adjustcolor(pal[as.factor(unique(fishinfo$Sex)[-3])],alpha.f = 0.6))
  dev.off()
}
