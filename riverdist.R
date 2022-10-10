#Attempting to quantify isolation by (river)distance in sauger

## reich fst estimator
## vectorized version
## input=genlight object
## FST will be calculated between pops in genlight object
## specify number of bootstraps using "bootstrap=100"
#from Jessi Rick
reich.fst <- function(gl, bootstrap=FALSE, verbose=TRUE) { 
  if (!require("matrixStats",character.only=T, quietly=T)) {
    install.packages("matrixStats")
    library(matrixStats, character.only=T)
  }
  
  nloc <- gl@n.loc
  npop <- length(levels(gl@pop))
  
  fsts <- matrix(nrow=npop,
                 ncol=npop,
                 dimnames=list(levels(gl@pop),levels(gl@pop)))
  
  if (bootstrap != FALSE){
    n.bs <- bootstrap
    bs <- data.frame(matrix(nrow=nrow(combinat::combn2(levels(gl@pop))),
                            ncol=n.bs+5))
  }
  
  k <- 0
  
  for (p1 in levels(gl@pop)){
    for (p2 in levels(gl@pop)){
      if (which(levels(gl@pop) == p1) < which(levels(gl@pop) == p2)) {
        k <- 1+k
        
        pop1 <- gl.keep.pop(gl, p1, mono.rm=FALSE, v=0)
        pop2 <- gl.keep.pop(gl, p2, mono.rm=FALSE, v=0)
        
        a1 <- colSums2(as.matrix(pop1),na.rm=T)
        a2 <- colSums2(as.matrix(pop2),na.rm=T)
        n1 <- apply(as.matrix(pop1),2,function(x) 2*sum(!is.na(x)))
        n2 <- apply(as.matrix(pop2),2,function(x) 2*sum(!is.na(x)))
        
        h1 <- (a1*(n1-a1))/(n1*(n1-1))
        h2 <- (a2*(n2-a2))/(n2*(n2-1))
        
        N <- (a1/n1 - a2/n2)^2 - h1/n1 - h2/n2
        D <- N + h1 + h2
        
        F <- sum(N, na.rm=T)/sum(D, na.rm=T)
        fsts[p2,p1] <- F
        if (verbose == TRUE) {
          print(paste("Pop1: ",p1,", Pop2: ",p2,", Reich FST: ",F,sep=""))
        }
        
        if (bootstrap != FALSE) {
          if (verbose == TRUE) {
            print("beginning bootstrapping")
          }
          
          bs[k,1:3] <- c(p2,p1,F)
          
          for (i in 1:n.bs){
            loci <- sample((1:nloc), nloc, replace=TRUE)
            
            pop1.bs <- as.matrix(pop1)[,loci]
            pop2.bs <- as.matrix(pop2)[,loci]
            
            a1 <- colSums2(as.matrix(pop1.bs),na.rm=T)
            a2 <- colSums2(as.matrix(pop2.bs),na.rm=T)
            n1 <- apply(as.matrix(pop1.bs),2,function(x) 2*sum(!is.na(x)))
            n2 <- apply(as.matrix(pop2.bs),2,function(x) 2*sum(!is.na(x)))
            
            h1 <- (a1*(n1-a1))/(n1*(n1-1))
            h2 <- (a2*(n2-a2))/(n2*(n2-1))
            
            N <- (a1/n1 - a2/n2)^2 - h1/n1 - h2/n2
            D <- N + h1 + h2
            
            F.bs <- sum(N, na.rm=T)/sum(D, na.rm=T)
            bs[k,i+5] <- F.bs
          }
          if (verbose == TRUE){
            print(paste("bootstrapping 95% CI: ",quantile(bs[k,6:(n.bs+5)],c(0.025),na.rm=T),"-",quantile(bs[k,3:n.bs+2],c(0.025),na.rm=T)))
          }
          
          bs[k,4:5] <- c(quantile(bs[k,6:n.bs+5],c(0.025),na.rm=T),
                         quantile(bs[k,6:n.bs+5],c(0.975),na.rm=T))
        }
        
      }
    }
  }
  
  #fsts[fsts < 0] <- 0
  colnames(bs)[1:5] <- c("pop1","pop2","fst_estimate","95_min_CI","95_max_CI")
  
  if (bootstrap != FALSE){
    fst.list <- list(fsts,bs)
    names(fst.list) <- c("fsts","bootstraps")
  } else {
    fst.list <- list(fsts)
    names(fst.list) <- "fsts"
  }
  
  return(fst.list)
}

#just Bighorn spawning sauger
#WRosenthal, 12/2/2019
library(sf)
library(sp)
library(riverdist)
library(ggmap)
library(GISTools)
library(rgdal)
library(mapview)
library(lwgeom)
library(ggpubr)

data <- readOGR(dsn="./sarwae_river.kml")
#writeOGR(data,dsn="temp",layer="sarwae_river",driver="ESRI Shapefile")
#projection uses Albers equal area
river <- line2network(path="./temp",layer="sarwae_river",reproject="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
plot(river)

sampling_locations <- read.csv("./sampling_locations_bighorn.csv")
colnames(sampling_locations)[1] <- "Location"
samp_lat <- as.character(sampling_locations$lat)
samp_lat <- as.numeric(char2dms(samp_lat,chd="d",chm="'",chs="\""))
samp_long <- as.character(sampling_locations$long)
samp_long <- as.numeric(char2dms(samp_long))
samp_long[29:33] <- -107.9472
samp_long[7:9] <- -108.1620

sampling_locations$lat <- samp_lat
sampling_locations$long <- samp_long
samp_locations <- sampling_locations
sampling_locations <- st_as_sf(sampling_locations,coords = c("long","lat"),crs = CRS("+proj=longlat +datum=WGS84 +no_defs "))

samp_loc_proj <- st_transform(sampling_locations,crs = 102003) # convert to Albers Equal Area
samp_loc_proj$X <- st_coordinates(samp_loc_proj)[,1]
samp_loc_proj$Y <- st_coordinates(samp_loc_proj)[,2]

snap_loc <- xy2segvert(x=samp_loc_proj$X,y=samp_loc_proj$Y,rivers=river)

distance_matrix <- riverdistancemat(seg=snap_loc$seg,vert = snap_loc$vert,rivers=river)

#now for the Fst calculation

library(vcfR)
library(dartR)

meta <- read.csv("./../entropy/results_meta_sarwae_0.6_0.05.csv",stringsAsFactors = F)
saug_meta <- meta[which(meta$q < 0.1),]
colnames(saug_meta)[1] <- "EGM_ID"

vcf <- read.vcfR("./../DAPC/sauger_0.5_0.03_nosex_filt.vcf")
names <- colnames(vcf@gt)
names <- as.vector(sapply(names,function(x) sub(".*aln_(EGM18_[[:digit:]]{4}).*","\\1",x)))[-1]
#names <- names[-569]
saug_meta <- saug_meta[which(saug_meta$EGM_ID %in% names),]

saug_meta <- saug_meta[grep("Bighorn",saug_meta$Location),]
saug_meta <- saug_meta[saug_meta$Temporal.Event == "Spawn",]
names2 <- names[names %in% saug_meta$EGM_ID]
saug_meta <- saug_meta[match(names2,saug_meta$EGM_ID),]


samp_loc <- paste0(samp_loc_proj$Location,samp_loc_proj$River.mile)
samp_loc <- gsub(" - ","-",samp_loc)


saug_loc <- paste0(saug_meta$Location,saug_meta$River_Mile)
saug_loc <- gsub("Popo Agie","Popo_Agie",saug_loc)
saug_loc <- gsub("Boysen","Boysen_Reservoir",saug_loc)
saug_loc <- gsub("Wind River","Wind_River",saug_loc)
saug_loc <- gsub("Little Wind_River","Little_Wind_River",saug_loc)
saug_loc <- gsub(" - ","-",saug_loc)
saug_loc <- gsub("51-50.5","50.5-51",saug_loc)
saug_loc <- gsub("70-68.5","68.5-70",saug_loc)
saug_loc <- gsub("48-46","46-48",saug_loc)
saug_loc <- gsub("52-50","50-52",saug_loc)
saug_loc <- gsub("54.5-52","52-54.5",saug_loc)
saug_loc <- gsub("68.5-67","67-68.5",saug_loc)

samp_loc <- unique(samp_loc[which(samp_loc %in% saug_loc ==TRUE)])

gl <- vcfR2genlight(vcf)
indNames(gl) <- as.vector(sapply(indNames(gl),function(x) sub(".*aln_(EGM18_[[:digit:]]{4}).*","\\1",x)))
gl <- gl.drop.ind(gl, indNames(gl)[which(indNames(gl) %in% saug_meta$EGM_ID == FALSE)])
pop(gl) <- saug_loc

fst <- reich.fst(gl,bootstrap = 200)

gl2 <- gl
pop(gl2) <- saug_meta$Location
fst2 <- reich.fst(gl2,bootstrap = 200)


dist_names <- paste0(samp_loc_proj$Location,samp_loc_proj$River.mile)[-14]
dist_names <- gsub(" - ","-",dist_names)

#saug_dist <- distance_matrix[-c(which(dist_names %in% samp_loc == FALSE),19),-c(which(dist_names %in% samp_loc == FALSE),19)]
saug_dist <- distance_matrix[which(dist_names %in% samp_loc),which(dist_names %in% samp_loc)]

ibd_info <- data.frame(saug_dist[which(lower.tri(saug_dist,diag=FALSE))]/1000,fst[[1]][which(lower.tri(fst[[1]],diag=FALSE))])
colnames(ibd_info) <- c("dist","fst")
ibd_info$fst_cor <- ibd_info$fst/(1-ibd_info$fst)

mantel_fst <- fst[[1]]/(1-fst[[1]])
diag(mantel_fst) <- 0
mantel_fst <- Matrix::forceSymmetric(mantel_fst,uplo = "L")

mantel_test_new <- ape::mantel.test(saug_dist,mantel_fst,graph=T)

ibdist <- lm(fst~dist,data=ibd_info)

ibdist_new <- lm(fst_cor~dist,data=ibd_info)

pal <- RColorBrewer::brewer.pal(n=3,"Dark2")
ggp <- ggscatter(ibd_info,x="dist",y="fst_cor",color=pal[1],shape=21,size=4,fill=adjustcolor(pal[1],alpha.f = 0.6))#add="reg.line",conf.int = TRUE,size=4,add.params = list(color="red",fill="gray60")) #+
  #stat_regline_equation(label.x = 175,label.y=-0.1,aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),size=5,fontface="bold")
ggp1 <- ggpar(ggp,xlab="Distance (river km)",ylab="Fst/(1-Fst)")
ggp2 <- ggtext(data.frame("x"=227,"y"=-0.12),x="x",y="y",label="Mantel test p-value = 0.917",ggp=ggp1,size=16)
ggp2



colnames(fst) <- samp_loc
rownames(fst) <- samp_loc
write.csv(fst,file="./bighorn_spawn_fst_maf.csv",row.names = T,quote=F)

#################################################
#Start conStruct analysis
#################################################

library(conStruct)

allele_freq <- counts/sample_sizes
saug_locations <- samp_locations[-14,]
saug_coords <- as.matrix(saug_locations[,3:4])
saug_coords[27,2] <- -107.9472
saug_dist <- saug_dist/1000

saug2 <- conStruct(spatial=T,K=2,freqs=allele_freq,coords=saug_coords,geoDist = saug_dist,n.iter=30000,n.chains=3,control = setNames(list(0.99,15),c("adapt_delta","max_treedepth")))
saug3 <- conStruct(spatial=T,K=3,freqs=allele_freq,coords=saug_coords,geoDist = saug_dist,n.iter=30000,n.chains=3,control = setNames(list(0.95,15),c("adapt_delta","max_treedepth")))
saug4 <- conStruct(spatial=T,K=4,freqs=allele_freq,coords=saug_coords,geoDist = saug_dist,n.iter=30000,n.chains=3,control = setNames(list(0.95,15),c("adapt_delta","max_treedepth")))
saug5 <- conStruct(spatial=T,K=5,freqs=allele_freq,coords=saug_coords,geoDist = saug_dist,n.iter=30000,n.chains=3,control = setNames(list(0.95,15),c("adapt_delta","max_treedepth")))



###############################################
#Fst calculation between Wind River Basin populations
###############################################


meta <- read.csv("./../entropy/results_meta_sarwae_0.6_0.05.csv",stringsAsFactors = F)
saug_meta <- meta[which(meta$q < 0.1),]
colnames(saug_meta)[1] <- "EGM_ID"

vcf <- read.vcfR("./../DAPC/sauger_0.5_0.03_nosex_filt.vcf")
names <- colnames(vcf@gt)
names <- as.vector(sapply(names,function(x) sub(".*aln_(EGM18_[[:digit:]]{4}).*","\\1",x)))[-1]
saug_meta <- saug_meta[which(saug_meta$EGM_ID %in% names),]

saug_meta <- saug_meta[which(saug_meta$Location %in% c("Wind River","Little Wind River","Popo Agie","Boysen")),]

samp_loc <- unique(saug_meta$Location)
saug_loc <- saug_meta$Location

gl <- vcfR2genlight(vcf)
indNames(gl) <- as.vector(sapply(indNames(gl),function(x) sub(".*aln_(EGM18_[[:digit:]]{4}).*","\\1",x)))
gl <- gl.drop.ind(gl, indNames(gl)[which(indNames(gl) %in% saug_meta$EGM_ID == FALSE)])
pop(gl) <- saug_loc


colnames(fst) <- samp_loc
rownames(fst) <- samp_loc
write.csv(fst,"./wind_basin_fst.csv",row.names=T,quote=F)

#between Wind River and Bighorn Basin samples
meta <- read.csv("./../entropy/results_meta_sarwae_0.6_0.05.csv",stringsAsFactors = F)
saug_meta <- meta[which(meta$q < 0.1),]
colnames(saug_meta)[1] <- "EGM_ID"

vcf <- read.vcfR("./../DAPC/sauger_0.5_0.03_nosex_filt.vcf")
names <- colnames(vcf@gt)
names <- as.vector(sapply(names,function(x) sub(".*aln_(EGM18_[[:digit:]]{4}).*","\\1",x)))[-1]
saug_meta <- saug_meta[which(saug_meta$EGM_ID %in% names),]
saug_meta <- saug_meta[match(names,saug_meta$EGM_ID),]

saug_meta$Location[grep("Bighorn",saug_meta$Location)] <- "Bighorn"
saug_meta$Location[grep("Wind",saug_meta$Location)] <- "Wind"
saug_meta$Location[grep("Boysen",saug_meta$Location)] <- "Wind"
saug_meta$Location[grep("Popo",saug_meta$Location)] <- "Wind"

samp_loc <- unique(saug_meta$Location)
saug_loc <- saug_meta$Location

gl <- vcfR2genlight(vcf)
indNames(gl) <- as.vector(sapply(indNames(gl),function(x) sub(".*aln_(EGM18_[[:digit:]]{4}).*","\\1",x)))
gl <- gl.drop.ind(gl, indNames(gl)[which(indNames(gl) %in% saug_meta$EGM_ID == FALSE)])
pop(gl) <- saug_loc
fst_basin <- reich.fst(gl,bootstrap = 200)

write.csv(fst,"./wind_bighorn_fst.csv",row.names = T,quote=F)




#IBD for all sauger
library(sf)
library(sp)
library(riverdist)
library(ggmap)
library(GISTools)
library(rgdal)
library(mapview)
library(lwgeom)

data <- readOGR(dsn="./sarwae_river.kml")
#writeOGR(data,dsn="temp",layer="sarwae_river",driver="ESRI Shapefile")
#projection uses Albers equal area
river <- line2network(path="./temp",layer="sarwae_river",reproject="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
plot(river)

sampling_locations <- read.csv("./saug_locations.csv")
colnames(sampling_locations)[1] <- "Location"
#samp_lat <- as.character(sampling_locations$lat)
#samp_lat <- as.numeric(char2dms(samp_lat,chd="d",chm="'",chs="\""))
#samp_long <- as.character(sampling_locations$long)
#samp_long <- as.numeric(char2dms(samp_long))
samp_lat <- sampling_locations$lat
samp_long <- sampling_locations$long
samp_long[27:33] <- -107.9472
samp_long[7:9] <- -108.1620

sampling_locations$lat <- samp_lat
sampling_locations$long <- samp_long
samp_locations <- sampling_locations
sampling_locations <- st_as_sf(sampling_locations,coords = c("long","lat"),crs = CRS("+proj=longlat +datum=WGS84 +no_defs "))

samp_loc_proj <- st_transform(sampling_locations,crs = 102003) # convert to Albers Equal Area
samp_loc_proj$X <- st_coordinates(samp_loc_proj)[,1]
samp_loc_proj$Y <- st_coordinates(samp_loc_proj)[,2]

snap_loc <- xy2segvert(x=samp_loc_proj$X,y=samp_loc_proj$Y,rivers=river)

distance_matrix <- riverdistancemat(seg=snap_loc$seg,vert = snap_loc$vert,rivers=river)

#now for the Fst calculation

library(BEDASSLE)
library(LEA)
library(vcfR)

meta <- read.csv("./../entropy/results_meta_sarwae_0.6_0.05.csv",stringsAsFactors = F)
saug_meta <- meta[which(meta$q < 0.1),]
colnames(saug_meta)[1] <- "EGM_ID"

vcf <- read.vcfR("./../DAPC/sauger_0.5_0.03_nosex_filt.vcf")
names <- colnames(vcf@gt)
names <- as.vector(sapply(names,function(x) sub(".*aln_(EGM18_[[:digit:]]{4}).*","\\1",x)))[-1]
#names <- names[-569]
saug_meta <- saug_meta[which(saug_meta$EGM_ID %in% names),]

#saug_meta <- saug_meta[grep("Bighorn",saug_meta$Location),]
#saug_meta <- saug_meta[saug_meta$Temporal.Event == "Spawn",]
names2 <- names[names %in% saug_meta$EGM_ID]
saug_meta <- saug_meta[match(names2,saug_meta$EGM_ID),]


samp_loc <- paste0(samp_loc_proj$Location,samp_loc_proj$River.mile)
samp_loc <- gsub(" - ","-",samp_loc)


saug_loc <- paste0(saug_meta$Location,saug_meta$River_Mile)
saug_loc <- gsub("Popo Agie","Popo_Agie",saug_loc)
saug_loc <- gsub("Boysen","Boysen_Reservoir",saug_loc)
saug_loc <- gsub("Wind River","Wind_River",saug_loc)
saug_loc <- gsub("Little Wind_River","Little_Wind_River",saug_loc)
saug_loc <- gsub(" - ","-",saug_loc)
saug_loc <- gsub("51-50.5","50.5-51",saug_loc)
saug_loc <- gsub("70-68.5","68.5-70",saug_loc)
saug_loc <- gsub("48-46","46-48",saug_loc)
saug_loc <- gsub("52-50","50-52",saug_loc)
saug_loc <- gsub("54.5-52","52-54.5",saug_loc)
saug_loc <- gsub("68.5-67","67-68.5",saug_loc)

samp_loc <- unique(samp_loc[which(samp_loc %in% saug_loc ==TRUE)])

gl <- vcfR2genlight(vcf)
indNames(gl) <- as.vector(sapply(indNames(gl),function(x) sub(".*aln_(EGM18_[[:digit:]]{4}).*","\\1",x)))
gl <- gl.drop.ind(gl, indNames(gl)[which(indNames(gl) %in% saug_meta$EGM_ID == FALSE)])
pop(gl) <- saug_loc

fst_all <- reich.fst(gl,bootstrap=100,verbose = FALSE)

dist_names <- paste0(samp_loc_proj$Location,samp_loc_proj$River.mile)
dist_names <- gsub(" - ","-",dist_names)

#saug_dist <- distance_matrix[-c(which(dist_names %in% samp_loc == FALSE),19),-c(which(dist_names %in% samp_loc == FALSE),19)]
saug_dist <- distance_matrix[which(dist_names %in% samp_loc),which(dist_names %in% samp_loc)]

ibd_info <- data.frame(saug_dist[which(lower.tri(saug_dist,diag=FALSE))]/1000,fst_all[[1]][which(lower.tri(fst_all[[1]],diag=FALSE))])
colnames(ibd_info) <- c("dist","fst")
ibd_info$fst_cor <- ibd_info$fst/(1-ibd_info$fst)

ibdist <- lm(fst~dist,data=ibd_info)

ibdist_new <- lm((fst/(1 - fst))~dist,data=ibd_info)

mantel_fst <- fst_all[[1]]/(1-fst_all[[1]])
diag(mantel_fst) <- 0
mantel_fst <- Matrix::forceSymmetric(mantel_fst,uplo = "L")

mantel_test_new <- ape::mantel.test(saug_dist,mantel_fst,graph=T)
mantel_test <- ape::mantel.test(saug_dist,(fst),graph=TRUE)

col <- RColorBrewer::brewer.pal(n=3,"Dark2")[1]

{png(file="./saug_IBD.png",width=6,height=6,units="in",res=300)
  plot(ibd_info$dist,(ibd_info$fst/(1-ibd_info$fst)),pch=21,col=col,bg=adjustcolor(col,alpha.f = 0.8),ylab="Fst/(1-Fst)",xlab="distance (km)",sub=paste0("Mantel test p-value=",mantel_test_new$p))
  abline(ibdist_new)
  dev.off()
}

pal <- RColorBrewer::brewer.pal(n=8,"Dark2")
ggp <- ggscatter(ibd_info,x="dist",y="fst_cor",color=pal[6],shape=21,fill=adjustcolor(pal[6],alpha.f = 0.6),add="reg.line",conf.int = TRUE,size=4,add.params = list(color="red",fill="gray60")) +
    stat_regline_equation(label.x = 203,label.y=-0.12,aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),size=5,fontface="bold")
ggp1 <- ggpar(ggp,xlab="Distance (river km)",ylab="Fst/(1 - Fst)")
ggp2 <- ggtext(data.frame("x"=288,"y"=-0.14),x="x",y="y",label="Mantel test p-value = 0.001",ggp=ggp1,size=16)
ggp2


colnames(fst) <- samp_loc
rownames(fst) <- samp_loc
write.csv(fst,file="./saug_fst.csv",row.names = T,quote=F)

