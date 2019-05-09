computerFlag="cluster"
computerFlag=""

if (computerFlag=="cluster") {
} else {
    dirSrc="/Users/royr/UCSF/"
    dirSrc2=dirSrc
    setwd(paste(dirSrc2,"JoeWiemels/leukMeth/normalize",sep=""))
}

####################################################################
####################################################################

cat("\n================== normalize_epic.R ===============\n\n")


##########################################
## Section 1

## ----------------------------------

library(minfi)

if (F) {
    library(minfiData)
    library(RColorBrewer)
    library(matrixStats)
    library(pheatmap)
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    library(IlluminaHumanMethylation450kmanifest)
    library(IlluminaHumanMethylationEPICmanifest)
    #library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
    library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    library(FlowSorted.Blood.450k)
    library(sva)
    library(limma)
    library(DMRcate)
    library(wateRmelon)
    library(ENmix)
    library(gplots)

    cat("------------------- gplots done ---------------\n\n")

    ## ----------------------------------

    if (T) {
        #library(FlowSorted.Blood.EPIC)
        library(ExperimentHub)
        hub <- ExperimentHub()
        query(hub, "FlowSorted.Blood.EPIC")
        FlowSorted.Blood.EPIC <- hub[["EH1136"]]
        FlowSorted.Blood.EPIC
        save(FlowSorted.Blood.EPIC,file="FlowSorted.Blood.EPIC.RData")
    }

    load(file="FlowSorted.Blood.EPIC.RData")


    ## ----------------------------------

    library(RPMM)

}

## ----------------------------------
## Parameters

cat("------------------- Parameters ---------------\n\n")

cohortFlag="_periDsal"

dirOut="periDsal/"

nBatch=4
nBatch=1
nBatch=10
nBatch=20

cexLeg=1

## ----------------------------------
#Reading in the data

##Set the data directory to the folder containing the Idat files as well as the Sample Excel sheet
if (computerFlag=="cluster") {
    dirMeth="/cbc2/data1/ritu/wiemelsJ_periDsal/idat/"
    dirClin="/home/royr/project/JoeWiemels/data/periDsal/"
    fNameClin="IDATsControls and other QC PERI-DSAL_UCB EPIC plates 3-12_6-20-2018HMH_.txt"
} else {
    dirMeth="/Users/royr/Downloads/IDATs/"
    dirMeth="/Users/royr/UCSF/JoeWiemels/leukMeth/docs/periDsal/idat/"
    dirClin="/Users/royr/UCSF/JoeWiemels/leukMeth/docs/periDsal/PERI-DSAL Sample QC 2018_UCB EPIC methylation plates 3-12/"
    fNameClin="IDATsControls and other QC PERI-DSAL_UCB EPIC plates 3-12_6-20-2018HMH_.txt"
}

clin=read.table(paste(dirClin,fNameClin,sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
names(clin)[match(c("Sample.Name","Subject.ID","Sentrix.Barcode","Sentrix.Position","array_pos","Plate...QC.batch..AorB.","SD1.Itemno.plated","X..Methylation.Detected.at.0.01","X..Methylation.Detected.at.0.05","Genomestudio.Predicted.Gender","Reported.Gender","MATCH.","PERI.DSAL","Case.control","Restoration","Staining.Green","Staining.Red","Extension.Green","Extension.Red","Hybridization.High.Medium","Hybridization.Medium.Low","Target.Removal.1","Target.Removal.2","Bisulfite.Conversion.1.Green","Bisulfite.Conversion.1.Background.Green","Bisulfite.Conversion.1.Red","Bisulfite.Conversion.1.Background.Red","Bisulfite.Conversion.2","Bisulfite.Conversion.2.Background","Specificity.1.Green","Specificity.1.Red","Specificity.2","Specificity.2.Background","Non.Polymorphic.Green","Non.Polymorphic.Red"),names(clin))]=
c("sampleName","subjectId","beadchip","position","arrayPos","plate","sd1","methDet0.01","methDet0.05","sexPredGenStudio","sex","sexMatched","periDsal","caco","restoration","stainingGreen","stainingRed","extensionGreen","extensionRed","hybHighMed","hybMedLow","targetRem1","targetRem2","biConv1Grn","biConv1BgndGrn","biConv1Red","biConv1BgndRed","biConv2","biConv2Bgnd","spec1Grn","spec1Red","spec2","spec2Bgnd","nonPolyGrn","nonPolyRed")
clin$id=paste("X",clin$sampleName,sep="")
clin$periDsal=tolower(clin$periDsal)
clin$periDsal[which(!clin$periDsal%in%c("dsal","peri"))]=NA
clin$caco[which(!clin$caco%in%c("case","control"))]=NA
clin$sexPredGenStudio=tolower(clin$sexPredGenStudio)
clin$sex=tolower(clin$sex)
clin$sex[which(!clin$sex%in%c("f","m"))]=NA
#clin$Basename=paste("/",clin$beadchip,"/",clin$beadchip,"_",clin$position,sep="")
clin$Basename=paste(dirMeth,clin$beadchip,"/",clin$beadchip,"_",clin$position,sep="")

clin1=clin

if (computerFlag=="cluster") {
    grpUniq=dir(dirMeth)
} else {
    grpUniq=1:118
    grpUniq=dir(dirMeth)
    nBatch=1
    grpUniq=grpUniq[1:2]
}

set.seed(54534)
grpId=sample(length(grpUniq),replace=F)
grpVec=round(seq(1,(length(grpId)+1),length=nBatch+1))
grpVec=cbind(grpVec[1:(length(grpVec)-1)],grpVec[2:length(grpVec)]-1)

cat("------------------- Section 1 End ---------------\n\n")

##########################################
## Section 2

if (F) {

for (gvId in 1:nrow(grpVec)) {
    clin=clin1[which(clin1$beadchip%in%c(grpUniq[grpVec[gvId,1]:grpVec[gvId,2]])),]
    grp=c()
    for (gId in 1:length(grpUniq)) {
        x=dir(paste(dirMeth,grpUniq[gId],sep=""))
        x=x[grep("_Grn",x)]
        grp=c(grp,sapply(x,function(x) paste(strsplit(x,"_")[[1]][1:2],collapse="_"),USE.NAMES=F))
    }
    #clin=clin[which(clin$beadchip%in%dir(dirMeth)),]
    clin=clin[which(paste(clin$beadchip,clin$position,sep="_")%in%grp),]

    targets=clin

    RGset <- read.metharray.exp(targets = targets, recursive=T,force = TRUE) #Reads the IDAT files based on the sample sheet and creates an RGset
    #RGset <- read.metharray.exp(base=dirMeth,targets = targets, recursive=T,force = TRUE) #Reads the IDAT files based on the sample sheet and creates an RGset
    #RGset=read.metharray.exp(base=dirMeth, recursive=T, targets=targets)
    save(RGset,file=paste("RGset_",gvId,".RData",sep=""))
}

#save.image("tmp11.RData")

## NOT USED
if (F) {
    if (computerFlag=="") {
        #targets=data.frame(Sample_Name=,Sample_Group=,BeadChip=,Position=,stringsAsFactors=F)
        targets=as.data.frame(t(sapply(colnames(RGset),function(x) {
            y=strsplit(x,"_")[[1]]
            beadchip=y[1]
            position=y[2]
            id=paste("X",y[1],sep="")
            c(id,beadchip,position)
        },USE.NAMES=T)),stringsAsFactors=F)
        names(targets)=c("id","beadchip","position")
        targets$Sample_Group="grp1"
        targets$Sample_Group[c(2,3,4)]="grp2"
        targets$Basename=paste("/",targets$beadchip,"/",targets$beadchip,"_",targets$position,sep="")
    }
}

x=rep(NA,nrow(grpVec))
for (gvId in 1:nrow(grpVec)) {
    load(file=paste("RGset_",gvId,".RData",sep=""))
    x[gvId]=ncol(RGset)
}

gvId=1
load(file=paste("RGset_",gvId,".RData",sep=""))
if (nBatch>1) {
    RGsetTmp=RGset
    for (gvId in 2:nrow(grpVec)) {
        load(file=paste("RGset_",gvId,".RData",sep=""))
        RGsetTmp=cbind(RGsetTmp,RGset)
    }
    RGset=RGsetTmp
    rm(RGsetTmp)
}
targets=clin1[match(sampleNames(RGset),paste(clin1$beadchip,clin1$position,sep="_")),]

RGset
#targets$id <- targets$Sample_Name
sampleNames(RGset) <- targets$id
pData(RGset)=DataFrame(targets)
sampleNames(RGset)=pData(RGset)$id
RGset@annotation[2]="ilm10b4.hg19"

save(RGset,file="RGset.RData")
rm(RGset)

for (gvId in 1:nrow(grpVec)) {
    clin=clin1[which(clin1$beadchip%in%c(grpUniq[grpVec[gvId,1]:grpVec[gvId,2]])),]
    grp=c()
    for (gId in 1:length(grpUniq)) {
        x=dir(paste(dirMeth,grpUniq[gId],sep=""))
        x=x[grep("_Grn",x)]
        grp=c(grp,sapply(x,function(x) paste(strsplit(x,"_")[[1]][1:2],collapse="_"),USE.NAMES=F))
    }
    clin=clin[which(paste(clin$beadchip,clin$position,sep="_")%in%grp),]
    
    targets=clin
    
    RGset1 <- read.metharray.exp(targets = targets , extended = TRUE, force = TRUE)
    #RGset1=read.metharray.exp(base=dirMeth, recursive=T, extended=T)
    save(RGset1,file=paste("RGset1_",gvId,".RData",sep=""))
}

gvId=1
load(file=paste("RGset1_",gvId,".RData",sep=""))
if (nBatch>1) {
    RGsetTmp=RGset1
    for (gvId in 2:nrow(grpVec)) {
        load(file=paste("RGset1_",gvId,".RData",sep=""))
        RGsetTmp=cbind(RGsetTmp,RGset1)
    }
    RGset1=RGsetTmp
    rm(RGsetTmp)
}
targets=clin1[match(sampleNames(RGset1),paste(clin1$beadchip,clin1$position,sep="_")),]

#targets$id <- targets$Sample_Name
sampleNames(RGset1) <- targets$id
pData(RGset1)=DataFrame(targets)
sampleNames(RGset1)=pData(RGset1)$id

save(RGset1,file="RGset1.RData")
rm(RGset1)

## NOT USED
if (F) {
    targets <- read.metharray.sheet(dataDirectory , pattern = "Sample_Sheet_Full.csv") #Reads in the Excel sheet
    head(targets)
    RGset <- read.metharray.exp(targets = targets, force = TRUE) #Reads the IDAT files based on the sample sheet and creates an RGset
    RGset1 <- read.metharray.exp(targets = targets , extended = TRUE, force = TRUE)
}


load(file="RGset.RData")
RGset
#names <- pData(RGset)$Sample_Name
#names <- pData(RGset)$id
groups <- pData(RGset)$periDsal
#beadchip <- pData(RGset)$beadchip
#manifest <- getManifest(RGset)


## ----------------------------------
## Pre-processing
## Step 1 - QC

cat("\n------------- Start QCinfo -------------\n\n")
library(ENmix)
load(file="RGset1.RData")
res=QCinfo(RGset1, detPthre=0.05, nbthre=3, CpGthre=0.05,samplethre=0.1,outlier=TRUE, distplot=T)
rm(RGset1)
cat("\n------------- End QCinfo -------------\n\n")

##########################################
## Section 3

cat("------------------- Section 3 Start ---------------\n\n")

load(file="RGset.RData")
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
save(ann850k,file="ann850k.RData")

library(pheatmap)
load(file="RGset.RData")
betaSnp=getSnpBeta(RGset)
save(betaSnp,file="betaSnp.RData")
pdf(paste(dirOut,"heatmap_snpBeta",cohortFlag,".pdf",sep=""))
pheatmap(betaSnp)
dev.off()

library(ENmix)
load(file="RGset.RData")
#pdf(paste(dirOut,"plotCtrl",cohortFlag,".pdf",sep=""))
plotCtrl(RGset,IDorder=NULL)
#dev.off()

#load(file="RGset.RData")
Mset <- preprocessRaw(RGset)
save(Mset,file="Mset.RData")
rm(RGset)
pdf(paste(dirOut,"plotQC",cohortFlag,".pdf",sep=""))
plotQC(getQC(Mset))
dev.off()
rm(Mset)

cat("------------------- Section 3 End ---------------\n\n")

}

##########################################
## Section 3.2

if (F) {
    
cat("------------------- Section 3.2 start ---------------\n\n")

load(file="RGset.RData")
beta=getBeta(RGset)
save(beta,file="beta_raw.RData")

load(file="RGset.RData")
res <- mapToGenome(RGset)
save(res,file="mapToGenome.RData")
rm(RGset)

load(file="Mset.RData")
cat("\n------------- Start QCfilter -------------\n\n")
RGset2 <- QCfilter(Mset,qcinfo=res,detPthre=0.000001,nbthre=3,samplethre=0.05,CpGthre=0.05,bisulthre=NULL,outlier=FALSE,outid=NULL, outCpG=NULL,plot=FALSE)
save(RGset2,file="RGset2.RData")
cat("\n------------- End QCfilter -------------\n\n")
rm(Mset,res,RGset2)

#save.image("tmp1.RData")

}

##########################################
## Section 4

cat("------------------- funNorm normalization start ---------------\n\n")

load(file="RGset.RData")
rgsetFN=preprocessFunnorm(RGset)
save(rgsetFN,file="rgsetFN.RData")
betaFN <- minfi::getBeta(rgsetFN)
save(betaFN,file="betaFN.RData")

#Plot raw and funNorm Normalized values in one plot
library(RColorBrewer)
load(file="RGset.RData")
targets=as.data.frame(pData(RGset))
load(file="betaFN.RData")
pdf(paste(dirOut,"densityPlot_funNorm",cohortFlag,".pdf",sep=""))
par(mfrow=c(1,2))
densityPlot(RGset , sampGroups = targets$periDsal , main = "Raw" , legend = FALSE, xlab = "Beta Values")
legend("top",cex=cexLeg,legend=levels(factor(targets$periDsal)) , text.col = brewer.pal(8,"Dark2"))
densityPlot(betaFN , sampGroups = targets$periDsal , main = "funNorm Normalized" , legend = FALSE, xlab = "Beta Values")
legend("top",cex=cexLeg,legend=levels(factor(targets$periDsal)) , text.col = brewer.pal(8,"Dark2"))
dev.off()

cat("------------------- funNorm normalization end ---------------\n\n")


if (F) {

cat("------------------- Normalization Step1 - Noob Normalization start ---------------\n\n")

load(file="RGset.RData")
targets=as.data.frame(pData(RGset))
RGset_Noob <- preprocessNoob(RGset , verbose = TRUE)
save(RGset_Noob,file="RGset_Noob.RData")
beta_Noob <- minfi::getBeta(RGset_Noob)
save(beta_Noob,file="beta_Noob.RData")
load(file="beta_Noob.RData")
sd_beta_Noob=apply(beta_Noob,1,sd,na.rm=T)
save(sd_beta_Noob,file="sd_beta_Noob.RData")

if (F) {
    ## Distribution of beta-values
    load(file="beta_Noob.RData")
    res=t(apply(beta_Noob,2,function(x) {
        c(quantile(x,probs=c(0,.25,.5,.75,1),na.rm=T),mean(x,na.rm=T))
    }))
    colnames(res)=c(paste("perc_",c(0,25,50,75,100),sep=""),"mean")
    save(res,file="stat_sample_beta_noob_sample.RData")

    load(file="stat_sample_beta_noob_sample.RData")
    clin=clin1[match(rownames(res),clin1$id),]
    j=which(!is.na(clin$periDsal))
    png("densityPlot_stat_sample_betaNoob_sample.png")
    par(mfrow=c(2,2))
    for (k in c("perc_25","perc_50","perc_75","mean")) {
        plot(density(res[j,k],na.rm=T),main=k,xlab="Noob-normalized beta-values")
    }
    dev.off()
}

## M-value matrix
offset=0.0000001
mVal=beta_Noob
rm(beta_Noob,sd_beta_Noob)
mVal[which(mVal==0)]=offset; mVal[which(mVal==1)]=1-offset
mVal=log2(mVal/(1-mVal))
save(mVal,file="mVal_noob.RData")
sd_mVal=apply(mVal,1,sd,na.rm=T)
save(sd_mVal,file="sd_mVal_noob.RData")
rm(mVal,sd_mVal)

# Plot raw and Noob Normalized values in one plot
library(RColorBrewer)
load(file="RGset.RData")
pdf(paste(dirOut,"densityPlot_noob",cohortFlag,".pdf",sep=""))
par(mfrow=c(1,2))
densityPlot(RGset , sampGroups = targets$periDsal , main = "Raw" , legend = FALSE, xlab = "Beta Values")
legend("top",cex=cexLeg,legend=levels(factor(targets$periDsal)) , text.col = brewer.pal(8,"Dark2"))
densityPlot(beta_Noob , sampGroups = targets$periDsal , main = "Noob Normalized" , legend = FALSE, xlab = "Beta Values")
legend("top",cex=cexLeg,legend=levels(factor(targets$periDsal)) , text.col = brewer.pal(8,"Dark2"))
dev.off()

# Plot raw and Noob Normalized values in one plot. Mark bad samples
load("results/stat_sample_beta_noob_sample.RData")
res=res[match(targets$id,rownames(res)),]
pdf(paste(dirOut,"plotMedian_noob",cohortFlag,"_QC.pdf",sep=""))
plot(sort(res[,"perc_50"]))
dev.off()
pdf(paste(dirOut,"plotIQR_noob",cohortFlag,"_QC.pdf",sep=""))
plot(sort(res[,"perc_75"]-res[,"perc_25"]))
abline(h=c(0.645,0.68),lty="dotted")
dev.off()
thres=c(0.65,0.60)
grp=rep("0. good",nrow(res)); grp[which(res[,"perc_50"]<thres[1])]="1. low 1"; grp[which(res[,"perc_50"]<thres)]="2. low 2"
pdf(paste(dirOut,"densityPlot_noob",cohortFlag,"_QC.pdf",sep=""))
par(mfrow=c(1,2))
#densityPlot(RGset , sampGroups = grp , main = "Raw" , legend = FALSE, xlab = "Beta Values")
#legend("top",cex=cexLeg,legend=levels(factor(grp)) , text.col = brewer.pal(8,"Dark2"))
densityPlot(beta_Noob , sampGroups = grp , main = "Noob Normalized" , legend = FALSE, xlab = "Beta Values")
legend("top",cex=cexLeg,legend=levels(factor(grp)) , text.col = brewer.pal(8,"Dark2"))
dev.off()

write.table(cbind(targets,res,quality=grp),file=paste(dirOut,"stat_sample_betaNoob",cohortFlag,"_QC.txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F,append=F)

## ------------
clin=read.table(paste(dirOut,"stat_sample_betaNoob",cohortFlag,"_QC.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
if (F) {
    clin$patientId=sapply(clin$sampleName,function(x) {
        samId=grep("-",x)
        if (substr(x,1,nchar("Jurkat"))=="Jurkat") {
            patId=x
        } else {
            patId=strsplit(x,"-")[[1]][1]
        }
        patId
    },USE.NAMES=F)
    j=which(clin$patientId%in% clin$patientId[duplicated(clin$patientId)])
    clin$hasReplicate="no"; clin$hasReplicate[j]="yes"
    nm=c("sampleName","patientId","hasReplicate")
    tbl=clin[order(2-as.integer(as.factor(clin$hasReplicate)),clin$patientId,decreasing=F),c(nm,names(clin)[!names(clin)%in%nm])]
    write.table(tbl,file=paste(dirOut,"stat_sample_betaNoob",cohortFlag,"_withReplFlag.txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F,append=F)
}

thres=118
thres=20
tmpC=character(ncol(clin))
varInfo=data.frame(variable=names(clin),varName=names(clin),class=tmpC,type=tmpC,stringsAsFactors=F)
varInfo$varName[match(c("sampleName","subjectId","beadchip","position","arrayPos","plate","sd1","methDet0.01","methDet0.05","sexPredGenStudio","sex","sexMatched","periDsal","caco","restoration","stainingGreen","stainingRed","extensionGreen","extensionRed","hybHighMed","hybMedLow","targetRem1","targetRem2","biConv1Grn","biConv1BgndGrn","biConv1Red","biConv1BgndRed","biConv2","biConv2Bgnd","spec1Grn","spec1Red","spec2","spec2Bgnd","nonPolyGrn","nonPolyRed"),varInfo$variable)]=c("Sample.Name","Subject.ID","Sentrix.Barcode","Sentrix.Position","array_pos","Plate...QC.batch..AorB.","SD1.Itemno.plated","X..Methylation.Detected.at.0.01","X..Methylation.Detected.at.0.05","Genomestudio.Predicted.Gender","Reported.Gender","MATCH.","PERI.DSAL","Case.control","Restoration","Staining.Green","Staining.Red","Extension.Green","Extension.Red","Hybridization.High.Medium","Hybridization.Medium.Low","Target.Removal.1","Target.Removal.2","Bisulfite.Conversion.1.Green","Bisulfite.Conversion.1.Background.Green","Bisulfite.Conversion.1.Red","Bisulfite.Conversion.1.Background.Red","Bisulfite.Conversion.2","Bisulfite.Conversion.2.Background","Specificity.1.Green","Specificity.1.Red","Specificity.2","Specificity.2.Background","Non.Polymorphic.Green","Non.Polymorphic.Red")
for (k in 1:ncol(clin)) {
    x=clin[!is.na(clin[,k]),k]
    varInfo$class[k]=class(x)
    if (is.numeric(x)) {
        if (varInfo$variable[k]!="beadchip" & sum(!duplicated(x))>thres) {
            varInfo$type[k]="continuous"
        } else {
            varInfo$type[k]="categorical"
        }
    } else {
        if (sum(!duplicated(x))<=thres) {
            varInfo$type[k]="categorical"
        }
    }
}
varInfo[which(varInfo$type==""),]
library(coin)
tbl1=NULL
nm=c("variable1","variable2","testType","pv")
pThres=0.05
samId=which(!is.na(clin$periDsal) & !is.na(clin$caco) & !duplicated(clin$subjectId))
varList1="quality"
varList2=varInfo$variable[which(varInfo$type=="continuous" & !varInfo$variable%in%varList1)]
for (vId1 in 1:length(varList1)) {
    #png(paste("boxplot_",varList1[vId1],"_%01d.png",sep=""),width=6*240,height=4*240)
    png(paste("boxplot_",varList1[vId1],".png",sep=""),width=6*240,height=4*240)
    par(mfrow=c(2,3))
    par(mar=c(5, 5, 4, .5) + 0.1)
    for (vId2 in 1:length(varList2)) {
        if (varList1[vId1]=="quality" & (substr(varList2[vId2],1,nchar("perc_"))=="perc_" | varList2[vId2]=="mean")) next
        x=as.factor(clin[samId,varList1[vId1]])
        pv=pvalue(kruskal_test(clin[samId,varList2[vId2]]~x))
        tbl2=data.frame(varList1[vId1],varList2[vId2],"kruskal",pv,stringsAsFactors=F)
        names(tbl2)=nm
        tbl1=rbind(tbl1,tbl2)
        header="Kuskal-wallis test pv "
        if (pv<pThres) boxplot(clin[samId,varList2[vId2]]~x,main=paste(header,signif(pv,2),sep=""),xlab=varInfo$varName[match(varList1[vId1],varInfo$variable)],ylab=varInfo$varName[match(varList2[vId2],varInfo$variable)],cex.main=3,cex.lab=2,cex.axis=2)
    }
    dev.off()
}
varList1="quality"
varList2=varInfo$variable[which(varInfo$type=="categorical" & !varInfo$variable%in%varList1)]
for (vId1 in 1:length(varList1)) {
    for (vId2 in 1:length(varList2)) {
        if (varList1[vId1]=="quality" & (substr(varList2[vId2],1,nchar("perc_"))=="perc" | varList2[vId2]=="mean")) next
        pv=fisher.test(clin[samId,varList1[vId1]],clin[samId,varList2[vId2]])$p.value
        tbl2=data.frame(varList1[vId1],varList2[vId2],"fisher",pv,stringsAsFactors=F)
        names(tbl2)=nm
        tbl1=rbind(tbl1,tbl2)
    }
}
tbl1[which(tbl1$pv<pThres),]
grp=paste(clin$periDsal,"/",clin$caco,sep="")
table(quality=clin$quality[samId],grp[samId])
tbl1$adjP=NA
k=which(!is.na(tbl1$pv))
tbl1$adjP[k]=p.adjust(tbl1$pv,method="BH")
tbl1[which(tbl1$adjP<pThres),]

## ----------------
## Distribution of bad samples

load(file="beta_Noob.RData")
x=beta_Noob[,match(targets$id,colnames(beta_Noob))]
rm(beta_Noob)
x=x[,which(grp!="0. good")]
out=apply(x,1,function(x) {
    median(x,na.rm=T)
})
y=res[match(targets$id,rownames(res)),][which(grp!="0. good"),]

load(file="ann850k.RData")
ann=as.data.frame(ann850k)
rm(ann850k)
i=which(!rownames(x)%in%ann$Name)
tmp=ann[1:length(i),]
for (k in 1:ncol(tmp)) tmp[,k]=NA
tmp$Name=rownames(x)[i]
ann=rbind(ann,tmp)
iA=match(rownames(x),ann$Name)
i=which(out>=.5 & out<=.6)
table(ann$chr[iA][i],exclude=NULL)

write.table(ann[iA[i],],file=paste(dirOut,"ann_hemiMeth_betaNoob",cohortFlag,"_QC.txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F,append=F)

}

## ----------------------------------
## Step 2 - Filtering

if (F) {

cat("------------------- Step 2 - Filtering start ---------------\n\n")

load(file="RGset.RData")

#filtering using p-values
load(file="RGset.RData")
detectP <- minfi::detectionP(RGset) # set detection cutoff.
save(detectP,file="detectP.RData")
rm(RGset)
failed.05 <- detectP > 0.05
failedProbes <- rownames(failed.05)[rowMeans(failed.05) > 0.2] # list of probes that failed in more than 20% of the sample
load("RGset_Noob.RData")
RGset_Noob_filt <- RGset_Noob[which(!rownames(RGset_Noob) %in% failedProbes)]
save(RGset_Noob_filt,file="RGset_Noob_filt_1.RData")
cat("The number of probes failed in more than 20% of samples:\n")
print(sum(rowMeans(failed.05)>0.2))
cat("Dimensions of the original dataset:\n")
print(dim(RGset_Noob))
cat("Dimensions after filtering out failed probes:\n")
print(dim(RGset_Noob_filt))
rm(RGset_Noob)

#save.image("tmp2.RData")

## 2nd Filtering Step
#Convert to GenomicRatio Set
load(file="RGset_Noob_filt_1.RData")
GRset_Noob_filt <- mapToGenome(RGset_Noob_filt)
save(GRset_Noob_filt,file="GRset_Noob_filt_1.RData")
rm(RGset_Noob_filt)
#Sex Prediction
load("GRset_Noob_filt_1.RData")
predictedSex <- getSex(GRset_Noob_filt, cutoff = -2)$predictedSex
load("GRset_Noob_filt_1.RData")
names(predictedSex)=sampleNames(GRset_Noob_filt)
save(predictedSex,file="predictedSex.RData")
#head(predictedSex)
print(table(observedSex=pData(GRset_Noob_filt)$sex,predictedSex,exclude=NULL))
rm(GRset_Noob_filt)
#check
#plotSex(getSex(GRset_Noob_filt, cutoff = -2))
if (F) {
    ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    save(ann850k,file="ann850k.RData")
} else {
    load(file="ann850k.RData")
}
load("RGset_Noob_filt_1.RData")
keepSex <- !(featureNames(RGset_Noob_filt) %in% ann850k$Name[ann850k$chr %in% c("chrX","chrY")])
rm(ann850k)
print(table(keepSex))
RGset_Noob_filt <- RGset_Noob_filt[keepSex,]
save(RGset_Noob_filt,file="RGset_Noob_filt_2.RData")
RGset_Noob_filt <- dropMethylationLoci(RGset_Noob_filt, dropRS = TRUE)
save(RGset_Noob_filt,file="RGset_Noob_filt_3.RData")
rm(RGset_Noob_filt)
#SNP Prediction and Filtering
load("GRset_Noob_filt_1.RData")
SNPs <- getSnpInfo(GRset_Noob_filt)
save(SNPs,file="SNPs.RData")
print(head(SNPs, 10))
rm(SNPs)
GRset_Noob_filt <- addSnpInfo(GRset_Noob_filt)
save(GRset_Noob_filt,file="GRset_Noob_filt_2.RData")
GRset_Noob_filt <- dropLociWithSnps(GRset_Noob_filt)
save(GRset_Noob_filt,file="GRset_Noob_filt_3.RData")
rm(GRset_Noob_filt)

load("RGset_Noob_filt_1.RData")
targets=as.data.frame(pData(RGset_Noob_filt))
beta_RGset_noob_filt <- minfi::getBeta(RGset_Noob_filt)
save(beta_RGset_noob_filt,file="beta_RGset_noob_filt_1.RData")
offset=0.0000001
mVal=beta_RGset_noob_filt
rm(beta_RGset_noob_filt)
mVal[which(mVal==0)]=offset; mVal[which(mVal==1)]=1-offset
mVal=log2(mVal/(1-mVal))
save(mVal,file="mVal_RGset_noob_filt_1.RData")
sd_mVal=apply(mVal,1,sd,na.rm=T)
save(sd_mVal,file="sd_mVal_RGset_noob_filt_1.RData")
rm(mVal,sd_mVal)
library(RColorBrewer)
pdf(paste(dirOut,"densityPlot_noob_filt",cohortFlag,".pdf",sep=""))
par(mfrow=c(1,2))
#densityPlot(RGset , sampGroups = targets$periDsal , main = "Raw" , legend = FALSE, xlab = "Beta Values")
#legend("top",cex=cexLeg,legend=levels(factor(targets$periDsal)) , text.col = brewer.pal(8,"Dark2"))
#densityPlot(beta_Noob , sampGroups = targets$periDsal , main = "Noob Normalized & filtered" , legend = FALSE, xlab = "Beta Values")
densityPlot(RGset_Noob_filt , sampGroups = targets$periDsal , main = "Noob Normalized & filtered" , legend = FALSE, xlab = "Beta Values")
legend("top",cex=cexLeg,legend=levels(factor(targets$periDsal)) , text.col = brewer.pal(8,"Dark2"))
dev.off()

#save.image("tmp3.RData")

}

## ----------------------------------
## Normalization Step2 - BMIQ Normalization

cat("------------------- Normalization Step2 - BMIQ Normalization start ---------------\n\n")

library(wateRmelon)

#load(file="GRset_Noob_filt_3.RData")
load(file="GRset_Noob_filt_1.RData")
load(file="ann850k.RData")
designType=as.character(ann850k$Type)
designType[which(designType=="I")]="1"
designType[which(designType=="II")]="2"
designType=as.integer(designType)
i=match(rownames(GRset_Noob_filt),rownames(ann850k))
rm(ann850k)
designType=designType[i]
beta=getBeta(GRset_Noob_filt)
save(beta,file="beta_Noob_filt_1.RData")
rm(GRset_Noob_filt,i)
fNameOut="_betaBmiq"
#betaBmiq=matrix(nrow=nrow(beta),ncol=ncol(beta),dimnames=list(rownames(beta),colnames(beta)))
#write.table(c("id",rownames(beta)),file=paste(fNameOut,".txt",sep=""), sep="\t", col.names=F, row.names=F, quote=F,append=F)
for (samId in 1:ncol(beta)) {
    #for (samId in 1:2) {
    #res=BMIQ(beta.v=beta[,samId],design.v=designType)
    #res=BMIQ(beta.v=beta[,samId],design.v=designType,plots=TRUE,sampleID=colnames(beta)[samId])
    res=BMIQ(beta.v=beta[,samId],design.v=designType,plots=F,sampleID=colnames(beta)[samId])
    #betaBmiq[,samId]=res$nbeta
    #write.table(c(colnames(beta)[samId],res$nbeta),file=fNameOut, sep="\t", col.names=F, row.names=F, quote=F,append=T)
    x=res$nbeta
    save(x,file=paste(colnames(beta)[samId],fNameOut,".RData",sep=""))
}
fNameOut="_betaBmiq"
#fileList=dir(pattern=paste(fNameOut,".RData",sep=""))
load(file="beta_Noob_filt_1.RData")
betaBmiq=matrix(nrow=nrow(beta),ncol=ncol(beta),dimnames=list(rownames(beta),colnames(beta)))
rm(beta)
for (samId in 1:ncol(betaBmiq)) {
    load(file=paste(colnames(betaBmiq)[samId],fNameOut,".RData",sep=""))
    betaBmiq[,samId]=x
}
save(betaBmiq,file="betaBmiq.RData")
cpgId=rownames(betaBmiq)
fName=paste(dirOut,"betaBmiq",cohortFlag,".txt",sep="")
write.table(paste(c("cpgId",colnames(betaBmiq)),collapse="\t"),file=fName,sep="\n",col.names=F,row.names=F,quote=F,append=F)
k1=seq(1,nrow(betaBmiq),length=1000); k2=k1[2:length(k1)]; k1=c(1,k1[2:(length(k1)-1)]+1)
print(length(k1))
for (k in 1:length(k1)) {
    if (k%%10==0) print(k)
    kk=k1[k]:k2[k]
    write.table(cbind(cpgId[kk],betaBmiq[kk,]),file=fName,sep="\t",col.names=F,row.names=F,quote=F,append=T)
}
rm(cpgId,k,k1,k2,kk,fName)

if (F) {
#tbl=read.table(paste(fNameOut,".txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
betaBmiq=t(as.matrix(tbl[-1]))
colnames(betaBmiq)=tbl$id
rownames(betaBmiq)=names(tbl)
save(betaBmiq,file="betaBmiq.RData")
}

## M-value matrix
load(file="betaBmiq.RData")
offset=0.0000001
mVal=betaBmiq
rm(betaBmiq)
mVal[which(mVal==0)]=offset; mVal[which(mVal==1)]=1-offset
mVal=log2(mVal/(1-mVal))
save(mVal,file="mValBmiq.RData")
sd_mVal=apply(mVal,1,sd,na.rm=T)
save(sd_mVal,file="sdMValBmiq.RData")
rm(mVal,sd_mVal)

# Plot raw, noob and bmiq normalized values in one plot. Mark bad samples
library(RColorBrewer)
#load(file="GRset_Noob_filt_1.RData")
#targets=as.data.frame(pData(GRset_Noob_filt))
#rm(GRset_Noob_filt)
targets=as.data.frame(pData(RGset))
#RGset_NoobBMIQ <- BMIQ(beta.v=GRset_Noob_filt,design.v=designType)
#res=BMIQ(beta.v=beta[i,samId],design.v=designType[i],nL=3,doH=TRUE,nfit=50000,th1.v=c(0.2,0.75),th2.v=NULL,niter=5,tol=0.001,plots=TRUE,sampleID=colnames(beta)[samId],fNameSuf=paste(datType,subsetFNName,subsetName,sep=""))
#Plot raw, Noob Normalized and Noob+BMIQ normalized data on a single plot frame
pdf(paste(dirOut,"densityPlot_bmiq",cohortFlag,".pdf",sep=""))
par(mfrow=c(1,3))
densityPlot(RGset , sampGroups = targets$periDsal , main = "Raw" , legend = FALSE, xlab = "Beta Values")
legend("top",cex=cexLeg,legend=levels(factor(targets$periDsal)) , text.col = brewer.pal(8,"Dark2"))
densityPlot(beta_Noob , sampGroups = targets$periDsal , main = "Noob Normalized" , legend = FALSE, xlab = "Beta Values")
legend("top",cex=cexLeg,legend=levels(factor(targets$periDsal)) , text.col = brewer.pal(8,"Dark2"))
#densityPlot(RGset_NoobBMIQ , sampGroups = targets$periDsal , main = "Noob+BMIQ Normalized" , legend = FALSE, xlab = "Beta Values")
densityPlot(betaBmiq , sampGroups = targets$periDsal , main = "Noob+BMIQ Normalized" , legend = FALSE, xlab = "Beta Values")
legend("top",cex=cexLeg,legend=levels(factor(targets$periDsal)) , text.col = brewer.pal(8,"Dark2"))
dev.off()

# Plot raw, noob and bmiq normalized values in one plot. Mark bad samples
library(RColorBrewer)
datadir="results/"
load(file=paste(datadir,"RGset.RData",sep=""))
targets=as.data.frame(pData(RGset))
rm(RGset)
load(paste(datadir,"stat_sample_betaBmiq_sample.RData",sep=""))
res=res[match(targets$id,rownames(res)),]
samId=which(res[,"perc_50"]>0.4)
pdf(paste(dirOut,"plotMedian_bmiq",cohortFlag,"_QC.pdf",sep=""))
plot(sort(res[samId,"perc_50"]))
dev.off()
pdf(paste(dirOut,"plotIQR_bmiq",cohortFlag,"_QC.pdf",sep=""))
plot(sort(res[samId,"perc_75"]-res[samId,"perc_25"]))
abline(h=c(0.645,0.68),lty="dotted")
dev.off()
load(paste(datadir,"stat_sample_beta_noob_sample.RData",sep=""))
thres=c(0.65,0.60)
grp=rep("0. good",nrow(res)); grp[which(res[,"perc_50"]<thres[1])]="1. low 1"; grp[which(res[,"perc_50"]<thres)]="2. low 2"
grp=paste(grp," from noob",sep="")
grp[-samId]="3. low from bmiq"
pdf(paste(dirOut,"densityPlot_bmiq",cohortFlag,"_QC.pdf",sep=""))
par(mfrow=c(1,3))
if (T) {
load(file=paste(datadir,"beta_raw.RData",sep=""))
beta=beta[,samId]
densityPlot(beta,sampGroups = grp[samId] , main = "Raw" , legend = FALSE, xlab = "Beta Values")
rm(beta)
legend("top",cex=cexLeg,legend=levels(factor(grp)) , text.col = brewer.pal(8,"Dark2"))
load(file=paste(datadir,"beta_Noob.RData",sep=""))
beta_Noob=beta_Noob[,samId]
densityPlot(beta_Noob , sampGroups = grp[samId] , main = "Noob Normalized" , legend = FALSE, xlab = "Beta Values")
rm(beta_Noob)
legend("top",cex=cexLeg,legend=levels(factor(grp)) , text.col = brewer.pal(8,"Dark2"))
}
load(file=paste(datadir,"betaBmiq.RData",sep=""))
betaBmiq=betaBmiq[,samId]
densityPlot(betaBmiq, sampGroups = grp[samId], main = "Noob+BMIQ Normalized" , legend = FALSE, xlab = "Beta Values")
rm(betaBmiq)
legend("top",cex=cexLeg,legend=levels(factor(grp)) , text.col = brewer.pal(8,"Dark2"))
dev.off()

write.table(cbind(targets,quality=grp),file=paste(dirOut,"stat_sample_betaNoobBmiq",cohortFlag,"_QC.txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F,append=F)

## Distribution of beta-values
load(file="betaBmiq.RData")
res=t(apply(betaBmiq,2,function(x) {
    c(quantile(x,probs=c(0,.25,.5,.75,1),na.rm=T),mean(x,na.rm=T))
}))
colnames(res)=c(paste("perc_",c(0,25,50,75,100),sep=""),"mean")
save(res,file="stat_sample_betaBmiq_sample.RData")

load(file="stat_sample_betaBmiq_sample.RData")
clin=clin1[match(rownames(res),clin1$id),]
j=which(!is.na(clin$periDsal))
png("densityPlot_stat_sample_betaBmiq_sample.png")
par(mfrow=c(2,2))
for (k in c("perc_25","perc_50","perc_75","mean")) {
    plot(density(res[j,k],na.rm=T),main=k,xlab="Noob+BMIQ-normalized beta-values")
}
dev.off()
jj=j[which(res[j,"mean"]<.4 | res[j,"mean"]>.6)]
length(jj)
jj=j[which(res[j,"perc_25"]<.1 | res[j,"perc_25"]>.2)]
length(jj)
jj=j[which(res[j,"perc_50"]<.4 | res[j,"perc_50"]>.6)]
length(jj)
jj=j[which(res[j,"perc_75"]<.85 | res[j,"perc_75"]>.95)]
length(jj)

#save.image("tmp4.RData")

## ----------------------------------
# Impute data

library(impute)

datadir="results/"
load(file=paste(datadir,"mValBmiq.RData",sep=""))
mVal=mVal[1:10000,]
mValI=impute.knn(mVal)$data
rm(mVal)
save(mValI,file="mValBmiqImp.RData")
betaBmiqImp=2^mValI
rm(mValI)
betaBmiqImp=betaBmiqImp/(1+betaBmiqImp)
save(betaBmiqImp,file="betaBmiqImp.RData")

## ----------------------------------
#Batch Effect analysis: PCA and ComBAT

#MDS plots to look at clusters - association between samples due to certain variables such as group or Beadchip
library(limma)
library(RColorBrewer)

datadir="results/"

load(file=paste(datadir,"RGset.RData",sep=""))
targets=as.data.frame(pData(RGset))
rm(RGset)
load(file=paste(datadir,"betaBmiqImp.RData",sep=""))
x=betaBmiqImp
x=x[1:10000,]
colnames(x)=as.character(1:ncol(x))

pdf(paste(dirOut,"plotMds_bmiq",cohortFlag,".pdf",sep=""))
#plotMDS(RGset_NoobBMIQ , gene.selection="common" ,col=pal[factor(targets$periDsal)] , main = "PCA plot with Group only",legend = legend("topright", legend=levels(factor(targets$periDsal)), text.col=pal, bg="white", cex=0.7))
#plotMDS(RGset_NoobBMIQ , gene.selection="common" ,col=pal[factor(targets$beadchip)] , main = "PCA plot with Group and beadchip",legend = legend("topright", legend=levels(factor(targets$beadchip)), text.col=pal, bg="white", cex=0.7))
pal <- brewer.pal(8,"Dark2")
plotMDS(x, gene.selection="common" ,col=pal[factor(paste(targets$periDsal,targets$caco))] , main = "PCA plot with Group only",legend = legend("topright", legend=levels(factor(paste(targets$periDsal,targets$caco))), text.col=pal, bg="white", cex=0.7))
pal <- rainbow(sum(!duplicated(targets$beadchip)))
plotMDS(x, gene.selection="common" ,col=pal[factor(targets$beadchip)] , main = "PCA plot with beadchip",legend = legend("topright", legend=levels(factor(targets$beadchip)), text.col=pal, bg="white", cex=0.7))
dev.off()

#PCobj <- prcomp(RGset_NoobBMIQ , retx = T , center = T , scale. = T)
PCobj <- prcomp(x , retx = T , center = T , scale. = T, na.action=na.omit)
save(PCobj,file="PCobj.RData")
attributes(PCobj)
PCs <- PCobj$x
summaryPC <- summary(PCobj)
pdf(paste(dirOut,"screePlot_bmiq",cohortFlag,".pdf",sep=""))
screeplot(PCobj , npcs = 5)
dev.off()
pdf(paste(dirOut,"pcaPlot_bmiq",cohortFlag,".pdf",sep=""))
plot(PCobj , type = "l")
dev.off()

## -------------------
colList=c("brown","red","orange","yellow","green","cyan","skyblue","blue","pink","peachpuff","purple","darkgreen","limegreen","salmon","gray","gold","antiquewhite","steelblue","aquamarine","lightcyan","turquoise","hotpink","black")
colList=c("skyblue","blue","yellow","purple","black","red","orange","green","cyan","darkgreen")
colList2=c("white","black")

fNameOut=paste("_bmiqImp",cohortFlag,sep="")
header="BMIQ normalized, imputed"

varList=c("sex","periDsal","caco")
varName=paste(varList," ",sep="")

annCol=targets
annColAll=annCol

if (is.null(varList)) {
    colCol=NULL
} else {
    colCol=matrix(nrow=length(varList),ncol=nrow(annCol))
    for (varId in 1:length(varList)) {
        if (is.numeric(annColAll[,varList[varId]]) & sum(!duplicated(annColAll[!is.na(annColAll[,varList[varId]]),varList[varId]]))>10) {
            j=match(annCol$id,annColAll$id)
            x=annColAll[,varList[varId]]
            if (varList[varId]%in%c("nanodrop_260.280","nanodrop_260.230")) {
                x=x*10
            } else if (varList[varId]=="pairedEndReads") {
                x=x/10^6
            } else if (varList[varId]%in%c("ruxdex_via","vorinrux_via")) {
                x=(x+6)*10
            }
            x=round(x)
            lim=range(x,na.rm=T)
            x[x<lim[1]]=lim[1]; x[x>lim[2]]=lim[2]
            x=x-min(x,na.rm=T)+1
            grpUniq=lim[1]:lim[2]
            colColUniq=maPalette(high=colList2[2],low=colList2[1],k=length(grpUniq))
            colCol[varId,]=colColUniq[x[j]]
        } else {
            x=as.character(annColAll[,varList[varId]])
            x[x==""]=NA; x=as.integer(as.factor(x))
            grpUniq=sort(unique(x))
            x=x[match(annCol$id,annColAll$id)]
            if (length(grpUniq)<=length(colList)) {
                colCol[varId,]=colList[x]
            } else {
                colCol[varId,]=rainbow(length(grpUniq))[x]
            }
        }
    }
    rownames(colCol)=varName
}

for (varId in 1:length(varList)) {
    pdf(paste(dirOut,"scorePlot_pca_",varList[varId],fNameOut,".pdf",sep=""))
    par(mfcol=c(2,2))
    plot(PCs[,"PC1"],PCs[,"PC2"],main=paste("PCA\n",header,sep=""),xlab="Score: PC1",ylab="Score: PC2",col=colCol[varId,])
    plot(PCs[,"PC1"],PCs[,"PC3"],main=paste("PCA\n",header,sep=""),xlab="Score: PC1",ylab="Score: PC3",col=colCol[varId,])
    plot(PCs[,"PC2"],PCs[,"PC3"],main=paste("PCA\n",header,sep=""),xlab="Score: PC2",ylab="Score: PC3",col=colCol[varId,])
    #sampleColorLegend(tls=colorInfo$grp,col=colorInfo$col,legendTitle=NULL)
    dev.off()
}

## ----------------------------------
## Get cell type estimates

library(minfi)

compCellType="CordBlood"; compCellTypeName="_cordBlood"; cellTypeVec=c("nRBC","CD8T","CD4T","NK","Bcell","Mono","Gran") ## default

cohortFlag="_periDsal"
datadir="results/"

ns.orig <- get(".harmonizeDataFrames", envir = asNamespace("minfi"))
.harmonizeDataFrames.quickfix <- function (x, y) {
    stopifnot(is(x, "DataFrame"))
    stopifnot(is(y, "DataFrame"))
    x.only <- setdiff(names(x), names(y))
    y.only <- setdiff(names(y), names(x))
    if (length(x.only) > 0) {
        df.add <- x[1, x.only]
        is.na(df.add[1, ]) <- TRUE
        y <- cbind(y, df.add)
    }
    if (length(y.only) > 0) {
        df.add <- data.frame(matrix(ncol = length(y.only), nrow = dim(x)[1]))
        names(df.add) <- y.only
        x <- cbind(x, df.add)
    }
    list(x = x, y = y[, names(x)])
}
environment(.harmonizeDataFrames.quickfix) <- environment(ns.orig)
attributes(.harmonizeDataFrames.quickfix) <- attributes(ns.orig)
assignInNamespace(".harmonizeDataFrames", .harmonizeDataFrames.quickfix, ns="minfi")

load(file=paste(datadir,"RGset.RData",sep=""))

if (computerFlag=="cluster") {
    dirClin="/home/royr/project/JoeWiemels/data/periDsal/"
} else {
    dirClin="/Users/royr/UCSF/JoeWiemels/leukMeth/docs/periDsal/"
}
fNameClin="sampleInfo_periDsal_20181119.txt"
phen=read.table(paste(dirClin,fNameClin,sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
phen=phen[which(phen$keep==1),]

j=match(sampleNames(RGset),phen$id); j1=which(!is.na(j)); j2=j[j1]
RGset=RGset[,j1]
phen=phen[j2,]

counts <- estimateCellCounts(RGset,compositeCellType=compCellType,cellTypes=cellTypeVec,meanPlot=FALSE)

save(counts,file=paste("counts_tmp",compCellTypeName,cohortFlag,".RData",sep=""))

write.table(cbind(id=rownames(counts),counts),file=paste("cellCount_minfi",compCellTypeName,cohortFlag,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)

## ------------
cohortFlag="_periDsal"

dirMethRef=""
if (computerFlag=="cluster") {
    dirClin=dirMeth="data/"
    switch(cohortFlag,
    "_periDsal"={
        dirMeth="/cbc2/data1/ritu/wiemelsJ_periDsal/idat/"
        dirClin="/home/royr/project/JoeWiemels/data/periDsal/"
        fNameMeth="betaBmiq"
        fNameClin="sampleInfo_periDsal_20181001.txt"
    }
    )
} else {
    switch(cohortFlag,
    "_periDsal"={
        dirMeth="/Users/royr/UCSF/JoeWiemels/leukMeth/docs/periDsal/idat/"
        dirClin="/Users/royr/UCSF/JoeWiemels/leukMeth/docs/periDsal/"
        fNameMeth="betaBmiq"
        fNameClin="sampleInfo_periDsal_20181001.txt"
    }
    )
}

phen=read.table(paste(dirClin,fNameClin,sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
tbl=read.table(paste("periDsal/stat_sample_betaNoobBmiq_periDsal_QC.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
phen=cbind(phen,quality=tbl$quality[match(phen$id,tbl$id)])
phen=phen[which(phen$keep==1),]
tbl=read.table(paste(dirClin,"cellCount_minfi_cordBlood_periDsal.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
j=match(phen$id,tbl$id); j1=which(!is.na(j)); j2=j[j1]
phen=cbind(phen,tbl[match(phen$id,tbl$id),which(!names(tbl)%in%names(phen))])

phen$caco2=as.character(phen$caco)
phen$caco2[which(phen$caco=="control")]="ctrl"

samId=which(phen$quality!="3. low from bmiq" & !is.na(phen$periDsal) & !is.na(phen$caco))
phen=phen[samId,]

grpUniq=c("nRBC","CD8T","CD4T","NK","Bcell","Mono","Gran")
lim=c(0,.7)
colList=c("red","blue","green","yellow","magenta","cyan","brown")
png(paste("boxplot_cellType",cohortFlag,"_%1d.png",sep=""),width=4*480,height=2*480)
#par(mar=c(5, 4, 4, 2) + 0.1)
par(mar=c(7, 4, 4, 2) + 0.1)
header="PERI/DSAL"
j=1:nrow(phen)
x=c(as.matrix(phen[j,grpUniq]))
grp=grp2=c()
colVec=c()
for (gId in 1:length(grpUniq)) {
    grp=c(grp,paste(grpUniq[gId],"\n",phen$periDsal[j],"/",phen$caco2[j],sep=""))
    grp2=c(grp2,phen$quality[j])
    colVec=c(colVec,rep(colList[gId],4))
}
#grp=rep(grpUniq,each=length(j))
boxplot(x~grp,main=header,ylab="Proportion",ylim=lim,pch=20,cex=0.5,cex.axis=1.5,las=3,col=colVec)
dev.off()

## -----------------------------
library(coin)

phen2=phen
phen2$caco=2-as.integer(as.factor(phen$caco))
phen2$periDsal=2-as.integer(as.factor(phen$periDsal))

varName1=c("periDsal","caco")
n=length(varName1)*length(grpUniq)
tmp=rep(NA,n)
tmpC=rep("",n)
out=data.frame(variable=rep(varName1,each=length(grpUniq)),cellType=rep(grpUniq,length(varName1)),pv=tmp)
for (vId1 in 1:length(varName1)) {
    for (gId in 1:length(grpUniq)) {
        k=which(out$variable==varName1[vId1] & out$cellType==grpUniq[gId])
        out$pv[k]=pvalue(wilcox_test(phen2[,grpUniq[gId]]~as.factor(phen2[,varName1[vId1]]),distribution="exact"))
    }
}

tbl=out
tbl$pv=round(tbl$pv,2)
fName="stat_cellType_periDsal.txt"
write.table(tbl,file=fName, sep="\t", col.names=T, row.names=F, quote=F,append=F)

## ----------------------------------
