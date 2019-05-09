
computerFlag="cluster"
computerFlag=""

if (computerFlag=="cluster") {
} else {
    dirSrc="/Users/royr/UCSF/"
    dirSrc2=dirSrc
    setwd(paste(dirSrc2,"JoeWiemels/leukMeth/epic",sep=""))
}

####################################################################
####################################################################

## ----------------------------------------------
candGeneFlag=list(adjP="nonSnp",gene="")
adjPFlag="qv"
covESFlag=""

transformFlag=""
transformFlag="_mVal"

## ----------------------------------------------
datadir="results/"
load(paste(datadir,"ann850k.RData",sep=""))
ann=as.data.frame(ann850k)
rownames(ann)=NULL
names(ann)[match(c("Name","chr","pos"),names(ann))]=c("IlmnID","CHR","MAPINFO")
ann$CHR[which(ann$CHR=="chrX")]="chr23"
ann$CHR[which(ann$CHR=="chrY")]="chr24"
ann$CHR=as.integer(sub("chr","",ann$CHR))
ann$snp=as.integer(!is.na(ann$Probe_rs))
ann$keep=0; ann$keep[which(ann$snp==0 & ann$CHR%in%1:22)]=1

ann$geneSym=sapply(toupper(ann$UCSC_RefGene_Name),function(x) {
    strsplit(x,";")[[1]][1]
},USE.NAMES=F)
ann$geneSym[is.na(ann$geneSym)]=""


####################################################################
library(data.table)

datadir="results/comparison/"
stat_s_p=read.table(paste(datadir,"stat_season_covPlate_subsetPeriCtrl_periDsal.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
stat_s_sp=read.table(paste(datadir,"stat_season_covSexPlate_subsetPeriCtrl_periDsal.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
stat_s_spc=read.table(paste(datadir,"stat_season_covSexPlateCelltype_subsetPeriCtrl_periDsal.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)

if (F) {
load(paste(datadir,"stat_season_covPlate_subsetPeriCtrl_periDsal.RData",sep=""))
stat_s_p=as.data.frame(all.results)
load(paste(datadir,"stat_season_covSexPlate_subsetPeriCtrl_periDsal.RData",sep=""))
stat_s_sp=as.data.frame(all.results)
load(paste(datadir,"stat_season_covSexPlateCelltype_subsetPeriCtrl_periDsal.RData",sep=""))
stat_s_spc=as.data.frame(all.results)
rm(all.results)
}

if (F) {
    for (k in 1:ncol(x1)) {
        if (any(x1[,k]!=x2[,k],na.rm=T) | any(is.na(x1[,k])!=is.na(x2[,k]))) print(k)
    }
}

datadir="../results/comparison/"
stat_s_p_12=read.table(paste(datadir,"stat_methResp_season_ctrlSubset_covPlate_allGuthSet1Set2_bmiq.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
stat_s_sp_12=read.table(paste(datadir,"stat_methResp_season_ctrlSubset_covSexPlate_allGuthSet1Set2_bmiq.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
stat_s_spc_12=read.table(paste(datadir,"stat_methResp_season_ctrlSubset_covSexPlateCelltype_allGuthSet1Set2_bmiq.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)

## -------------------

datadir="results/comparison/"
stat_cc_p=read.table(paste(datadir,"stat_caco_covSexPlateCelltype_subsetPeri_periDsal.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
stat_cc_d=read.table(paste(datadir,"stat_caco_covSexPlateCelltype_subsetDsal_periDsal.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
stat_dp_co=read.table(paste(datadir,"stat_dsalPeri_covSexPlateCelltype_subsetCtrl_periDsal.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
stat_dp_ca=read.table(paste(datadir,"stat_dsalPeri_covSexPlateCelltype_subsetCase_periDsal.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)

## -------------------

datadir="results/comparison/"
datadir="results/comparison/subset/"
stat_cc2_p_m=read.table(paste(datadir,"stat_methResp_caco_periSubset_covSexPlate_covPrinComp1234_covEpStr_periDsal_bmiq_mVal.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
stat_cc2_d_m=read.table(paste(datadir,"stat_methResp_caco_dsalSubset_covSexPlate_covPrinComp1234_covEpStr_periDsal_bmiq_mVal.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
stat_dp2_co_m=read.table(paste(datadir,"stat_methResp_dsalPeri_ctrlSubset_covSexPlate_covPrinComp1234_covEpStr_periDsal_bmiq_mVal.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
stat_dp2_ca_m=read.table(paste(datadir,"stat_methResp_dsalPeri_caseSubset_covSexPlate_covPrinComp1234_covEpStr_periDsal_bmiq_mVal.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
stat_cc22_p_m=read.table(paste(datadir,"stat_methResp_caco_periSubset_covSexPlate_covPrinComp1234_periDsal_bmiq_mVal.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
stat_cc22_d_m=read.table(paste(datadir,"stat_methResp_caco_dsalSubset_covSexPlate_covPrinComp1234_periDsal_bmiq_mVal.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)

datadir="results/comparison/"
stat_cc2_p=read.table(paste(datadir,"stat_methResp_caco_periSubset_covSexPlate_covPrinComp1234_covEpStr_periDsal_bmiq.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
stat_cc2_d=read.table(paste(datadir,"stat_methResp_caco_dsalSubset_covSexPlate_covPrinComp1234_covEpStr_periDsal_bmiq.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
datadir="results/comparison/subset/"
stat_dp2_co=read.table(paste(datadir,"stat_methResp_dsalPeri_ctrlSubset_covSexPlate_covPrinComp1234_covEpStr_periDsal_bmiq.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
stat_dp2_ca=read.table(paste(datadir,"stat_methResp_dsalPeri_caseSubset_covSexPlate_covPrinComp1234_covEpStr_periDsal_bmiq.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
#stat_cc22_p=read.table(paste(datadir,"stat_methResp_caco_periSubset_covSexPlate_covPrinComp1234_periDsal_bmiq.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
#stat_cc22_d=read.table(paste(datadir,"stat_methResp_caco_dsalSubset_covSexPlate_covPrinComp1234_periDsal_bmiq.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)

####################################################################

library(qvalue)
source(paste(dirSrc,"functions/TTest.9.1.7.R",sep=""))

colIdPV="pv"

compList=c("sea_p","sea_sp","sea_spc")
compList=paste(rep(c("sea_p","sea_sp","sea_spc"),each=2),c("_12",""),sep="")
compList=c("cc_p","cc_d","dp_co","dp_ca")
compList=c("dp2_co")
compList=paste(c("cc2_d","cc22_d"),"_m",sep="")
compList=paste(c("cc2_p","cc2_d","dp2_co","dp2_ca"),"_m",sep="")
compList=paste(rep(c("cc2_p","cc2_d","dp2_co","dp2_ca","cc22_p","cc22_d"),each=2),c("","_m"),sep="")
compList=c("cc2_p","cc2_d","dp2_co","dp2_ca")
for (compId in compList) {
	cat("\n\n==================",compId,"==================\n")
	switch(compId,
        "sea_p"={
            stat2=stat_s_p
            colIdPV=names(stat2)[grep("pv_",names(stat2))]
        },
        "sea_sp"={
            stat2=stat_s_sp
            colIdPV=names(stat2)[grep("pv_",names(stat2))]
        },
        "sea_spc"={
            stat2=stat_s_spc
            colIdPV=names(stat2)[grep("pv_",names(stat2))]
        },
        "sea_p_12"={
            stat2=stat_s_p_12
            colIdPV=names(stat2)[grep("pv_",names(stat2))]
        },
        "sea_sp_12"={
            stat2=stat_s_sp_12
            colIdPV=names(stat2)[grep("pv_",names(stat2))]
        },
        "sea_spc_12"={
            stat2=stat_s_spc_12
            colIdPV=names(stat2)[grep("pv_",names(stat2))]
        },
        "cc_p"={
            stat2=stat_cc_p
            colIdPV=names(stat2)[grep("pv_",names(stat2))]
        },
        "cc_d"={
            stat2=stat_cc_d
            colIdPV=names(stat2)[grep("pv_",names(stat2))]
        },
        "dp_co"={
            stat2=stat_dp_co
            colIdPV=names(stat2)[grep("pv_",names(stat2))]
        },
        "dp_ca"={
            stat2=stat_dp_ca
            colIdPV=names(stat2)[grep("pv_",names(stat2))]
        },
        "cc2_p"={
            stat2=stat_cc2_p
            colIdPV=names(stat2)[grep("pv_",names(stat2))]
        },
        "cc2_d"={
            stat2=stat_cc2_d
            colIdPV=names(stat2)[grep("pv_",names(stat2))]
        },
        "dp2_co"={
            stat2=stat_dp2_co
            colIdPV=names(stat2)[grep("pv_",names(stat2))]
        },
        "dp2_ca"={
            stat2=stat_dp2_ca
            colIdPV=names(stat2)[grep("pv_",names(stat2))]
        },
        "cc22_p"={
            stat2=stat_cc22_p
            colIdPV=names(stat2)[grep("pv_",names(stat2))]
        },
        "cc22_d"={
            stat2=stat_cc22_d
            colIdPV=names(stat2)[grep("pv_",names(stat2))]
        },
        "cc2_p_m"={
            stat2=stat_cc2_p_m
            colIdPV=names(stat2)[grep("pv_",names(stat2))]
        },
        "cc2_d_m"={
            stat2=stat_cc2_d_m
            colIdPV=names(stat2)[grep("pv_",names(stat2))]
        },
        "dp2_co_m"={
            stat2=stat_dp2_co_m
            colIdPV=names(stat2)[grep("pv_",names(stat2))]
        },
        "dp2_ca_m"={
            stat2=stat_dp2_ca_m
            colIdPV=names(stat2)[grep("pv_",names(stat2))]
        },
        "cc22_p_m"={
            stat2=stat_cc22_p_m
            colIdPV=names(stat2)[grep("pv_",names(stat2))]
        },
        "cc22_d_m"={
            stat2=stat_cc22_d_m
            colIdPV=names(stat2)[grep("pv_",names(stat2))]
        }
        )
    colInfo=data.frame(var1=c("probeID",paste(c("coef","se","pv"),"_",rep(paste("season",c("spring","summer","winter"),sep=""),each=3),sep="")),var2=c("cpgId",paste(c("coef","se","pv"),"_",rep(paste("seasonOfBirth",c("Spring","Summer","Winter"),sep=""),each=3),sep="")),stringsAsFactors=F)
    k=match(colInfo$var1,names(stat2)); k1=which(!is.na(k)); k2=k[k1]
    if (length(k1)!=0) names(stat2)[k2]=colInfo$var2[k1]
    k=match(colInfo$var1,colIdPV); k1=which(!is.na(k)); k2=k[k1]
    if (length(k1)!=0) colIdPV[k2]=colInfo$var2[k1]
    iA2=match(stat2$cpgId,ann$IlmnID)
    for (cId in 1:length(colIdPV)) {
        stat2$qv=NA
        if (candGeneFlag$adjP=="nonSnp") {
            i1=list(which(ann$keep[iA2]==1))
        } else {
            i1=list(1:nrow(stat2))
        }
        for (iThis in 1:length(i1)) {
            ii=!is.na(stat2[i1[[iThis]],colIdPV[cId]])
            stat2$qv[i1[[iThis]]][ii]=qvalue(stat2[i1[[iThis]][ii],colIdPV[cId]])$qvalue
            if (F) {
                stat2$qv[i1[[iThis]]]=getAdjustedPvalue(stat2[i1[[iThis]],colIdPV[cId]],method=adjPFlag,strict=T)
                if (all(is.na(stat2$qv[i1[[iThis]]]))) {
                    cat("Using BH method !!!\n")
                    stat2$qv[i1[[iThis]]]=getAdjustedPvalue(stat2[i1[[iThis]],colIdPV[cId]],method=adjPFlag,strict=F)
                }
            }
        }
        i=which(stat2$qv<.05)
        if (any(stat2$qv[i]<stat2[i,colIdPV[cId]])) {cat("Q-value < p-value !!!\n")}
        names(stat2)[which(names(stat2)=="qv")]=sub("pv",adjPFlag,colIdPV[cId])
        
    }
	switch(compId,
        "sea_p"={
            stat_s_p=stat2
        },
        "sea_sp"={
            stat_s_sp=stat2
        },
        "sea_spc"={
            stat_s_spc=stat2
        },
        "sea_p_12"={
            stat_s_p_12=stat2
        },
        "sea_sp_12"={
            stat_s_sp_12=stat2
        },
        "sea_spc_12"={
            stat_s_spc_12=stat2
        },
        "cc_p"={
            stat_cc_p=stat2
        },
        "cc_d"={
            stat_cc_d=stat2
        },
        "dp_co"={
            stat_dp_co=stat2
        },
        "dp_ca"={
            stat_dp_ca=stat2
        },
        "cc2_p"={
            stat_cc2_p=stat2
        },
        "cc2_d"={
            stat_cc2_d=stat2
        },
        "dp2_co"={
            stat_dp2_co=stat2
        },
        "dp2_ca"={
            stat_dp2_ca=stat2
        },
        "cc22_p"={
            stat_cc22_p=stat2
        },
        "cc22_d"={
            stat_cc22_d=stat2
        },
        "cc2_p_m"={
            stat_cc2_p_m=stat2
        },
        "cc2_d_m"={
            stat_cc2_d_m=stat2
        },
        "dp2_co_m"={
            stat_dp2_co_m=stat2
        },
        "dp2_ca_m"={
            stat_dp2_ca_m=stat2
        },
        "cc22_p_m"={
            stat_cc22_p_m=stat2
        },
        "cc22_d_m"={
            stat_cc22_d_m=stat2
        }
        )
}

save.image("tmp2.RData")

## -------------------
stat1=stat_s_p_12
stat2=stat_s_p
i=match(stat1$cpgId,stat2$cpgId); iS1=which(!is.na(i)); iS2=i[iS1]

####################################################################
####################################################################
####################################################################
####################################################################
## ----------------------------------------------

plotFlag=""
plotFlag=c("","manhattanPlot")
plotFlag="_qqPlot"
plotFlag="_manhattanPlot"
plotFlag="_volcanoPlot"
plotFlag="_onePlot"
plotFlag="_agreePlot"

parList=list(ylimM=c(0,3))
parList=list()
if (candGeneFlag$gene%in%c("telaml","telaml38","telaml118")) {
    parList=list(xlimM=c(x[1]-diff(x)/2,x[2]+diff(x)/2))
    grp=candGeneFlag$geneSym
    grpUniq=unique(grp)
    x=c()
    for (gId in 1:length(grpUniq)) {
        x2=range(ann$MAPINFO[which(ann$IlmnID%in%candGeneFlag$IlmnID[which(grp==grpUniq[gId])])],na.rm=T)
        x2=c(x2[1]-diff(x2)/2,x2[2]+diff(x2)/2)
        x=c(x,x2)
    }
    parList=list(chr=c(12,21),xlimM=matrix(x,ncol=2,byrow=T),ylimM=c(0,7),xlabM=grpUniq)
}

geneSumFlag=T
geneSumFlag=F

outlierFlag=F
outlierFlag=T

subsetFFlag=""
subsetFFlag="_autosomesNoChr21"
subsetFFlag="_chr21"

pThresVec=99
pThresVec=c(0.001,0.0001,0.00001)
pThresVec=0.05

for (pThresList in pThresVec) {
    pThres2=pThresList
## -------------------
compList=c("sea_p","sea_sp","sea_spc")
compList=c("cc_p","cc_d","dp_co","dp_ca")
compList=paste(c("cc2_p","cc2_d","dp2_co","dp2_ca"),"_m",sep="")
compList=paste(c("dp2_co"),"_m",sep="")
compList=paste(rep(c("cc22_p","cc22_d"),each=2),c("","_m"),sep="")
compList=c("dp2_co","dp2_ca")
compList=paste(c("cc2_d","cc22_d"),"_m",sep="")
compList=paste(rep(c("cc2_p","cc2_d","dp2_co","dp2_ca","cc22_p","cc22_d"),each=2),c("","_m"),sep="")
compList=c("cc2_p","cc2_d","dp2_co","dp2_ca")
compList=c("dp2_ca")
compList=c("cc2_d")
if (plotFlag[1]%in%c("_qqPlot","_histogram","_volcanoPlot","_manhattanPlot","_agreePlot")) {
    #png(paste(sub("_","",plotFlag[1]),"_%1d.png",sep=""),width=3*240, height=2*240)
    #par(mfcol=c(2,3))
    #png(paste(sub("_","",plotFlag[1]),"_periDsal_pv",formatC(pThres2,format="fg"),".png",sep=""),width=3*240, height=3*240)
    png(paste(sub("_","",plotFlag[1]),"_periDsal_qv",formatC(pThres2,format="fg"),subsetFFlag,".png",sep=""),width=3*240, height=3*240)
    #par(mfcol=c(3,3))
}
for (compId in compList) {
    if (exists("stat1")) rm(stat1)
    switch(compId,
        "sea_p"={
            stat2=stat_s_p; fName0=paste("_covPlate_subsetPeriCtrl_periDsal",sep=""); compName1=paste("Beta-value based\nPeri ctrl: meth ~ varPred\nCov: plate",sep="")
            if (exists("stat_s_p_12")) stat1=stat_s_p_12
        },
        "sea_sp"={
            stat2=stat_s_sp; fName0=paste("_covSexPlate_subsetPeriCtrl_periDsal",sep=""); compName1=paste("Beta-value based\nPeri ctrl: meth ~ varPred\nCov: Sex, plate",sep="")
            if (exists("stat_s_sp_12")) stat1=stat_s_sp_12
        },
        "sea_spc"={
            stat2=stat_s_spc; fName0=paste("_covSexPlateCelltype_subsetPeriCtrl_periDsal",sep=""); compName1=paste("Beta-value based\nPeri ctrl: meth ~ varPred\nCov: Sex, plate, cellType",sep="")
            if (exists("stat_s_spc_12")) stat1=stat_s_spc_12
        },
        "cc_p"={
            stat2=stat_cc_p; fName0=paste("_covSexPlateCelltype_subsetPeri_periDsal",sep=""); compName1=paste("Beta-value based\nPeri: meth ~ caco\nCov: Sex, plate, cellType",sep="")
            if (exists("stat_cc_d")) {
                stat1=stat_cc_d
                cohortName1=c("Dsal","eri")
            }
        },
        "cc_d"={
            stat2=stat_cc_d; fName0=paste("_covSexPlateCelltype_subsetDsal_periDsal",sep=""); compName1=paste("Beta-value based\nDsal: meth ~ caco\nCov: Sex, plate, cellType",sep="")
            if (exists("stat_cc_p")) {
                stat1=stat_cc_p
                cohortName1=c("Peri","Dsal")
            }
        },
        "dp_co"={
            stat2=stat_dp_co; fName0=paste("_covSexPlateCelltype_subsetCtrl_periDsal",sep=""); compName1=paste("Beta-value based\nCtrl: meth ~ dsalPeri\nCov: Sex, plate, cellType",sep="")
            if (exists("stat_dp_ca")) {
                stat1=stat_dp_ca
                cohortName1=c("Case","Ctrl")
            }
        },
        "dp_ca"={
            stat2=stat_dp_ca; fName0=paste("_covSexPlateCelltype_subsetCase_periDsal",sep=""); compName1=paste("Beta-value based\nCase: meth ~ dsalPeri\nCov: Sex, plate, cellType",sep="")
        },
        "cc2_p"={
            stat2=stat_cc2_p; fName0=paste("_periSubset_covSexPlate_covPrinComp1234_covEpStr_periDsal",sep=""); compName1=paste("Beta-value based\nPeri: meth ~ caco\nCov: Sex, plate, ReFACTor, epistructure",sep="")
            if (exists("stat_cc2_d")) {
                stat1=stat_cc2_d
                cohortName1=c("Dsal","Peri")
            }
        },
        "cc2_d"={
            stat2=stat_cc2_d; fName0=paste("_dsalSubset_covSexPlate_covPrinComp1234_covEpStr_periDsal",sep=""); compName1=paste("Beta-value based\nDsal: meth ~ caco\nCov: Sex, plate, ReFACTor, epistructure",sep="")
            if (exists("stat_cc2_p")) {
                stat1=stat_cc2_p
                cohortName1=c("Peri","Dsal")
            }
        },
        "dp2_co"={
            stat2=stat_dp2_co; fName0=paste("_ctrlSubset_covSexPlate_covPrinComp1234_covEpStr_periDsal",sep=""); compName1=paste("Beta-value based\nCtrl: meth ~ dsalPeri\nCov: Sex, plate, ReFACTor, epistructure",sep="")
            if (exists("stat_dp2_ca")) {
                stat1=stat_dp2_ca
                cohortName1=c("Case","Ctrl")
            }
        },
        "dp2_ca"={
            stat2=stat_dp2_ca; fName0=paste("_caseSubset_covSexPlate_covPrinComp1234_covEpStr_periDsal",sep=""); compName1=paste("Beta-value based\nCase: meth ~ dsalPeri\nCov: Sex, plate, ReFACTor, epistructure",sep="")
            if (exists("stat_dp2_co")) {
                stat1=stat_dp2_co
                cohortName1=c("Ctrl","Case")
            }
        },
        "cc22_p"={
            stat2=stat_cc22_p; fName0=paste("_periSubset_covSexPlate_covPrinComp1234_periDsal",sep=""); compName1=paste("Beta-value based\nPeri: meth ~ caco\nCov: Sex, plate, ReFACTor",sep="")
            if (exists("stat_cc22_d")) {
                stat1=stat_cc22_d
                cohortName1=c("Dsal","Peri")
            }
        },
        "cc22_d"={
            stat2=stat_cc22_d; fName0=paste("_dsalSubset_covSexPlate_covPrinComp1234_periDsal",sep=""); compName1=paste("Beta-value based\nDsal: meth ~ caco\nCov: Sex, plate, ReFACTor",sep="")
            if (exists("stat_cc22_p")) {
                stat1=stat_cc22_p
                cohortName1=c("Peri","Dsal")
            }
        },
        "cc2_p_m"={
            stat2=stat_cc2_p_m; fName0=paste("_periSubset_covSexPlate_covPrinComp1234_covEpStr_periDsal_mVal",sep=""); compName1=paste("M-value based\nPeri: meth ~ caco\nCov: Sex, plate, ReFACTor, epistructure",sep="")
            if (exists("stat_cc2_d_m")) {
                stat1=stat_cc2_d_m
                cohortName1=c("Dsal","Peri")
            }
        },
        "cc2_d_m"={
            stat2=stat_cc2_d_m; fName0=paste("_dsalSubset_covSexPlate_covPrinComp1234_covEpStr_periDsal_mVal",sep=""); compName1=paste("M-value based\nDsal: meth ~ caco\nCov: Sex, plate, ReFACTor, epistructure",sep="")
            if (exists("stat_cc2_p_m")) {
                stat1=stat_cc2_p_m
                cohortName1=c("Peri","Dsal")
            }
        },
        "dp2_co_m"={
            stat2=stat_dp2_co_m; fName0=paste("_ctrlSubset_covSexPlate_covPrinComp1234_covEpStr_periDsal_mVal",sep=""); compName1=paste("M-value based\nCtrl: meth ~ dsalPeri\nCov: Sex, plate, ReFACTor, epistructure",sep="")
            if (exists("stat_dp2_ca_m")) {
                stat1=stat_dp2_ca_m
                cohortName1=c("Case","Ctrl")
            }
        },
        "dp2_ca_m"={
            stat2=stat_dp2_ca_m; fName0=paste("_caseSubset_covSexPlate_covPrinComp1234_covEpStr_periDsal_mVal",sep=""); compName1=paste("M-value based\nCase: meth ~ dsalPeri\nCov: Sex, plate, ReFACTor, epistructure",sep="")
        },
        "cc22_p_m"={
            stat2=stat_cc22_p_m; fName0=paste("_periSubset_covSexPlate_covPrinComp1234_periDsal_mVal",sep=""); compName1=paste("M-value based\nPeri: meth ~ caco\nCov: Sex, plate, ReFACTor",sep="")
            if (exists("stat_cc22_d_m")) {
                stat1=stat_cc22_d_m
                cohortName1=c("Dsal","Peri")
            }
        },
        "cc22_d_m"={
            stat2=stat_cc22_d_m; fName0=paste("_dsalSubset_covSexPlate_covPrinComp1234_periDsal_mVal",sep=""); compName1=paste("M-value based\nDsal: meth ~ caco\nCov: Sex, plate, ReFACTor",sep="")
            if (exists("stat_cc22_p_m")) {
                stat1=stat_cc22_p_m
                cohortName1=c("Peri","Dsal")
            }
        }
    )
    fName0=paste(fName0,subsetFFlag,sep="")
    x=strsplit(compId,"_")[[1]]
    predVarFlag=x[1]
    switch(predVarFlag,
        "sea"={
            if (exists("stat1")) cohortName1=c("Set1+set2 ctrl","Peri ctrl")
        },
        "cc"={
            if (exists("stat1")) cohortName1=c("Peri","Dsal")
        },
        "dp"={
            if (exists("stat1")) cohortName1=c("Case","Ctrl")
        }
    )
    colListPV=names(stat2)[grep("pv_",names(stat2))]
    for (cId in 1:length(colListPV)) {
        if (length(colListPV)==1) {
            compName=compName1
            fName1=fName0
            colIdEst="coef"; colIdPV=c("pv",adjPFlag); pThres=0.05
            colIdEst="coef"; colIdPV=c("pv","pv"); pThres=99
            colIdEst="coef"; colIdPV=c("pv","pv"); pThres=0.001
            
            colIdEst=names(stat2)[grep("coef_",names(stat2))][cId]
            x=sub("coef_","",colIdEst)
            compName=sub("varPred",x,compName1)
            fName1=paste("_",x,fName0,sep="")
            colIdPV=c(colListPV[cId],sub("pv",adjPFlag,colListPV[cId])); pThres=0.05
            colIdPV=rep(sub("pv",adjPFlag,colListPV[cId]),2); pThres=0.05
            #colIdPV=rep(colListPV[cId],2); pThres=0.05
        } else {
            colIdEst=names(stat2)[grep("coef_",names(stat2))][cId]
            x=sub("coef_","",colIdEst)
            compName=sub("varPred",x,compName1)
            fName1=paste("_",x,fName0,sep="")
            #colIdPV=c(colListPV[cId],colListPV[cId]); pThres=0.05
            colIdPV=c(colListPV[cId],sub("pv",adjPFlag,colListPV[cId])); pThres=0.05
            #colIdPV=c(colListPV[cId],sub("pv",adjPFlag,colListPV[cId])); pThres=0.001
        }
        switch(predVarFlag,
            "sea"={
                if (exists("stat1")) cohortName=c(sub("Peri ctrl:","",compName),cohortName1)
            },
            "cc"={
                if (exists("stat1")) cohortName=c(sub("Dsal:","",compName),cohortName1)
            },
            "dp"={
                if (exists("stat1")) cohortName=c(sub("Ctrl:","",compName),cohortName1)
            },
            "cc2"={
                if (exists("stat1")) cohortName=c(sub("Dsal:","",compName),cohortName1)
            },
            "dp2"={
                if (exists("stat1")) cohortName=c(sub("Ctrl:|Case:","",compName),cohortName1)
            }
        )
            
        }
        pThresName=ifelse(pThres>1,1,pThres)
        nm=c()
        if (!colListPV[cId]%in%names(stat2)) {cat(colListPV[cId]," not a column in stat file!!!\n",sep=""); next}
        if (length(nm)==0) nm=names(stat2)
        i=match(stat2$cpgId,ann$IlmnID)
        if (candGeneFlag$gene%in%c("mgmt","mgmt13","telaml","telaml38","telaml118")) {
            compName=sub("Coefficient:","Coef:",sub("M-value based",paste(toupper(candGeneFlag$gene),". Mval based",sep=""),compName))
            fName=paste("_",candGeneFlag$gene,fName1,sep="")
            stat2=stat2[which(tolower(ann[i,candGeneFlag$variable])%in%candGeneFlag[[candGeneFlag$variable]]),]
        } else {
            fName=fName1
            stat2=stat2[which(ann$keep[i]==1),]
        }
        iA2=match(stat2$cpgId,ann$IlmnID)
        
        ann2=ann[iA2,]
        
        ####################################################################
        
        if (all(is.na(stat2[,colIdEst]))) {cat("No stats computed !!!\n"); next}
        
        if (candGeneFlag$gene=="") {
            i=which(ann$keep[iA2]==1)
            switch(subsetFFlag,
                "_chr21"={
                    compName[1]=paste("Chr 21. ",compName[1],sep="")
                    i=i[which(ann$CHR[iA2][i]%in%21)]
                },
                "_autosomesNoChr21"={
                    compName[1]=paste("Autosomes, no chr 21. ",compName[1],sep="")
                    i=i[which(ann$CHR[iA2][i]%in%c(1:20,22))]
                }
            )
        } else {
            i=1:nrow(stat2)
        }
        stat=stat2
        if (exists("stat1")) {
            stat11=stat1[match(stat$cpgId,stat1$cpgId),]
            rm(stat1)
        }
        
        cat("\n\n",compName,"\n",sep="")
        if (length(colListPV)!=1) {
            fName=paste(fName1,sub("pv_","_",colIdPV[1][grep("pv_",colIdPV[1])]),sep="")
        }
        
        x1=matrix(0,nrow=2,ncol=2,dimnames=list(c("dn","up"),paste(colIdPV[2],c(">=","<"),pThresName,sep="")))
        ii=i[which(stat[i,colIdEst]!=0)]
        x2=table(stat[ii,colIdEst]>0,stat[ii,colIdPV[2]]<pThres)
        x1[match(rownames(x2),c("FALSE","TRUE")),match(colnames(x2),c("FALSE","TRUE"))]=x2
        print(x1)
        
        if (plotFlag[1]%in%c("_agreePlot")) {
            x1=matrix(0,nrow=2,ncol=2,dimnames=list(paste(colIdPV[2],c(">=","<"),pThresName,sep=""),paste(colIdPV[2],c(">=","<"),pThresName,sep="")))
            ii=i[which(stat[i,colIdEst]!=0)]
            x2=table(stat11[ii,colIdPV[2]]<pThres,stat[ii,colIdPV[2]]<pThres)
            x1[match(rownames(x2),c("FALSE","TRUE")),match(colnames(x2),c("FALSE","TRUE"))]=x2
            #x1=as.table(x1,dnn=cohortName[2:3])
            x=paste(rep(" ",nchar(paste(cohortName[3],": ",sep=""))),collapse="")
            rownames(x1)[1]=paste(cohortName[2],": ",rownames(x1)[1],sep="")
            rownames(x1)[2]=paste(x,rownames(x1)[2],sep="")
            colnames(x1)[1]=paste(cohortName[3],": ",colnames(x1)[1],sep="")
            #x1=rbind(c(cohortName[3],""),x1)
            #x1=cbind(c(cohortName[2],"",""),x1)
            print(x1)
        }
        
        if (plotFlag[1]=="_onePlot") {
            png(paste("plots",fName,".png",sep=""),width=3*240, height=1*240)
            par(mfcol=c(1,3))
        }
        
        if (plotFlag[1]%in%c("","_onePlot","_qqPlot")) {
            if (plotFlag[1]=="") {
                png(paste("qqPlot",fName,".png",sep=""))
                header=compName
            } else if (plotFlag[1]=="_qqPlot") {
                header=compName
            } else {
                header=""
            }
            pvs <- sort(na.exclude(stat[i,colIdPV[1]]))
            qqplot(-log10(runif(length(pvs),0,1)),-log10(pvs),xlab="Expected -log10(p-values) by random",ylab="Observed -log10(p-values)",pch=19,cex.axis=1.5,cex.lab=1.5,main=header)
            abline(0,1)
            if (plotFlag[1]=="") {
                dev.off()
            }
        }
        
        if (plotFlag[1]%in%c("","_onePlot","_histogram")) {
            if (plotFlag[1]=="") {
                png(paste("histogram",fName,".png",sep=""))
                header=compName
            } else if (plotFlag[1]=="_histogram") {
                header=compName
            } else {
                header=""
            }
            hist(stat[i,colIdPV[1]],xlab="P-value",pch=19,cex.axis=1.5,cex.lab=1.5,main=header)
            if (plotFlag[1]=="") {
                dev.off()
            } else if (plotFlag[1]=="_onePlot") {
                title(main=compName)
            }
        }
        
        if (plotFlag[1]%in%c("","_onePlot","_volcanoPlot")) {
            if (plotFlag[1]=="") {
                png(paste("volcanoPlot",fName,"_",colIdPV[2],pThresName,".png",sep=""))
                header=compName
            } else if (plotFlag[1]=="_volcanoPlot") {
                header=compName
            } else {
                header=""
            }
            iThis=i
            if (!outlierFlag) {
                x=quantile(abs(stat[iThis,colIdEst]),probs=.95,na.rm=T)
                iThis=iThis[which(abs(stat[iThis,colIdEst])<=x)]
            }
            if ("ylimM"%in%names(parList)) yLim=parList$ylimM else yLim=NULL
            plot(stat[iThis,colIdEst],-log10(stat[iThis,colIdPV[1]]),ylim=yLim,xlab="Estimate",ylab="-log10(p-value)",pch=19,cex.axis=1.5,cex.lab=1.5,main=header,col="grey")
            ii=iThis[which(stat[iThis,colIdPV[2]]<pThres)]
            points(stat[ii,colIdEst],-log10(stat[ii,colIdPV[1]]),pch=19,col="red")
            if (plotFlag[1]=="") {
                dev.off()
            }
        }
        
        if (plotFlag[1]%in%c("_manhattanPlot")) {
            if (plotFlag[1]=="") {
                png(paste("manhattanPlot",fName,"_",colIdPV[2],pThresName,".png",sep=""))
                header=compName
            } else if (plotFlag[1]=="_manhattanPlot") {
                header=compName
            } else {
                header=""
            }
            iThis=i
            if ("ylimM"%in%names(parList)) yLim=parList$ylimM else yLim=NULL
            xLab=""; xLim=NULL
            if ("xlabM"%in%names(parList)) {xLab=parList$xlabM; xLim=parList$xlimM}
            for (p in 1:length(xLab)) {
                #plot(ann[iA2[iThis],"MAPINFO"],-log10(stat[iThis,colIdPV[1]]),xlab=paste(xLab[p],": Chromosome ",ann[iA2[iThis][1],"CHR"],sep=""),xlim=xLim[p,],ylim=yLim,ylab="-log10(p-value)",pch=19,cex.axis=1.5,cex.lab=1.5,main=header,col="grey")
                iThis=i[which(ann$CHR[iA2[i]]==parList$chr[p] & ann$MAPINFO[iA2[i]]>=parList$xlimM[p,1] & ann$MAPINFO[iA2[i]]<=parList$xlimM[p,2])]
                plot(ann[iA2[iThis],"MAPINFO"],-log10(stat[iThis,colIdPV[1]]),xlab=paste(xLab[p],": Chromosome ",ann[iA2[iThis][1],"CHR"],sep=""),xlim=xLim[p,],ylim=yLim,ylab="-log10(p-value)",pch=19,cex.axis=1.5,cex.lab=1.5,main=header,col="grey")
                ii=iThis[which(stat[iThis,colIdPV[2]]<pThres)]
                points(ann[iA2[ii],"MAPINFO"],-log10(stat[ii,colIdPV[1]]),pch=19,col="red")
            }
            if (plotFlag[1]=="") {
                dev.off()
            }
        }
        
        if (plotFlag[1]%in%c("_agreePlot")) {
            if (plotFlag[1]=="") {
                png(paste("agreePlot",fName,"_",colIdPV[2],pThresName,".png",sep=""))
                header=cohortName[1]
            } else if (plotFlag[1]=="_agreePlot") {
                header=cohortName[1]
            } else {
                header=""
            }
            iThis=i
            iThis=iThis[which(stat[iThis,colIdPV[1]]<pThres2)]
            iThis=iThis[which(stat11[iThis,colIdPV[1]]<pThres2)]
            if (length(iThis)==0) {
                iThis=i
                plot(stat11[iThis,colIdEst],stat[iThis,colIdEst],type="n")
            } else {
                header=cohortName[1]
                #header=paste(sub("Beta-value based\n ","",header),"\n",cohortName[3]," ",strsplit(colIdPV[1],"_")[[1]][1],"<",formatC(pThres2,format="fg"),", N ",length(iThis),", prop agree ",round(mean(sign(stat11[iThis,colIdEst])==sign(stat[iThis,colIdEst]),na.rm=T),2),sep="")
                header=paste(sub("Beta-value based\n ","",header),"\n",strsplit(colIdPV[1],"_")[[1]][1],"<",formatC(pThres2,format="fg")," in both, N ",length(iThis),", prop agree ",round(mean(sign(stat11[iThis,colIdEst])==sign(stat[iThis,colIdEst]),na.rm=T),2),sep="")
                lim=range(stat[iThis,colIdEst],stat11[iThis,colIdEst],na.rm=T)
                plot(stat11[iThis,colIdEst],stat[iThis,colIdEst],xlim=lim,ylim=lim,xlab=paste(cohortName[2],": Estimate",sep=""),ylab=paste(cohortName[3],": Estimate",sep=""),main=header,pch=19,cex=0.5)
                abline(c(0,1),lty="dotted")
                abline(c(0,-1),lty="dotted")
                abline(h=0,lty="dotted")
                abline(v=0,lty="dotted")
            }
            if (plotFlag[1]=="") {
                dev.off()
            }
        }
        
        #plot(stat2$coef,-log10(stat2$pv)); iii=which(stat2$qv<.05); points(stat2$coef[iii],-log10(stat2$pv)[iii],col="red")
        
        if (plotFlag[1]=="_onePlot") {
            #par(mfrow=c(1,1))
            #title(main=compName)
            dev.off()
        }
        
        if (F) {
            png(paste("plots",fName,".png",sep=""),width=3*240, height=1*240)
            par(mfcol=c(1,3))
            plot(1:6)
            plot(1:6)
            title(main=sub(" ReFACTor","\nReFACTor",compName))
            plot(1:6)
            dev.off()
        }
        
        ii=i[order(stat[i,colIdPV[2]])]
        ii=ii[which(stat[ii,colIdPV[2]]<pThres)]
        if (length(ii)!=0) {
            #tbl=cbind(ann[iA2,][ii,],stat[ii,c(colIdEst,unique(colIdPV))])
            colId=c("IlmnID","CHR","MAPINFO","UCSC_RefGene_Name","Relation_to_Island","snp","Probe_rs")
            tbl=cbind(ann[iA2,colId][ii,],stat[ii,which(!nm%in%c("cpgId","gene_genotype","cpgId_ahrr"))])
            names(tbl)=c(colId,nm[which(!nm%in%c("cpgId","gene_genotype","cpgId_ahrr"))])
            if ("gene_genotype"%in%names(stat)) tbl=cbind(gene_genotype=stat$gene_genotype[ii],tbl)
            if ("cpgId_ahrr"%in%names(stat)) {
                nm2=c(paste(colId,"_ahrr"),names(tbl))
                tbl=cbind(ann[match(stat$cpgId_ahrr[ii],ann$IlmnID),colId],tbl)
                names(tbl)=nm2
            }
            if (pThresName<1) fName2=paste("_",colIdPV[2],pThresName,sep="") else fName2=""
            write.table(tbl, file=paste("stat",fName1,fName2,".txt",sep=""), append=F,col.names=T,row.names=F, sep="\t",quote=F)
        }
        
        ####################################################################
        ## Gene level summarization of p-values
        ####################################################################
        if (geneSumFlag) {
            pvFlag=""
            
            k=which(colnames(stat2)==colIdPV[2])
            kk=grep("Mean",colnames(stat2)[k])
            if (length(kk)!=0) k=k[-kk]
            
            fName2=fName
            # pThres=.001; prId=which(ann2$CHR%in%1:22 & !ann2$IlmnID%in%snpVec & !is.na(stat2[,k]) & stat2[,k]<pThres)
            # prId=which(ann2$CHR%in%1:22 & !ann2$IlmnID%in%snpVec & !is.na(stat2[,k]) & stat2[,k]<pThres)
            if (candGeneFlag$gene=="") {
                prId=which(ann2$keep==1 & !is.na(stat2[,k]) & stat2[,k]<pThres)
            } else {
                prId=which(!is.na(stat2[,k]) & stat2[,k]<pThres)
            }
            
            if (length(prId)!=0) {
                annotSel=ann2[prId,]
                pval=stat2[,k][prId]
                coef=stat2[,colIdEst[1]][prId]
                signGiven=stat2$sLEU[prId]
                
                
                if (F) {
                    # Assign polycomb status to each CpG
                    tmp1 = strsplit(annotSel$UCSC_RefGene_Accession,";")
                    names(tmp1)=paste(1:length(tmp1),":",sep="")
                    tmp2 = unlist(tmp1)
                    tmp3 = as.numeric(gsub("[:][[:digit:]]*$","",names(tmp2)))
                    tmp4 = data.frame(UCSC_REFGENE_ACCESSION=tmp2,
                    rowid=tmp3, stringsAsFactors=FALSE)
                    load("PolycombComplete-120109.RData")
                    tmp5 = merge(polycombTab,tmp4)
                    tmp6 = unique(tmp5$rowid)
                    annotSel$PcG = rep(0, dim(annotSel)[1])
                    annotSel$PcG[tmp6] = 1
                }
                
                # Get an index of CpGs by gene region
                tmp1 = strsplit(annotSel$UCSC_RefGene_Name,";")
                names(tmp1)=paste(1:length(tmp1),":",sep="")
                tmp2 = unlist(tmp1)
                tmp3 = as.numeric(gsub("[:][[:digit:]]*$","",names(tmp2)))
                tmp4 = strsplit(annotSel$UCSC_RefGene_Group,";")
                names(tmp4)=paste(1:length(tmp4),":",sep="")
                tmp5 = unlist(tmp4)
                tmp6 = as.numeric(gsub("[:][[:digit:]]*$","",names(tmp5)))
                all(tmp3==tmp6)
                
                if (F) {
                    library(qvalue)
                    qval = qvalue(pval)
                    qval$pi0
                    qThresh = max(pval[qval$qv<=0.05])
                    qThresh
                }
                
                GeneAnnotation = unique(data.frame(UCSC_REFGENE_NAME=tmp2,UCSC_REFGENE_GROUP=tmp5,rowid=tmp3, stringsAsFactors=FALSE))
                GeneIndex = split(GeneAnnotation$rowid,with(GeneAnnotation,paste(UCSC_REFGENE_NAME,UCSC_REFGENE_GROUP,sep=":")))
                GeneIndexN = sapply(GeneIndex, length)
                
                if (length(GeneIndexN)!=0) {
                    medPval = sapply(GeneIndex, function(u) median(pval[u],na.rm=T))
                    propHit = sapply(GeneIndex, function(u) mean(pval[u]<pThres,na.rm=T))
                    
                    #isPcG = sapply(GeneIndex, function(u) min(annotSel$PcG[u]))
                    #isNearPcG = sapply(GeneIndex, function(u) max(annotSel$PcG[u]))
                    
                    tmp1 = sapply(GeneIndex, function(u) u[which.min(annotSel$MAPINFO[u])])
                    tmp2 = sapply(GeneIndex, function(u) u[which.max(annotSel$MAPINFO[u])])
                    GeneSym1 = annotSel$UCSC_RefGene_Name[tmp1]
                    GeneSym2 = annotSel$UCSC_RefGene_Name[tmp2]
                    
                    annotSelMap <- as.numeric(annotSel$MAPINFO)
                    tmpF <- function(u) {
                        sgn = sign(coef[u])
                        pv = pval[u]
                        sgnChar = ifelse(pv>pThres, ".", ifelse(sgn<0,"-","+"))
                        paste(sgnChar[order(annotSel$CHR[u], annotSelMap[u])],collapse="")
                    }
                    #Check
                    #if (length(GeneIndex)>2) print(tmpF(GeneIndex[[3]]))
                    
                    hypohyper =sapply(GeneIndex, tmpF)
                    
                    combineSigns <- function(u) {
                        sgn = signGiven[u]
                        sgnChar = ifelse(sgn==0, ".", ifelse(sgn<0,"-","+"))
                        paste(sgnChar[order(annotSel$CHR[u], annotSelMap[u])],collapse="")
                    }
                    hypohyper2 =sapply(GeneIndex, combineSigns)
                    
                    tmp = strsplit(names(GeneIndexN),":")
                    #GeneResults = data.frame(Gene=sapply(tmp,function(u)u[1]),Region=sapply(tmp,function(u)u[2]), nCpG=GeneIndexN, Sign=hypohyper, SignLeuk=hypohyper2, medPval, propHit, stringsAsFactors=FALSE, GeneSymsFirst=GeneSym1, GeneSymsLast=GeneSym2)
                    GeneResults = data.frame(Gene=sapply(tmp,function(u)u[1]),Region=sapply(tmp,function(u)u[2]), nCpG=GeneIndexN, Sign=hypohyper, medPval, propHit, stringsAsFactors=FALSE, GeneSymsFirst=GeneSym1, GeneSymsLast=GeneSym2)
                    
                    rownames(GeneResults)=NULL
                    ord = order(medPval)
                    
                    if (pThres>1) {
                        fName=paste("geneSummary_top500_",ifelse(pvFlag=="",colIdPV[2],"pvPerm"),fName2,".txt",sep="")
                        write.table(GeneResults[ord[1:500],which(!names(GeneResults)%in%c("isPcG","isNearPcG"))], file=fName, append=FALSE,col.names=T,row.names=FALSE, sep="\t",quote=FALSE)
                        
                        fName=paste("geneSummary_med",ifelse(pvFlag=="",colIdPV[2],"PvPerm"),".001",fName2,".txt",sep="")
                        write.table(GeneResults[ord[medPval[ord]<.001],which(!names(GeneResults)%in%c("isPcG","isNearPcG"))], file=fName, append=FALSE,col.names=T,row.names=FALSE, sep="\t",quote=FALSE)
                        
                        fName=paste("geneSummary_med",ifelse(pvFlag=="",colIdPV[2],"PvPerm"),".05",fName2,".txt",sep="")
                        write.table(GeneResults[ord[medPval[ord]<.05],which(!names(GeneResults)%in%c("isPcG","isNearPcG"))], file=fName, append=FALSE,col.names=T,row.names=FALSE, sep="\t",quote=FALSE)
                    } else {
                        fName=paste("geneSummary",fName2,"_",ifelse(pvFlag=="",colIdPV[2],"pvPerm"),pThres,".txt",sep="")
                        write.table(GeneResults[ord,which(!names(GeneResults)%in%c("isPcG","isNearPcG"))], file=fName, append=FALSE,col.names=T,row.names=FALSE, sep="\t",quote=FALSE)
                    }
                    
                    if (F) {
                        tbl=GeneResults[ord,which(!names(GeneResults)%in%c("isPcG","isNearPcG"))]
                        
                        i=prId[grep("FAM5C",ann2$UCSC_RefGene_Name[prId])]
                        i=i[grep("5'UTR",ann2$UCSC_RefGene_Group[i])]
                        tbl[which(tbl$Gene=="FAM5C"),]
                        stat2[i,]
                    }
                }
                
            }
        }
        ####################################################################
    }
}
if (plotFlag[1]%in%c("_qqPlot","_histogram","_volcanoPlot","_manhattanPlot","_agreePlot")) {
    dev.off()
}

####################################################################
####################################################################
####################################################################
####################################################################
