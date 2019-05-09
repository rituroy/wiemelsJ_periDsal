## Wiemels peri/dsal project

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

cat("\n================== heatmap.R ===============\n\n")

dirDat="results/"
phen=read.table(paste(dirDat,"clin_periDsal_20180801.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
tbl=read.table(paste("periDsal/stat_sample_betaNoobBmiq_periDsal_QC.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
phen$quality=tbl$quality[match(phen$id,tbl$id)]
load(paste(dirDat,"ann850k.RData",sep=""))
datObj=list(ann=as.data.frame(ann850k),phen=phen)
datObj$ann$chr[which(datObj$ann$chr=="chrX")]="chr23"
datObj$ann$chr[which(datObj$ann$chr=="chrY")]="chr24"
datObj$ann$chr=as.integer(sub("chr","",datObj$ann$chr))
load(paste(dirDat,"betaSnp.RData",sep=""))
tmp=datObj$ann[1:nrow(betaSnp),]
for (k in 1:ncol(tmp)) {
    tmp[,k]=NA
}
rownames(tmp)=tmp$Name=rownames(betaSnp)
rm(ann850,phen,tbl,betaSnp)

####################################################


########################################################################
########################################################################
########################################################################

library(marray)
source(paste(dirSrc,"functions/heatmap.5.7.R",sep=""))
source(paste(dirSrc,"functions/heatmapAcgh.7.4.R",sep=""))
source(paste(dirSrc,"functions/miscFuncs.1.3.R",sep=""))

outFormat="pdf"
outFormat="png"

absFlag=T
absFlag=F

nClust=c(2,2)
nClust=c(2,15)

distMethod="euclidean"; absFlag=F
linkMethod="ward.D2"

datList="_betaSnp"
datList="_mValNoob"
datList="_mValNoobFilt"
datList="_mValNoobFiltSnp"
datList="_mValNoobFiltNoSnp"
datList=c("_mValNoobFiltSnp","_mValNoobFiltNoSnp")
datList=paste(c("_mValNoobFiltSnp","_mValNoobFiltNoSnp"),rep(c("Chr21","NoChr21"),each=2),sep="")
datList=paste(c("_mValBmiqSnp","_mValBmiqNoSnp"),rep(c("Chr21","NoChr21"),each=2),sep="")
datList="_mValBmiqNoSnpChr21"

colHmap=c("blue", "red", "white")
colHmap=c("red", "blue", "white")
colList=c("brown","red","orange","yellow","green","cyan","skyblue","blue","pink","peachpuff","purple","darkgreen","limegreen","salmon","gray","gold","antiquewhite","steelblue","aquamarine","lightcyan","turquoise","hotpink","black")
colList=c("skyblue","blue","yellow","purple","black","red","orange","green","cyan","darkgreen")
colList2=c("white","black")

centerFlag=""
centerFlag="_mean"
centerFlag="_med"

scaleFlag="_sd" ## mad = 0 for some cases so not used
scaleFlag=""

subsetList=c("")
subsetList=c("_periDsal")

subsetFList=""; clustFlag=c("cluster","cluster")
subsetFList="_top2000"; clustFlag=c("cluster","cluster")

clustFlag=c("cluster","cluster")
clustFlag=c("","cluster")

load(paste(dirDat,"predictedSex.RData",sep=""))
annColAll=cbind(datObj$phen[,c("id","sex","periDsal","caco","quality")],predictedSex=predictedSex[match(datObj$phen$id,names(predictedSex))])
rm(predictedSex)
annRowAll=cbind(datObj$ann,sdCnt=rep(NA,nrow(datObj$ann)))

subDir="heatmap"
if (!file.exists(subDir)) dir.create(file.path(subDir))

subDir="legend"
if (!file.exists(subDir)) dir.create(file.path(subDir))
dirL=paste(subDir,"/",sep="")

for (datFlag in datList) {
    subDir=paste("heatmap/",sub("_","",datFlag),sep="")
    if (!file.exists(subDir)) dir.create(file.path(subDir))
    
    header1=""
    sdCntVec=NULL
    switch(datFlag,
        "_betaSnp"={
            header1="SNP (raw M-value): "
            load(paste(dirDat,"betaSnp.RData",sep=""))
            varPred=betaSnp
            rm(betaSnp)
            offset=0.0000001
            varPred[which(varPred==0)]=offset; varPred[which(varPred==1)]=1-offset
            varPred=log2(varPred/(1-varPred))
        },
        "_mValNoob"={
            header1="Noob-normalized M-value: "
            load(paste(dirDat,"mVal_noob.RData",sep=""))
            annRow=annRowAll[match(rownames(mVal),annRowAll$Name),]
            prId=which(!rownames(mVal)%in%annRow$Name[which(annRow$chr%in%1:22)])
            varPred=mVal[prId,]
            load(paste(dirDat,"sd_mVal_noob.RData",sep=""))
            sdCntVec=sd_mVal[prId]
            rm(mVal,sd_mVal)
        },
        "_mValNoobFilt"={
            header1="Noob-normalized filtered M-value (autosomes only): "
            load(paste(dirDat,"mVal_RGset_noob_filt_1.RData",sep=""))
            annRow=annRowAll[match(rownames(mVal),annRowAll$Name),]
            prId=which(!rownames(mVal)%in%annRow$Name[which(annRow$chr%in%1:22)])
            varPred=mVal[prId,]
            load(paste(dirDat,"sd_mVal_RGset_noob_filt_1.RData",sep=""))
            sdCntVec=sd_mVal[prId]
            rm(mVal,sd_mVal)
        },
        "_mValNoobFiltSnp"={
            header1="Noob-normalized filtered M-value (autosomal SNPs only): "
            load(paste(dirDat,"mVal_RGset_noob_filt_1.RData",sep=""))
            annRow=annRowAll[match(rownames(mVal),annRowAll$Name),]
            prId=which(!is.na(annRow$Probe_rs) & annRow$chr%in%1:22)
            varPred=mVal[prId,]
            load(paste(dirDat,"sd_mVal_RGset_noob_filt_1.RData",sep=""))
            sdCntVec=sd_mVal[prId]
            rm(mVal,sd_mVal)
        },
        "_mValNoobFiltNoSnp"={
            header1="Noob-normalized filtered M-value (no SNPs, autosomes only): "
            load(paste(dirDat,"mVal_RGset_noob_filt_1.RData",sep=""))
            annRow=annRowAll[match(rownames(mVal),annRowAll$Name),]
            prId=which(is.na(annRow$Probe_rs) & annRow$chr%in%1:22)
            varPred=mVal[prId,]
            load(paste(dirDat,"sd_mVal_RGset_noob_filt_1.RData",sep=""))
            sdCntVec=sd_mVal[prId]
            rm(mVal,sd_mVal)
        },
        "_mValNoobFiltSnpChr21"={
            header1="Noob-normalized filtered M-value (chr 21 SNPs only): "
            load(paste(dirDat,"mVal_RGset_noob_filt_1.RData",sep=""))
            annRow=annRowAll[match(rownames(mVal),annRowAll$Name),]
            prId=which(!is.na(annRow$Probe_rs) & annRow$chr%in%21)
            varPred=mVal[prId,]
            load(paste(dirDat,"sd_mVal_RGset_noob_filt_1.RData",sep=""))
            sdCntVec=sd_mVal[prId]
            rm(mVal,sd_mVal)
        },
        "_mValNoobFiltNoSnpChr21"={
            header1="Noob-normalized filtered M-value (no SNPs, chr 21 only): "
            load(paste(dirDat,"mVal_RGset_noob_filt_1.RData",sep=""))
            annRow=annRowAll[match(rownames(mVal),annRowAll$Name),]
            prId=which(is.na(annRow$Probe_rs) & annRow$chr%in%21)
            varPred=mVal[prId,]
            load(paste(dirDat,"sd_mVal_RGset_noob_filt_1.RData",sep=""))
            sdCntVec=sd_mVal[prId]
            rm(mVal,sd_mVal)
        },
        "_mValNoobFiltSnpNoChr21"={
            header1="Noob-normalized filtered M-value (autosomal (no chr 21) SNPs only): "
            load(paste(dirDat,"mVal_RGset_noob_filt_1.RData",sep=""))
            annRow=annRowAll[match(rownames(mVal),annRowAll$Name),]
            prId=which(!is.na(annRow$Probe_rs) & annRow$chr%in%c(1:20,22))
            varPred=mVal[prId,]
            load(paste(dirDat,"sd_mVal_RGset_noob_filt_1.RData",sep=""))
            sdCntVec=sd_mVal[prId]
            rm(mVal,sd_mVal)
        },
        "_mValNoobFiltNoSnpNoChr21"={
            header1="Noob-normalized filtered M-value (no SNPs, no chr 21 autosomes only): "
            load(paste(dirDat,"mVal_RGset_noob_filt_1.RData",sep=""))
            annRow=annRowAll[match(rownames(mVal),annRowAll$Name),]
            prId=which(is.na(annRow$Probe_rs) & annRow$chr%in%c(1:20,22))
            varPred=mVal[prId,]
            load(paste(dirDat,"sd_mVal_RGset_noob_filt_1.RData",sep=""))
            sdCntVec=sd_mVal[prId]
            rm(mVal,sd_mVal)
        },
        "_mValBmiqSnpChr21"={
            header1="BMIQ-normalized M-value (chr 21 SNPs only): "
            load(paste(dirDat,"mvalBmiq.RData",sep=""))
            annRow=annRowAll[match(rownames(mVal),annRowAll$Name),]
            prId=which(!is.na(annRow$Probe_rs) & annRow$chr%in%21)
            varPred=mVal[prId,]
            load(paste(dirDat,"sdMValBmiq.RData",sep=""))
            sdCntVec=sd_mVal[prId]
            rm(mVal,sd_mVal)
        },
        "_mValBmiqNoSnpChr21"={
            header1="BMIQ-normalized M-value (no SNPs, chr 21 only): "
            load(paste(dirDat,"mvalBmiq.RData",sep=""))
            annRow=annRowAll[match(rownames(mVal),annRowAll$Name),]
            prId=which(is.na(annRow$Probe_rs) & annRow$chr%in%21)
            varPred=mVal[prId,]
            load(paste(dirDat,"sdMValBmiq.RData",sep=""))
            sdCntVec=sd_mVal[prId]
            rm(mVal,sd_mVal)
        },
        "_mValBmiqSnpNoChr21"={
            header1="BMIQ-normalized M-value (no chr 21 autosomal (no chr 21) SNPs only): "
            load(paste(dirDat,"mvalBmiq.RData",sep=""))
            annRow=annRowAll[match(rownames(mVal),annRowAll$Name),]
            prId=which(!is.na(annRow$Probe_rs) & annRow$chr%in%c(1:20,22))
            varPred=mVal[prId,]
            load(paste(dirDat,"sdMValBmiq.RData",sep=""))
            sdCntVec=sd_mVal[prId]
            rm(mVal,sd_mVal)
        },
        "_mValBmiqNoSnpNoChr21"={
            header1="BMIQ-normalized M-value (no SNPs,  no chr 21 autosomes only): "
            load(paste(dirDat,"mvalBmiq.RData",sep=""))
            annRow=annRowAll[match(rownames(mVal),annRowAll$Name),]
            prId=which(is.na(annRow$Probe_rs) & annRow$chr%in%c(1:20,22))
            varPred=mVal[prId,]
            load(paste(dirDat,"sdMValBmiq.RData",sep=""))
            sdCntVec=sd_mVal[prId]
            rm(mVal,sd_mVal)
        }
    )
    if (is.null(sdCntVec)) sdCntVec=rep(NA,nrow(varPred))
    sdCntVecAll=sdCntVec
    for (subsetFFlag in subsetFList) {
        for (subsetFlag in subsetList) {
            fName=paste("_periDsal",subsetFFlag,subsetFlag,datFlag,sep="")
            header=paste(header1,sep="")
            cat("\n\n",fName,"\n\n")

            i=match(rownames(varPred),annRowAll$Name)
            annCol=annColAll[match(colnames(varPred),annColAll$id),]
            annRow=annRowAll[i,]
            sdCntVec=sdCntVecAll[i]
            
            #prId=1:nrow(varPred)
            prId=order(annRow$chr*10^9+annRow$pos)
            #samId=which(datObj$phen$sexMatched=="yes")
            samId=which(tolower(annColAll$sex)==tolower(annColAll$predictedSex))
            if (length(prId)<2 | length(samId)<2) next
            
            if (F) {
                x=strsplit(datFlag,"_")[[1]]
                if (any(x=="chr21")) {
                    x=as.integer(sub("chr","",x[k]))
                    prId=prId[which(annRow$chr[prId]==x)]
                    header=paste(header," chr",x,sep="")
                }
                if (any(x=="noChr21")) {
                    x=as.integer(sub("noChr","",x[k]))
                    prId=prId[which(annRow$chr[prId]!=x)]
                    header=paste(header," no chr",x,sep="")
                }
            }
            x=strsplit(subsetFFlag,"_")[[1]]
            if (any(substr(x,1,2)=="sd")) {
                k=which(substr(x,1,2)=="sd")
                prId=prId[which(sdCntVec[prId]>as.numeric(sub("sd","",x[k])))]
                header=paste(header," ",sub("_","",x[k]),sep="")
            }
            if (any(substr(x,1,3)=="top")) {
                k=which(substr(x,1,3)=="top")
                y=min(as.integer(sub("top","",x[k])),length(prId))
                header=paste(header," ",ifelse(length(prId)>y,sub("_","",x[k]),y),sep="")
                prId=prId[order(sdCntVec[prId],decreasing=T)][1:y]
            }
            if (substr(subsetFFlag,1,nchar("_rnd"))=="_rnd") {
                x=as.integer(sub("_rnd","",subsetFFlag))
                set.seed(5393)
                prId=sample(prId,size=x,replace=F)
                header=paste(header," random ",x," genes",sep="")
            }
            if (length(prId)==0) {
                cat("No genes !!!\n")
                next
            }

            if (subsetFlag=="_periDsal") {
                samId=samId[which(!is.na(annCol$periDsal[samId]))]
                header=paste(header,", peri/dsal",sep="")
            } else {
                x=strsplit(subsetFlag,"_")
            }
            if (length(samId)==0) {
                cat("No samples !!!\n")
                next
            }
            
            varList=varName=NULL
            varList=c("quality","sex","periDsal","caco")
            varName=paste(varList," ",sep="")

            varFList=varFName=NULL

            varListAll=varList; varNameAll=varName
            varFListAll=varFList; varFNameAll=varFName

            annCol=annCol[samId,]
            annRow=annRow[prId,]
            arrayData=varPred[prId,samId]
            sdCntVec=sdCntVec[prId]

            if (centerFlag=="_med") {
                arrayData=arrayData-apply(arrayData,1,median,na.rm=T)
            }
            if (centerFlag=="_mean") {
                arrayData=arrayData-apply(arrayData,1,mean,na.rm=T)
            }
            if (scaleFlag=="_sd") {
                arrayData=arrayData/apply(arrayData,1,sd,na.rm=T)
            }

            nameRow=annRow$geneId
            nameRow=rep("",length(prId))
            if (is.null(varFList)) {
                rowCol=NULL
            } else {
                rowCol=matrix(nrow=length(varFList),ncol=nrow(annRow))
                for (varId in 1:length(varFList)) {
                    if (sum(!duplicated(annRowAll[!is.na(annRowAll[,varFList[varId]]),varFList[varId]]))>10) {
                        if (varFList[varId]=="slope") {
                            lim=100*(limSl+1)
                            x=round(100*(annRowAll[,varFList[varId]]+1))
                            x[x<lim[1]]=lim[1]; x[x>lim[2]]=lim[2]
                        } else {
                            x=round(annRowAll[,varFList[varId]])+1
                            x2=as.integer(as.factor(x))
                            j=match(annRow$id,annRowAll$id); j1=which(!is.na(j)); j2=j[j1]
                            lim=range(x,na.rm=T)
                        }
                        grpUniq=lim[1]:lim[2]
                        rowColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                        rowCol[varId,j1]=rowColUniq[x2[j2]]
                    } else {
                        x=as.character(annRowAll[,varFList[varId]])
                        x[x==""]=NA; x=as.integer(as.factor(x))
                        grpUniq=sort(unique(x))
                        x=x[match(annRow$id,annRowAll$id)]
                        if (length(grpUniq)<=length(colList)) {
                            rowCol[varId,]=colList[x]
                        } else {
                            rowCol[varId,]=rainbow(length(grpUniq))[x]
                        }
                    }
                }
                rownames(rowCol)=varFName
            }
            lineClone=seq(5,nrow(arrayData),by=5)
            lineClone=NULL

            nameCol=rep("",length(samId))
            if (length(samId)<=10) nameCol=annCol$id[samId]
            #nameCol=annCol$id[samId]
            lineSam=NULL

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

            if (clustFlag[1]=="cluster") {
                switch(distMethod,
                    "pearson" = {
                        x=cor(t(arrayData),method=distMethod,use="complete.obs")
                        if (absFlag) x=abs(x)
                        distMat <- as.dist(1-x)
                        clustR <- hclust(distMat, method=linkMethod)},
                    "spearman" = {
                        x=cor(t(arrayData),method=distMethod,use="complete.obs")
                        if (absFlag) x=abs(x)
                        distMat <- as.dist(1-x)
                        clustR <- hclust(distMat, method=linkMethod)},
                    "kappa" = {
                        distMat <- getKappaDist(t(arrayData),absolute=absFlag)
                        clustR <- hclust(distMat, method=linkMethod)},
                    "euclidean" = {
                        distMat <- dist(arrayData, method=distMethod)
                        clustR <- hclust(distMat, method=linkMethod)}
                )
            } else {
                clustR=NA
                nClust[1]=NA
            }
            
            if (clustFlag[2]=="cluster") {
                switch(distMethod,
                    "pearson" = {
                        x=cor(arrayData,method=distMethod,use="complete.obs")
                        if (absFlag) x=abs(x)
                        distMat <- as.dist(1-x)
                        clustC <- hclust(distMat, method=linkMethod)},
                    "spearman" = {
                        x=cor(arrayData,method=distMethod,use="complete.obs")
                        if (absFlag) x=abs(x)
                        distMat <- as.dist(1-x)
                        clustC <- hclust(distMat, method=linkMethod)},
                    "kappa" = {
                        distMat <- getKappaDist(arrayData,absolute=absFlag)
                        clustC <- hclust(distMat, method=linkMethod)},
                    "euclidean" = {
                        distMat <- dist(t(arrayData), method=distMethod)
                        clustC <- hclust(distMat, method=linkMethod)}
                )
            } else {
                clustC=NA
                nClust[2]=NA
            }

            print("dim(arrayData)")
            print(dim(arrayData))
            print("summary(c(arrayData))")
            print(summary(c(arrayData)))

            summary(c(arrayData))
            dat2=arrayData
            limit=c(-1,1)
            limit=c(-0.5,0.5)
            limit=c(-0.25,0.25)
            limit=c(-0.2,0.2)

            margins=c(4,.2)
            margins=c(4,4)
            margins=c(10,1)
            main=header

            if (clustFlag[1]=="cluster") {
                while (F) {
                    x=table(cutree(clustR,k=nClust[1]))
                    if(max(x)<5) {
                        nClust[1]=NA
                        break
                    }
                    if (sum(x>=5)==3) break
                    nClust[1]=nClust[1]+1
                }
            }
            if (clustFlag[2]=="cluster") {
                while (F) {
                    x=table(cutree(clustC,k=nClust[2]))
                    if (sum(x>=5)==3) break
                    nClust[2]=nClust[2]+1
                }
            }

            dirH=paste("heatmap/",sep="")
            subDir=paste(dirH,sub("_","",datFlag),"/",ifelse(subsetFlag=="","allSamples",sub("_","",subsetFlag)),sep="")
            if (!file.exists(subDir)) dir.create(file.path(subDir))
            dirH=paste(subDir,"/",sep="")
            if (centerFlag!="_med" & scaleFlag!="") {
                subDir=paste(dirH,sub("_","",paste(centerFlag,ifelse(centerFlag=="","","Centered"),scaleFlag,ifelse(centerFlag=="","","Scaled"),sep="")),sep="")
                if (!file.exists(subDir)) dir.create(file.path(subDir))
                dirH=paste(subDir,"/",sep="")
            }
            if (outFormat=="png") {
                png(paste(dirH,"heatmap",fName,".png",sep=""),width=480*2,height=480*2); cexRow=cexCol=2
            } else {
                pdf(paste(dirH,"heatmap",fName,".pdf",sep="")); cexRow=cexCol=1
            }
            par(cex.main=.7)
            hcc <- heatmap3(x=dat2, Rowv=clustR, Colv=clustC, distfun=distMethod, hclustfun=hclust, symm=F, ColSideColors=colCol, RowSideColors=rowCol, labCol=nameCol, labRow=nameRow, lineRow=lineClone, ncr=nClust[1], ncc=nClust[2], scale="none", na.rm=F, margins=margins, main=main, xlab=NULL, ylab=NULL, zlm=limit, cexCol=cexCol, cexRow=cexRow, high=colHmap[1], low=colHmap[2], mid=colHmap[3])
            dev.off()

            fName2=""

            if (!is.null(colCol)) {
                for (varId in 1:length(varListAll)) {
                    if (sum(!duplicated(annColAll[!is.na(annColAll[,varListAll[varId]]),varListAll[varId]]))>10) {
                        if (outFormat=="png") {
                            png(paste(dirL,"heatmapSampleColorBarLegend_",varListAll[varId],fName2,".png",sep=""),width=480,height=140)
                        } else {
                            pdf(paste(dirL,"heatmapSampleColorBarLegend_",varListAll[varId],fName2,".pdf",sep=""))
                        }
                        x=annColAll[,varListAll[varId]]
                        if (varListAll[varId]%in%c("nanodrop_260.280","nanodrop_260.230")) {
                            x=x*10
                        } else if (varListAll[varId]=="pairedEndReads") {
                            x=x/10^6
                        } else if (varListAll[varId]%in%c("ruxdex_via","vorinrux_via")) {
                            x=(x+6)*10
                        }
                        x=round(x)
                        header=varNameAll[varId]
                        lim=range(x,na.rm=T)
                        if (varListAll[varId]%in%c("nanodrop_260.280","nanodrop_260.230")) {
                            lim=lim/10
                        } else if (varListAll[varId]=="pairedEndReads") {
                            header=paste(header," (in million)",sep="")
                        } else if (varListAll[varId]%in%c("ruxdex_via","vorinrux_via")) {
                            lim=(lim/10)-6
                        }
                        heatmapColorBar(limit=lim,cols=rev(colList2),main=header)
                        if (F) {
                            grpUniq=lim[1]:lim[2]
                            colColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                            heatmapColorBar(limit=lim,cols=c(colColUniq[c(length(colColUniq),1)],median(1:length(colColUniq))),main=varNameAll[varId])
                        }
                    } else {
                        x=as.character(annColAll[,varListAll[varId]]); x[x==""]=NA
                        x2=as.character(annCol[,varList[varId]]); x2[x2==""]=NA
                        if (all(is.na(x2))) next
                        x2=table(x2)
                        x2=names(x2)
                        grp=table(x)
                        grpUniq=names(grp)
                        k=1:length(grpUniq)
                        ttl=grpUniq[k]
                        if (varListAll[varId]%in%c("il7","stat5_bcl2")) {
                            ttl=sapply(ttl,function(x) strsplit(x,": ")[[1]][2],USE.NAMES=F)
                        }
                        #ttl=paste(grpUniq[k]," (",grp[k],")",sep="")
                        k=match(x2,grpUniq)
                        if (!is.null(varListAll) && varListAll[varId]%in%c("sd","mutPerc")) {
                            width = 480; height = 140
                        } else {
                            if (length(grpUniq)<6) {
                                width = 480; height = 480
                            } else {
                                width = 560; height = 960
                            }
                        }
                        if (outFormat=="png") {
                            png(paste(dirL,"heatmapSampleColorBarLegend_",varListAll[varId],fName2,".png",sep=""),width=width,height=height)
                        } else {
                            pdf(paste(dirL,"heatmapSampleColorBarLegend_",varListAll[varId],fName2,".pdf",sep=""))
                        }
                        cexThis=1.5
                        if (outFormat=="pdf" & (length(grpUniq)>15 | max(nchar(grpUniq))>20)) cexThis=1
                        if (outFormat=="pdf") cexThis=1
                        if (length(grpUniq)<=length(colList)) {
                            sampleColorLegend(tls=ttl[k],col=colList[k],legendTitle=varNameAll[varId],cex=cexThis)
                        } else {
                            sampleColorLegend(tls=ttl[k],col=rainbow(length(grpUniq))[k],legendTitle=varNameAll[varId],cex=cexThis)
                        }
                    }
                    dev.off()
                }
            }
            if (!is.null(rowCol)) {
                for (varId in 1:length(varFListAll)) {
                    if (is.numeric(annColAll[,varList[varId]]) & sum(!duplicated(annColAll[!is.na(annColAll[,varList[varId]]),varList[varId]]))>10) {
                        if (outFormat=="png") {
                            png(paste(dirL,"heatmapFeatureColorBarLegend_",varFListAll[varId],fName2,".png",sep=""),width=480,height=140)
                        } else {
                            pdf(paste(dirL,"heatmapFeatureColorBarLegend_",varFListAll[varId],fName2,".pdf",sep=""))
                        }
                        x=round(annRowAll[,varFListAll[varId]])+1
                        lim=range(x,na.rm=T)
                        grpUniq=lim[1]:lim[2]
                        rowColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                        lim=lim-1
                        heatmapColorBar(limit=lim,cols=c(rowColUniq[c(length(rowColUniq),1,median(1:length(rowColUniq)))]))
                    } else {
                        x=as.character(annRowAll[,varFListAll[varId]]); x[x==""]=NA
                        x2=as.character(annRow[,varFList[varId]]); x2[x2==""]=NA
                        if (all(is.na(x2))) next
                        x2=table(x2)
                        x2=names(x2)
                        grp=table(x)
                        grpUniq=names(grp)
                        k=1:length(grpUniq)
                        ttl=grpUniq[k]
                        ttl=paste(grpUniq[k]," (",grp[k],")",sep="")
                        k=match(x2,grpUniq)
                        if (!is.null(varFListAll) && varFListAll[varId]%in%c("sd","mutPerc")) {
                            width = 480; height = 140
                        } else {
                            if (length(grpUniq)<6) {
                                width = 480; height = 480
                            } else {
                                width = 560; height = 960
                            }
                        }
                        if (outFormat=="png") {
                            png(paste(dirL,"heatmapFeatureColorBarLegend_",varFListAll[varId],fName2,".png",sep=""),width=width,height=height)
                        } else {
                            pdf(paste(dirL,"heatmapFeatureColorBarLegend_",varFListAll[varId],fName2,".pdf",sep=""))
                        }
                        cexThis=NULL
                        if (outFormat=="pdf" & (length(grpUniq)>15 | max(nchar(grpUniq))>20)) cexThis=1
                        if (outFormat=="pdf") cexThis=1
                        if (length(grpUniq)<=length(colList)) {
                            sampleColorLegend(tls=ttl[k],col=colList[k],legendTitle=varFNameAll[varId],cex=cexThis)
                        } else {
                            sampleColorLegend(tls=ttl[k],col=rainbow(length(grpUniq))[k],legendTitle=varFNameAll[varId],cex=cexThis)
                        }
                    }
                    dev.off()
                }
            }

            if (is.na(nClust[1])) {
                tbl=cbind(annRow,order=1:nrow(annRow))
            } else {
                n=max(c(nClust[1],5))-1
                tbl=matrix(nrow=nrow(annRow),ncol=n)
                colnames(tbl)=paste("clustId_",1:ncol(tbl)+1,sep="")
                for (kk in 1:ncol(tbl)) {
                    clustId=cutree(clustR,k=kk+1)[clustR$order]
                    k1=which(!duplicated(clustId))
                    for (k in 1:length(k1)) {
                        clustId[which(clustId==clustId[k1[k]])]=paste("cluster",k,sep="")
                    }
                    tbl[,kk]=clustId
                    
                }
                tbl=cbind(annRow[clustR$order,],tbl,order=1:nrow(annRow))
                if (F) {
                    clustId=cutree(clustR,k=nClust[1])[clustR$order]
                    k1=which(!duplicated(clustId))
                    for (k in 1:length(k1)) {
                        clustId[which(clustId==clustId[k1[k]])]=paste("cluster",k,sep="")
                    }
                    tbl=cbind(annRow[clustR$order,],clustId,order=1:nrow(annRow))
                }
            }
            write.table(tbl, paste(dirH,"clusterInfoFeature",fName,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)

            if (is.na(nClust[2])) {
                tbl=cbind(annCol,order=1:nrow(annCol))
            } else {
                n=max(c(nClust[2],5))-1
                tbl=matrix(nrow=nrow(annCol),ncol=n)
                colnames(tbl)=paste("clustId_",1:ncol(tbl)+1,sep="")
                for (kk in 1:ncol(tbl)) {
                    clustId=cutree(clustC,k=kk+1)[clustC$order]
                    k1=which(!duplicated(clustId))
                    for (k in 1:length(k1)) {
                        clustId[which(clustId==clustId[k1[k]])]=paste("cluster",k,sep="")
                    }
                    tbl[,kk]=clustId
                    
                }
                tbl=cbind(annCol[clustC$order,],tbl,order=1:nrow(annCol))
                if (F) {
                    clustId=cutree(clustC,k=nClust[2])[clustC$order]
                    k1=which(!duplicated(clustId))
                    for (k in 1:length(k1)) {
                        clustId[which(clustId==clustId[k1[k]])]=paste("cluster",k,sep="")
                    }
                    tbl=cbind(annCol[clustC$order,],clustId,order=1:nrow(annCol))
                }
            }
            write.table(tbl, paste(dirH,"clusterInfoSample",fName,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
        }
    }
}
if (outFormat=="png") {
    png(paste(dirL,"heatmapColorRange.png",sep=""),width=480,height=140)
} else {
    pdf(paste(dirL,"heatmapColorRange.pdf",sep=""))
}
heatmapColorBar(limit=limit,cols=colHmap)
dev.off()
if (outFormat=="png") {
    png(paste(dirL,"heatmapSampleColorBarLegend",fName2,".png",sep=""))
} else {
    pdf(paste(dirL,"heatmapSampleColorBarLegend",fName2,".pdf",sep=""))
}
cexThis=1.5
ttl=c("low","high")
if (length(grpUniq)<=length(colList)) {
    sampleColorLegend(tls=ttl,col=colList2,legendTitle=NULL,cex=cexThis)
} else {
    sampleColorLegend(tls=ttl,col=colList2,legendTitle=NULL,cex=cexThis)
}
dev.off()
