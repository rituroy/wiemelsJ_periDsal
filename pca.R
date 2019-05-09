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

cat("\n================== pca.R ===============\n\n")

dirDat="results/"
phen=read.table(paste(dirDat,"clin_periDsal_20180801.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
load(paste(dirDat,"ann850k.RData",sep=""))
datObj=list(ann=as.data.frame(ann850k),phen=phen)
load(paste(dirDat,"betaSnp.RData",sep=""))
tmp=datObj$ann[1:nrow(betaSnp),]
for (k in 1:ncol(tmp)) {
    tmp[,k]=NA
}
rownames(tmp)=tmp$Name=rownames(betaSnp)
rm(ann850,phen,betaSnp)

####################################################


phen_0=datObj$phen
cutoff=3

########################################################################
########################################################################
########################################################################

outFormat="pdf"
outFormat="png"

absFlag=T
absFlag=F

nClust=c(2,2)

distMethod="euclidean"; absFlag=F
linkMethod="ward.D2"

datList="_betaSnp"

colHmap=c("blue", "red", "white")
colHmap=c("red", "blue", "white")
colList=c("brown","red","orange","yellow","green","cyan","skyblue","blue","pink","peachpuff","purple","darkgreen","limegreen","salmon","gray","gold","antiquewhite","steelblue","aquamarine","lightcyan","turquoise","hotpink","black")
colList=c("skyblue","blue","yellow","purple","black","red","orange","green","cyan","darkgreen")
colList2=c("white","black")

colVec=rep("black",nrow(datObj$phen))
colVec[which(datObj$phen$periDsal=="peri")]="skyblue"
colVec[which(datObj$phen$periDsal=="dsal")]="blue"

centerFlag=""
centerFlag="_mean"
centerFlag="_med"

scaleFlag="_sd" ## mad = 0 for some cases so not used
scaleFlag=""

clustFlag=c("","cluster")
clustFlag=c("cluster","cluster")

subsetList=c("")
subsetList=c("_periDsal")

subsetFList=""; clustFlag=c("cluster","cluster")

annColAll=datObj$phen[,c("id","sex","periDsal","caco")]
annRowAll=datObj$ann

subDir="heatmap"
if (!file.exists(subDir)) dir.create(file.path(subDir))

subDir="legend"
if (!file.exists(subDir)) dir.create(file.path(subDir))
dirL=paste(subDir,"/",sep="")

for (datFlag in datList) {
    subDir=paste("heatmap/",sub("_","",datFlag),sep="")
    if (!file.exists(subDir)) dir.create(file.path(subDir))
    
    header1=""
    switch(datFlag,
        "_betaSnp"={
            header1="SNPs (beta-values): "
            load(paste(dirDat,"betaSnp.RData",sep=""))
            varPred=betaSnp
            rm(betaSnp)
            offset=0.0000001
            varPred[which(varPred==0)]=offset; varPred[which(varPred==1)]=1-offset
            varPred=log2(varPred/(1-varPred))
        }
    )
    annCol=annColAll[match(colnames(varPred),annColAll$id),]
    annRow=annRowAll[match(rownames(varPred),annRowAll$Name),]
    for (subsetFFlag in subsetFList) {
        for (subsetFlag in subsetList) {
            fNameOut=paste("_periDsal",subsetFFlag,subsetFlag,datFlag,sep="")
            header=paste(header1,sep="")
            cat("\n\n",fNameOut,"\n\n")
            prId=1:nrow(varPred)
            samId=which(datObj$phen$sexMatched=="yes")
            if (length(prId)<2 | length(samId)<2) next
            
            x=strsplit(subsetFFlag,"_")[[1]]
            if (any(substr(x,1,2)=="sd")) {
                k=which(substr(x,1,2)=="sd")
                prId=prId[which(sdCntVec[prId]>as.numeric(sub("sd","",x[k])))]
                header=paste(header," ",sub("_","",x[k]),sep="")
            }
            if (any(substr(x,1,3)=="top")) {
                k=which(substr(x,1,3)=="top")
                prId=prId[order(sdCntVec[prId],decreasing=T)][1:as.integer(sub("top","",x[k]))]
                header=paste(header," ",sub("_","",x[k]),sep="")
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
            }
            if (length(samId)==0) {
                cat("No samples !!!\n")
                next
            }
            
            annCol=annCol[samId,]
            annRow=annRowAll[prId,]
            arrayData=varPred[prId,samId]
            
            varList=varName=NULL
            varList=c("sex","periDsal","caco")
            varName=paste(varList," ",sep="")
            
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
            
            ## --------------------
            
            fit=prcomp(t(arrayData), center=T, scale=T)
            
            png(paste("screePlot_pca",fNameOut,".png",sep=""))
            #screeplot(fit,npcs=length(fit$sdev),main=paste("PCA screeplot\n",header,sep=""))
            screeplot(fit,npcs=20,main=paste("PCA screeplot\n",header,sep=""))
            #if (plotCutoffFlag) abline(v=cutoff+1,col="red",lty="dotted")
            dev.off()
            for (varId in 1:length(varList)) {
                png(paste("scorePlot_pca_",varList[varId],fNameOut,".png",sep=""),width=1.5*480,height=1.5*480)
                par(mfcol=c(2,2))
                plot(fit$x[,"PC1"],fit$x[,"PC2"],main=paste("PCA\n",header,sep=""),xlab="Score: PC1",ylab="Score: PC2",col=colCol[varId,])
                plot(fit$x[,"PC1"],fit$x[,"PC3"],main=paste("PCA\n",header,sep=""),xlab="Score: PC1",ylab="Score: PC3",col=colCol[varId,])
                plot(fit$x[,"PC2"],fit$x[,"PC3"],main=paste("PCA\n",header,sep=""),xlab="Score: PC2",ylab="Score: PC3",col=colCol[varId,])
                #sampleColorLegend(tls=colorInfo$grp,col=colorInfo$col,legendTitle=NULL)
                dev.off()
            }
            
            tbl=cbind(id=rownames(fit$x),as.data.frame(fit$x[,1:cutoff],stringsAsFactors=F))
            names(tbl)=c("id",paste("pc",1:cutoff,sep=""))
            rownames(tbl)=NULL
            pcaTbl=tbl
            
            nm=c(names(tbl),names(datObj$phen))
            tbl=cbind(tbl,datObj$phen[samId,])
            tbl=tbl[,!duplicated(nm)]
            write.table(tbl, paste("prinComp",fNameOut,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
            
            pcaTbl$pcOutlier=as.integer(pcaTbl$pc1>=400)
            
            ## --------------------
 
        }
    }
}
