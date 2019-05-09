dirSrc="/Users/royr/UCSF/"
dirSrc2=dirSrc
setwd(paste(dirSrc2,"JoeWiemels/leukMeth/normalize",sep=""))

##############################################
## Cell type

dirClin="/Users/royr/UCSF/JoeWiemels/leukMeth/docs/periDsal/"
fNameClin="sampleInfo_periDsal_20181119.txt"

phen=read.table(paste(dirClin,fNameClin,sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
phen=phen[which(phen$keep==1),]

varName1=c("nRBC","CD8T","CD4T","NK","Bcell","Mono","Gran")
x=as.matrix(phen[,varName1])
offset=min(c(x)[which(c(x)!=0)])/10
phen$CD4TByCD8T=(phen$CD4T+offset)/(phen$CD8T+offset)
phen$nlr=phen$Gran/(phen$CD4T+phen$CD8T+phen$NK+phen$Bcell)

##############################################
samId=which(!is.na(phen$periDsal))
dat=phen[samId,]
varName1=c("nRBC","CD8T","CD4T","NK","Bcell","Mono","Gran","CD4TByCD8T","nlr")
for (vId1 in 1:length(varName1)) {
    j=1:nrow(dat)
    if (varName1[vId1]=="CD4TByCD8T") j=which(dat[,varName1[vId1]]>=0)
    png(paste("densityPlot_",varName1[vId1],".png",sep=""))
    plot(density(dat[j,varName1[vId1]],na.rm=T),main=varName1[vId1])
    dev.off()
}

##############################################

samId=which(!is.na(phen$periDsal))
dat=phen[samId,]
dat$periDsal=as.integer(dat$periDsal=="dsal")
dat$caco=as.integer(dat$caco=="case")

varName1=c("nRBC","CD8T","CD4T","NK","Bcell","Mono","Gran","CD4TByCD8T","nlr")
varInfo=data.frame(variable=c("periDsal","caco"),subset=c("_ctrl","_peri"),stringsAsFactors=F)
n=length(varName1)*nrow(varInfo)
n=length(varName1)*nrow(varInfo)+length(varName1)
tmp=rep(NA,n); tmpC=rep("",n)
out=data.frame(subset=tmpC,model=tmpC,test=tmpC,est=tmp,se=tmp,stat=tmp,pv=tmp,stringsAsFactors=F)
k=1
for (vId2 in 1:nrow(varInfo)) {
    for (vId1 in 1:length(varName1)) {
        modelThis=paste(varInfo$variable[vId2],"~",varName1[vId1],"+ sex")
        #modelThis=paste(varInfo$variable[vId2],"~",varName1[vId1],"+ sex + gestWeek")
        j=1:nrow(dat)
        switch(varInfo$subset[vId1],
            "_peri"={
                j=which(dat$periDsal==0)
            },
            "_ctrl"={
                j=which(dat$caco==0)
            }
        )
        res=glm(as.formula(modelThis),family="binomial",data=dat[j,])
        res=summary(res)$coef
        out$subset[k]=sub("_","",varInfo$subset[vId2])
        out$model[k]=modelThis
        out$test[k]="logistic regression"
        out$est[k]=res[2,"Estimate"]
        out$se[k]=res[2,"Std. Error"]
        out$stat[k]=res[2,"z value"]
        out$pv[k]=res[2,"Pr(>|z|)"]
        k=k+1
    }
}
for (vId1 in 1:length(varName1)) {
    modelThis=paste(varName1[vId1],"~ periDsal*caco + sex")
    #modelThis=paste(varName1[vId1],"~ periDsal*caco + sex + gestWeek")
    j=1:nrow(dat)
    res=lm(as.formula(modelThis),data=dat[j,])
    res=summary(res)$coef
    p=nrow(res)
    out$model[k]=modelThis
    out$test[k]="logistic regression"
    out$est[k]=res[p,"Estimate"]
    out$se[k]=res[p,"Std. Error"]
    out$stat[k]=res[p,"t value"]
    out$pv[k]=res[p,"Pr(>|t|)"]
    k=k+1
}

tbl=out
for (k in c("est","se","stat")) {
    tbl[,k]=round(tbl[,k],2)
}
for (k in c("pv")) {
    tbl[,k]=signif(tbl[,k],2)
}

colId=which(!names(tbl)%in%c("se","stat"))
pThres=0.05
tbl[out$pv<pThres,colId]

##############################################
## Among the 7 probes selected as DMRs by Kerkel et al. (cg07991621, cg08822227, cg09554443, cg05590257, cg14972143, cg00983520, cg21053323) all but one (cg05590257)

candGene=data.frame(cpgId=c("cg07991621","cg08822227","cg09554443","cg05590257","cg14972143","cg00983520","cg21053323","cg05590257"),stringsAsFactors=F)

##############################################
