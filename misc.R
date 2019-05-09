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

cat("\n================== misc.R ===============\n\n")

## ----------------------------------
## Parameters

cat("------------------- Parameters ---------------\n\n")

cohortFlag="_periDsal"

dirOut="periDsal/"

dateFlag="_20181025"

## ----------------------------------
#Reading in the data

##Set the data directory to the folder containing the Idat files as well as the Sample Excel sheet
if (computerFlag=="cluster") {
    dirCom="/home/royr/project/JoeWiemels/data/periDsal/"
    dirMeth="/cbc2/data1/ritu/wiemelsJ_periDsal/idat/"
    dirClin="/home/royr/project/JoeWiemels/data/periDsal/"
    fNameClin="IDATsControls and other QC PERI-DSAL_UCB EPIC plates 3-12_6-20-2018HMH_.txt"
} else {
    dirCom="/Users/royr/UCSF/JoeWiemels/leukMeth/docs/periDsal/"
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

clin2=read.table(paste(dirCom,"perinatalBlindIds_SeasonOfBirth.csv",sep=""), sep=",", h=T, quote="", comment.char="",as.is=T,fill=T)
names(clin2)[match(c("Blindid","Y_birth_month","season"),names(clin2))]=c("subjectId","birthMonth","season")
j=which(!clin$subjectId%in%clin2$subjectId)
tmp=clin2[1:length(j),]
for (k in 1:ncol(tmp)) tmp[,k]=NA
tmp$subjectId=clin$subjectId[j]
clin2=rbind(clin2,tmp)
j=match(clin$subjectId,clin2$subjectId); j1=which(!is.na(j)); j2=j[j1]
clin=cbind(clin,clin2[j2,which(!names(clin2)%in%names(clin))])

tbl=read.table(paste("periDsal/stat_sample_betaNoobBmiq_periDsal_QC.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
clin=cbind(clin,quality=tbl$quality[match(clin$id,tbl$id)])
tbl=read.table(paste(dirCom,"cellCount_minfi_cordBlood_periDsal.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
clin=cbind(clin,tbl[match(clin$id,tbl$id),which(!names(tbl)%in%names(clin))])

write.table(clin,file=paste("sampleInfo",cohortFlag,dateFlag,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F,append=F)

####################################################################
####################################################################

dateFlag="_20181119"
cohortFlag="_periDsal"

dirClin="/Users/royr/UCSF/JoeWiemels/leukMeth/docs/periDsal/"
fNameClin="sampleInfo_periDsal_20181025.txt"
clin=read.table(paste(dirClin,fNameClin,sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)

clin2=read.table(paste("/Users/royr/UCSF/JoeWiemels/leukMeth/docs/pandeyP/all_unique_covariates_for_prenatal.csv",sep=""), sep=",", h=T, quote="", comment.char="",as.is=T,fill=T)
names(clin2)[match(c("blindid","Age.at.blood.collection..hours.","peri_dsal","sex","caco","season","Y_birth_year","Y_birth_month","Y_birthweight","Y_gestation_week","Y_gestation_day","Y_delivery","Y_csection_type","CensusTract","MonthPrenatalCareBegan","NumberOfPrenatalCareVisits","ExpectedPrincipalPaymentPrenatal","ExpectedPrincipalPaymentDelivery","AGE"),names(clin2))]=
c("id","ageBloodCollect","periDsal","sex","caco","season","birthYear","birthMonth","birthWt","gestWeek","gestDay","delivery","csecType","censusTract","monthPrenatStart","numPrenatVisit","payPrenatal","payDelivery","age")

names(clin2)[!names(clin2)%in%names(clin)]

out=as.data.frame(t(sapply(clin$sampleName,function(x) {
    y=strsplit(x,"-")[[1]]
    if (length(y)==1) y=c(y,"")
    if (y[1]=="Jurkat") y=c(x,"")
    y
},USE.NAMES=F)))
names(out)=c("subjectId","replicate")
clin$replicate=out$replicate
j=which(!clin$subjectId%in%clin$subjectId[clin$id%in%clin2$id])
clin$keep=0
clin$keep[which(!duplicated(clin$subjectId))]=1
clin$keep[which(clin$id%in%clin2$id)]=0

j=match(clin$id,clin2$id); j1=which(!is.na(j)); j2=j[j1]
k=match(names(clin),names(clin2)); k1=which(!is.na(k)); k2=k[k1]
for (k in 1:length(k1)) {
    if (any(clin[j1,k1[k]]!=clin2[j2,k2[k]],na.rm=T) | any(is.na(clin[j1,k1[k]])!=is.na(clin2[j2,k2[k]]))) print(k)
}

j2=which(!clin2$id%in%clin$id)
tmp=clin[1:length(j2),]
for (k in 1:ncol(tmp)) tmp[,k]=NA
tmp$id=clin2$id[j2]
k=match(names(tmp),names(clin2)); k1=which(!is.na(k)); k2=k[k1]
for (k in 1:length(k1)) tmp[,k1[k]]=clin2[j2,k2[k]]
phen=rbind(clin,tmp)
phen$keepFromPriya=0
phen$keepFromPriya[which(phen$id%in%clin2$id)]=1
phen$keep[which(phen$id%in%clin2$id)]=1

j1=which(!phen$id%in%clin2$id)
tmp=clin2[1:length(j1),]
for (k in 1:ncol(tmp)) tmp[,k]=NA
tmp$id=phen$id[j1]
k=match(names(tmp),names(phen)); k1=which(!is.na(k)); k2=k[k1]
for (k in 1:length(k1)) tmp[,k1[k]]=phen[j1,k2[k]]
phen2=rbind(clin2,tmp)

j=match(phen$id,phen2$id); j1=which(!is.na(j)); j2=j[j1]
phen=cbind(phen[j1,],phen2[j2,which(!names(phen2)%in%names(phen))])

phen$id=gsub("-","_",phen$id)

write.table(phen,file=paste("sampleInfo",cohortFlag,dateFlag,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F,append=F)

####################################################################
####################################################################
