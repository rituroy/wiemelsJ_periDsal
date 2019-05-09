
if (T) {
argv=commandArgs(TRUE)
modelId=argv[1]
subsetFlag=argv[2]
}

if (T) {
modelId="season_covSexPlate"
modelId="season_covPlate"
modelId="season_covSexPlateCelltype"
subsetFlag="_periCtrl"

modelId="dsalPeri_covSexPlateCelltype"
subsetFlag="_case"

modelId="dsalPeri_covSexPlateCelltype"
subsetFlag="_ctrl"

modelId="caco_covSexPlateCelltype"
subsetFlag="_peri"

modelId="caco_covSexPlateCelltype"
subsetFlag="_dsal"
}

##############################################

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

# libraries
library(R.utils)
library(data.table)# to process results
library(MASS) # rlm function for robust linear regression
library(sandwich) # Huberís estimation of the standard error
library(lmtest) # to use coeftest

## --------------
workdir="results/"
#typeFlag=""
#projectFlag=""
#studyFlag=""

## --------------

#subset1Flag="ctrl"

## --------------

#predTypeFlag="_season"
#covFlag=""


x=strsplit(modelId,"_")[[1]]
predVarFlag=predVar2Flag=x[1]
cov2Flag=paste("_",x[2],sep="")
modelThis=predVar2Flag
switch(predVar2Flag,
    "season"={
        modelThis="seasonOfBirth"
    },
    "caco"={
    },
    "dsalPeri"={
        predVarFlag="periDsal"
    }
)
switch(cov2Flag,
    "_covPlate"={
        modelThis=paste(modelThis," + Batch",sep="")
    },
    "_covariate"={
        modelThis=paste(modelThis," + Sex + Mat_age + Gest_age + Msmk + SocioEc + Batch",sep="")
    },
    "_covSexPlateCelltype"={
        modelThis=paste(modelThis," + Sex + Batch + nRBC + CD8T + CD4T + NK + Bcell + Mono + Gran",sep="")
    },
    "_covSexPlate"={
        modelThis=paste(modelThis," + Sex + Batch",sep="")
    }
)
if (F) {
    switch(modelId,
        "season_covPlate"={
            predVarFlag="season"; cellFlag="NoHouse"
            cov2Flag <- "_covPlate"
        },
        "season_covariate"={
            predVarFlag="season"; cellFlag="NoHouse"
            cov2Flag <- "_covariate"
        },
        "season_covSexPlate"={
            predVarFlag="season"; cellFlag="House"
            cov2Flag <- "_covSexPlate"
        },
        "season_covSexPlateCelltype"={
            predVarFlag="season"; cellFlag="House"
            cov2Flag <- "_covSexPlateCelltype"
        },
        "caco_covSexPlateCelltype"={
            predVarFlag="caco"; cellFlag="House"
            cov2Flag <- "_covSexPlateCelltype"
        },
        "caco_covSexPlateCelltype"={
            predVarFlag="caco"; cellFlag="House"
            cov2Flag <- "_covSexPlateCelltype"
        },
        "dsalPeri_covSexPlateCelltype"={
            predVarFlag="periDsal"; cellFlag="House"
            cov2Flag <- "_covSexPlateCelltype"
        },
        "dsalPeri_covSexPlateCelltype"={
            predVarFlag="periDsal"; cellFlag="House"
            cov2Flag <- "_covSexPlateCelltype"
        }
    )

    compCellTypeName="_cordBlood"; cellMixCutoff=0.1
}

## --------------

## ---------------------------------

nProbe=-1
nProbe=101

## ---------------------------------

cohortFlag="_allGuthSet1Set2"
cohortFlag="_periDsal"

normFlag=""
normFlag="_funNorm"
normFlag="_bmiq"

## ---------------------------------

##############################################

heading=paste(c(", ",subsetFlag,", ",", ",cohortFlag),collapse="")
cat("\n\n============================ EPIC: lingReg_pace ===========================\n\n")
cat("\n\n============================",subsetFlag,", ",cohortFlag,", ",modelId,"===========================\n\n")

timeStamp1=Sys.time()
cat(format(timeStamp1, "%x %X"),"\n",sep="")

dirMethRef=""
if (computerFlag=="cluster") {
    dirClin=dirMeth="data/"
    switch(cohortFlag,
        "_periDsal"={
            dirMeth="/cbc2/data1/ritu/wiemelsJ_periDsal/idat/"
            dirClin="/home/royr/project/JoeWiemels/data/periDsal/"
            fNameMeth="betaBmiq"
            #fNameClin="sampleInfo_periDsal_20181025.txt"
            fNameClin="sampleInfo_periDsal_20181119.txt"
        }
    )
} else {
    switch(cohortFlag,
        "_periDsal"={
            dirMeth="/Users/royr/UCSF/JoeWiemels/leukMeth/docs/periDsal/idat/"
            dirClin="/Users/royr/UCSF/JoeWiemels/leukMeth/docs/periDsal/"
            fNameMeth="betaBmiq"
            #fNameClin="sampleInfo_periDsal_20181025.txt"
            fNameClin="sampleInfo_periDsal_20181119.txt"
        }
    )
}

switch(cohortFlag,
    "_periDsal"={
        phen=read.table(paste(dirClin,fNameClin,sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
        #tbl=read.table(paste(dirClin,"cellCount_minfi_cordBlood_periDsal.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
        #phen=cbind(phen,tbl[match(phen$id,tbl$id),which(!names(tbl)%in%names(phen))])
        phen=phen[which(phen$keep==1),]
        
        if (nProbe==-1) {
            load(file=paste(workdir,fNameMeth,".RData",sep=""))
            if (nProbe!=-1) betaBmiq=betaBmiq[1:nProbe,]
            
            # transpose betas so that rows are samples and columns are probes
            beta_matrix <- t(betaBmiq)
            rm(betaBmiq)
            save(beta_matrix,file="beta_matrix.RData")
        } else {
            load(file="beta_matrix.RData")
        }
        rownames(beta_matrix)=gsub("-","_",rownames(beta_matrix))
        j=match(rownames(beta_matrix),phen$id); j1=which(!is.na(j)); j2=j[j1]
        beta_matrix=beta_matrix[j1,]
        phen=phen[j2,]
        
        samId=1:nrow(phen)
        samId=which(phen$quality!="3. low from bmiq" & !is.na(phen$periDsal) & !is.na(phen$caco))
        switch(subsetFlag,
            "_case"={
                samId=samId[which(phen$caco[samId]=="case")]
                subsetName="_subsetCase"
            },
            "_ctrl"={
                samId=samId[which(phen$caco[samId]=="control")]
                subsetName="_subsetCtrl"
            },
            "_peri"={
                samId=samId[which(phen$periDsal[samId]=="peri")]
                subsetName="_subsetPeri"
            },
            "_dsal"={
                samId=samId[which(phen$periDsal[samId]=="dsal")]
                subsetName="_subsetDsal"
            },
            "_periCtrl"={
                samId=samId[which(phen$periDsal[samId]=="peri" & phen$caco[samId]=="control")]
                subsetName="_subsetPeriCtrl"
            }
        )
        beta_matrix=beta_matrix[samId,]
        phen=phen[samId,]
    }
)


##############################################


# set variables for design formula
# NOTE: season names for birth and collection must be named differently or you will get an error in the model.
# If you include additional covariates, add them below and in the design object further below.
phen$predVar=phen[,predVarFlag]
levInfo=data.frame(lev1=sort(unique(phen$predVar)),lev2=sort(unique(phen$predVar)),stringsAsFactors=F)
switch(predVar2Flag,
    "season"={
        phen$predVar[which(phen$season=="autumn")]="1_autumn"
        phen$predVar[which(phen$season=="winter")]="2_winter"
        phen$predVar[which(phen$season=="spring")]="3_spring"
        phen$predVar[which(phen$season=="summer")]="4_summer"
        seasonOfBirth=factor(phen$predVar,levels=c("1_autumn","2_winter","3_spring","4_summer"))
        levInfo=data.frame(lev1=sort(unique(phen$season)),lev2=sort(unique(seasonOfBirth)),stringsAsFactors=F)
    },
    "caco"={
        caco=as.integer(phen$predVar=="case")
    },
    "dsalPeri"={
        dsalPeri=as.integer(phen$predVar=="dsal")
    }
)

#  Change as needed, but keep levels for season the same.
Sex <- factor(phen$sex)

if (F) {
    Mat_age <- phen$mo_age
    Gest_age <- phen$gestage
    Msmk <- factor(phen$smoke_mo_preg)
    SocioEc <- phen$income3
    Batch <- phen$Beadchip
    SocioEc <- factor(phen$income3)
    Batch <- factor(paste(phen$set,phen$Batch))
}

#Batch <- factor(phen$beadchip)
Batch <- factor(phen$plate)

if (T) {
    nRBC=phen$nRBC
    CD8T=phen$CD8T
    CD4T=phen$CD4T
    NK=phen$NK
    Bcell=phen$Bcell
    Mono=phen$Mono
    Gran=phen$Gran
}

if (F) {
    switch(modelId,
        "season_covPlate"={
            modelThis="seasonOfBirth + Batch"
        },
        "season_covariate"={
            modelThis="seasonOfBirth + Sex + Mat_age + Gest_age + Msmk + SocioEc + Batch"
        },
        "season_covSexPlateCelltype"={
            modelThis="seasonOfBirth + Sex + Batch + nRBC + CD8T + CD4T + NK + Bcell + Mono + Gran"
        },
        "season_covSexPlate"={
            modelThis="seasonOfBirth + Sex + Batch"
        }
)
}

design <- as.formula(paste("beta_matrix[,methcol] ~ ",modelThis,sep=""))
cat("\n\n============================",modelThis,"===========================\n\n")

# Add function for running the model.  This function uses the object “design”.  If you change the model, adjust that object above and make sure the additional covariates exist in your workspace.  This function will take the first 3 covariates which should be birth during winter, spring, and summer with autumn as the reference.  Occasionally, you may run into an error with the coeftest function.  This function is designed to identify this error and produce NAs as output when this occurs instead of halting the function.

RLMtest.my=function(methcol,design,numCol) {
    mod = try(rlm(design,maxit=200))
    if ("try-error"%in%class(mod)){
        print(paste("error thrown by column", methcol))
        return(invisible(rep(NA, numCol*3)))
    } else {cf = try(coeftest(mod, vcov=vcovHC(mod, type="HC0")))}
    if ("try-error"%in%class(cf)){
        print(paste("error in coeftest by column", methcol,"setting values to NA"))
        return(invisible(rep(NA,numCol*3)))
    } else {
        keep <- c("Estimate","Std. Error","Pr(>|z|)")
        #  Pull out winter, spring, and summer.
        output=c()
        for (k in 1:length(numCol)) output=c(output,cf[k+1,keep])
        names1 <- rownames(cf)[(1:numCol)+1]
        names2 <- c("BETA","SE","P_VAL")
        names_full=matrix(NA,nrow=length(names2),ncol=length(names1))
        for (i in 1:length(names1)){
            for (j in 1:length(names2)){
                names_full[j,i] <- paste(names1[i],names2[j],sep="_")
            }
        }
        names_full <- as.vector(names_full)
        names(output) <- names_full
        return(output)
    }
}

numCol=(nrow(levInfo)-1)
colId=c("coef","se","pv")
switch(predVar2Flag,
    "season"={
        colId=paste(colId,"_seasonOfBirth",rep(c("Winter","Spring","Summer"),each=length(colId)),sep="")
    },
    "caco"={
        colId=paste(colId,"_",predVar2Flag,sep="")
    },
    "dsalPeri"={
        colId=paste(colId,"_",predVar2Flag,sep="")
    }
)
all.results=matrix(nrow=ncol(beta_matrix),ncol=numCol*length(colId),dimnames=list(colnames(beta_matrix),colId))
for (methcol in 1:ncol(beta_matrix)) {
    all.results[methcol,]=RLMtest.my(methcol,design,numCol)
}
all.results=data.table(probeID=colnames(beta_matrix),all.results)

colnames(all.results)

datadir="tmp/"
fName=paste("_",predVar2Flag,cov2Flag,subsetName,cohortFlag,sep="")
#save(all.results,file=paste("stat",fName,".RData",sep=""))
write.table(all.results, paste(datadir,"stat",fName,".txt",sep=""),na="NA",row.names=FALSE,sep="\t",quote=FALSE)
gzip(paste(datadir,"stat",fName,".txt",sep=""))

#table(is.na(all.results[,"pv_seasonOfBirthSummer"]))

##############################################
