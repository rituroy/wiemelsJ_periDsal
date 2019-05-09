####### Cell-mixture estimation : NEW METHOD ##########

## ---------------------------------

computerFlag="cluster"
computerFlag=""

nProbe=-1
nProbe=10001

## ---------------------------------

if (computerFlag=="cluster") {
    #setwd("/home/royr/project/JoeWiemels")
} else {
    dirSrc="/Users/royr/UCSF/"
    dirSrc2=dirSrc
    setwd(paste(dirSrc2,"JoeWiemels/leukMeth/epic",sep=""))
}

## ---------------------------------
# Packages
# ----------
# library(limma)
# library(qqman)
# library(scales)
# library(bigpca)
library(FactoMineR)

# Functions
#------------------
# return the variable position in dataset

number <- function(data) {
    data_matrix <- matrix(c(1:dim(data)[2],(names(data))),,2)
    data_matrix
}

# Remove NAs if more than 3% of the row
remove_NAs <- function(x) {
    if (table(is.na(x))/ length(x) < 0.97) {TRUE} else {FALSE}
}

# Impute mean values on missing data

impute_mean <- function(x) {
    ifelse (is.na(x), mean(x, na.rm = TRUE), x)
}

# Load the clinical data
# --------------------

#load("/Users/sgonseth/Desktop/methylation_ucsf/reFACTor/all_phenotypical_variables.RData")

# Load the annotations file

cat("\n\n============================ cell-mixt_reFACTor_BMIQ_V4_acc_to_matlab_code_RR.R ===========================\n\n",sep="")

## ----------------------------------------------
if (computerFlag=="cluster") {
    datadir="results/"
} else {
    datadir="results/"
}
load(paste(datadir,"ann850k.RData",sep=""))
ann=as.data.frame(ann850k)
rm(ann850k)
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

## ----------------------------------------------

datType="_aml"
datType="_allGuthSet2"
datType="_allGuthSet1"
datType="_allGuthSet1Set2"
datType="_allGuthSet1Set2_ctrlSubset"
datType="_allGuthSet1Set2_ctrlSubsetFunNorm_ctrlSubset"
datType="_allGuthSet1Set2_caseSubsetFunNorm_caseSubset"
datType="_periDsal"
datType="_periDsal_periSubset"
datType="_periDsal_dsalSubset"

cat("\n\n============================",datType,"===========================\n\n")

if (T) {
    #for (datType in c("_allGuthSet1","_allGuthSet2","_aml")) {
    
    if (computerFlag=="cluster") {
        dirClin=dirMeth="data/"
        switch(datType,
            "_allGuthSet2"={
                dirClin=dirMeth="data/set2/"
                fNameMeth=paste("beta_bmiq",datType,".txt",sep="")
                fNameClin=paste("clin_guthrieSet2_20140619.txt",sep="")
            },
            "_allGuthSet1"={
                dirClin=dirMeth="data/set1/"
                fNameMeth=paste("beta_bmiq",datType,".txt",sep="")
                fNameClin=paste("final.txt",sep="")
            },
            "_allGuthSet1Set2"={
                dirClin=dirMeth="data/set1set2/"
                fNameMeth=paste("beta_bmiq",datType,".txt",sep="")
                fNameClin=paste("clin_allGuthSet1Set2_20160523.txt",sep="")
            },
            "_allGuthSet1Set2_ctrlSubset"={
                dirClin=dirMeth="data/set1set2/"
                fNameMeth=paste("beta_bmiq",datType,".txt",sep="")
                fNameClin=paste("clin_allGuthSet1Set2_20160523.txt",sep="")
            },
            "_allGuthSet1Set2_ctrlSubsetFunNorm_ctrlSubset"={
                dirClin=dirMeth="data/set1set2/"
                fNameMeth=paste("beta_bmiq",datType,".txt",sep="")
                fNameClin=paste("clin_allGuthSet1Set2_20160523.txt",sep="")
            },
            "_allGuthSet1Set2_caseSubsetFunNorm_caseSubset"={
                dirClin=dirMeth="data/set1set2/"
                fNameMeth=paste("beta_bmiq",datType,".txt",sep="")
                fNameClin=paste("clin_allGuthSet1Set2_20160523.txt",sep="")
            },
            "_aml"={
                dirClin=dirMeth="data/aml/"
                fNameMeth=paste("beta_bmiq_aml.txt",sep="")
                fNameClin=paste("clin_aml_20150114.txt",sep="")
            },
            "_periDsal"={
                dirMeth="results/"
                dirClin="../data/periDsal/"
                fNameMeth="betaBmiq.RData"
                fNameClin="sampleInfo_periDsal_20181119.txt"
            },
            "_periDsal_periSubset"={
                dirMeth="results/"
                dirClin="../data/periDsal/"
                fNameMeth="betaBmiq.RData"
                fNameClin="sampleInfo_periDsal_20181119.txt"
            },
            "_periDsal_dsalSubset"={
                dirMeth="results/"
                dirClin="../data/periDsal/"
                fNameMeth="betaBmiq.RData"
                fNameClin="sampleInfo_periDsal_20181119.txt"
            }
        )
    } else {
            switch(datType,
            "_allGuthSet2"={
                dirClin=dirMeth="docs/all/set2/"
                fNameMeth=paste("beta_bmiq",datType,".txt",sep="")
                fNameClin=paste("clin_guthrieSet2_20140619.txt",sep="")
            },
            "_allGuthSet1"={
                dirClin=dirMeth="docs/all/set1/"
                fNameMeth=paste("beta_bmiq",datType,".txt",sep="")
                fNameClin=paste("final.txt",sep="")
            },
            "_allGuthSet1Set2"={
                dirClin=dirMeth="docs/all/set1set2/"
                fNameMeth=paste("beta_bmiq",datType,".txt",sep="")
                fNameClin=paste("clin_allGuthSet1Set2_20160523.txt",sep="")
            },
            "_aml"={
                dirClin=dirMeth="docs/aml/"
                fNameMeth=paste("beta_bmiq_aml.txt",sep="")
                fNameClin=paste("clin_aml_20150114.txt",sep="")
            },
            "_periDsal"={
                dirMeth="/Users/royr/UCSF/JoeWiemels/leukMeth/epic/results/"
                dirClin="/Users/royr/UCSF/JoeWiemels/leukMeth/docs/periDsal/"
                fNameMeth="betaBmiq.RData"
                fNameClin="sampleInfo_periDsal_20181119.txt"
            },
            "_periDsal_periSubset"={
                dirMeth="/Users/royr/UCSF/JoeWiemels/leukMeth/epic/results/"
                dirClin="/Users/royr/UCSF/JoeWiemels/leukMeth/docs/periDsal/"
                fNameMeth="betaBmiq.RData"
                fNameClin="sampleInfo_periDsal_20181119.txt"
            },
            "_periDsal_dsalSubset"={
                dirMeth="/Users/royr/UCSF/JoeWiemels/leukMeth/epic/results/"
                dirClin="/Users/royr/UCSF/JoeWiemels/leukMeth/docs/periDsal/"
                fNameMeth="betaBmiq.RData"
                fNameClin="sampleInfo_periDsal_20181119.txt"
            }
        )
    }
    phen=read.table(paste(dirClin,fNameClin,sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    if (length(grep("txt",fNameMeth))==1) {
        meth=read.table(paste(dirMeth,fNameMeth,sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=nProbe)
        rownames(meth)=meth$probeId
        meth=as.matrix(meth[,-1])
    } else {
        load(paste(dirMeth,fNameMeth,sep=""))
        meth=betaBmiq
        rm(betaBmiq)
        if (nProbe!=-1) meth=meth[1:nProbe,]
    }
    switch(datType,
    "_allGuthSet2"={
        names(phen)[match(c("subjectID","birth_wt"),names(phen))]=c("subjectId","birthWt")
        phen$id=paste("X",phen$guthrieId,sep="")
    },
    "_allGuthSet1"={
        names(phen)[match(c("Subject_ID","birth_weight"),names(phen))]=c("subjectId","birthWt")
        phen$id=paste("X",phen$TargetID,sep="")
    },
    "_aml"={
        phen$id=paste("X",phen$guthrieId,sep="")
    },
    "_periDsal_periSubset"={
        phen=phen[which(phen$periDsal=="peri"),]
    },
    "_periDsal_dsalSubset"={
        phen=phen[which(phen$periDsal=="dsal"),]
    }
    )
    j=match(colnames(meth),phen$id); j1=which(!is.na(j)); j2=j[j1]
    meth=meth[,j1]
    phen=phen[j2,]
    rownames(phen)=phen$id
    if (F) {
        switch(datType,
        "_allGuthSet2"={
            clin_set1=phen
            beta_set2=meth
        },
        "_allGuthSet1"={
            clin_set2=phen
            beta_set1=meth
        },
        "_aml"={
            clin_aml=phen
            beta_aml=meth
        }
        )
    }
    #}
    clinThis=phen
    betaThis=meth
    #save(betaThis,clinThis,file="tmp.RData")
    rm(phen,meth)
    save.image(file=paste("tmp",datType,".RData",sep=""))
    #rm(clin_aml,beta_aml)
    
    #require(stringr)
    #colnames(set1) <- str_sub(colnames(set1), 2,end=100)
    
    
}


# Load outliers found by Elior
# ------------------------


if (T) {
    ## Find sample outliers
    
    load(file=paste("tmp",datType,".RData",sep=""))
    i=match(rownames(betaThis),ann$IlmnID); i1=which(!is.na(i)); i2=i[i1]
    betaThis=betaThis[i1,]
    ann=ann[i2,]
    
    IDs_to_keep=1:ncol(betaThis)
    x=apply(betaThis,1,function(x) mean(!is.na(x)))
    CpG_to_keep <- which(ann$snp==0 & ann$CHR%in%1:22 & x>=0.03)
    
    betaThis <- betaThis[CpG_to_keep,IDs_to_keep]
    betaThis=t(betaThis)
    R_est <- PCA(betaThis, graph = F, ncp =6)
    save(R_est,file=paste("R_est_init",datType,".RData",sep=""))
    Refactor_dat <- R_est$ind$coord[,c(1:6)]
    thPCA=c(-200,200)
    if (datType%in%c("_allGuthSet1Set2_ctrlSubsetFunNorm_ctrlSubset","_allGuthSet1Set2_caseSubsetFunNorm_caseSubset")) thPCA=c(-500,500)
    if (datType%in%c("_periDsal")) thPCA=c(-1000,1000)
    png(paste("pca_1",datType,".png",sep=""))
    plot(R_est,xlim=c(-2000,1000),ylim=c(-2000,1000))
    abline(c(0,1),lty="dotted")
    dev.off()
    ## Outliers are those outside the dotted lines
    png(paste("pca_1_1",datType,".png",sep=""))
    plot(R_est$ind$coord[,1],R_est$ind$coord[,2],xlim=c(-2000,1000),ylim=c(-2000,1000))
    i=which(abs(R_est$ind$coord[,1])>thPCA[2] | abs(R_est$ind$coord[,2])>thPCA[2])
    points(R_est$ind$coord[i,1],R_est$ind$coord[i,2],,col="green")
    abline(c(0,1),lty="dotted"); abline(h=0,lty="dotted"); abline(v=0,lty="dotted")
    abline(h=thPCA,lty="dotted")
    abline(v=thPCA,lty="dotted")
    dev.off()
    png(paste("pca_1_2",datType,".png",sep=""))
    plot(R_est,xlim=c(-2000,1000),ylim=c(-2000,1000),label="none")
    abline(c(0,1),lty="dotted")
    dev.off()
    
    load(file=paste("R_est_init",datType,".RData",sep=""))
    table((abs(R_est$ind$coord[,1])>thPCA[2] | abs(R_est$ind$coord[,2])>thPCA[2]))
    ## Outliers
    paste(rownames(R_est$ind$coord)[which(abs(R_est$ind$coord[,1])>thPCA[2] | abs(R_est$ind$coord[,2])>thPCA[2])],collapse=",")
}

if (F) {

    load(file=paste("tmp",datType,".RData",sep=""))
    i=match(rownames(betaThis),ann$IlmnID); i1=which(!is.na(i)); i2=i[i1]
    betaThis=betaThis[i1,]
    ann=ann[i2,]

    ## Exclude outliers
    IDs_to_keep=1:ncol(betaThis)
    IDs_to_keep=which(!colnames(betaThis)%in%c("X1288G","X1466G","X0153G","X0077G","X0201G","X1191G","X1264G","X1219G","X0873G","X0927G","X1075G","X1766G","X0691G"))
    switch(datType,
        "_allGuthSet1"={
            IDs_to_keep=which(!colnames(betaThis)%in%c("X1218G","X0580G","X0404G","X1336G","X1225G","X1202G","X1162G","X1224G","X0527G","X0665G","X0451G","X0557G","X1664G","X0498G","X0419G","X1446G","X2006G","X1577G","X0486G","X1419G","X1739G","X1540G","X1736G","X1988G","X1718G","X0525G","X2221G","X1749G","X0420G","X1769G","X0511G","X0563G","X1650G","X1237G","X1328G","X1297G","X0993G","X1114G","X1042G","X1116G","X1870G","X1240G","X1967G","X1011G","X1439G","X1906G","X1070G","X1933G","X0935G","X1244G","X1639G","X1384G","X2045G","X1246G","X1579G","X1205G","X0489G","X0280G","X0183G","X0924G","X0470G","X0162G","X0481G","X1288G","X1578G","X1037G","X1157G","X0384G","X0283G","X0716G","X1323G","X1965G","X0804G","X1392G","X1221G","X0284G","X0194G","X0264G","X0285G","X0459G","X1141G","X1098G","X0594G","X0586G","X1280G","X0401G","X0136G","X0504G","X0332G","X0381G","X0584G","X0541G"))
        },
        "_allGuthSet1Set2_ctrlSubset"={
            IDs_to_keep=which(!colnames(betaThis)%in%c("X1218G","X0580G","X1336G","X1204G","X1225G","X1202G","X1252G","X1162G","X1224G","X0527G","X0665G","X0465G","X0471G","X1477G","X1562G","X0486G","X1718G","X0525G","X1769G","X1650G","X0311G","X0469G","X1248G","X1639G","X1384G","X1579G","X0885G","X1341G","X1008G","X1205G","X1117G","X1405G","X1288G","X1578G","X1260G","X0874G","X1700G","X1326G","X1164G","X1359G","X1392G","X1147G","X1221G","X1141G","X1098G","X0594G","X0586G","X1280G","X1404G","X0552G","X2046G","X1770G","X1824G","X2029G","X0247G","X1658G","X0080G","X1360G","X1024G","X1173G","X0229G","X1021G","X1994G","X0132G","X0091G","X0137G","X0197G","X0531G","X1881G","X1781G","X0784G","X0184G","X1285G","X0053G","X1230G","X1920G","X0380G","X0231G","X0876G","X0125G","X0815G","X0036G","X0643G","X1628G","X1092G","X1901G","X0827G","X1169G","X1807G","X1884G","X1777G","X0513G","X0270G","X0437G","X1834G","X1946G","X1277G","X1843G","X1859G","X0173G","X0107G","X1123G","X1996G","X2188G","X0065G","X0225G","X0149G","X1184G"))
        },
        "_allGuthSet1Set2_ctrlSubsetFunNorm_ctrlSubset"={
            IDs_to_keep=which(!colnames(betaThis)%in%c("X1218G","X1336G","X1225G","X1202G","X1224G","X1650G","X1205G","X1288G","X1392G","X1221G","X0594G","X2046G","X0137G","X1881G","X0053G","X1230G","X0125G","X1277G","X1184G"))
        },
        "_allGuthSet1Set2_caseSubsetFunNorm_caseSubset"={
            IDs_to_keep=which(!colnames(betaThis)%in%c("X1237G","X1967G","X1906G","X1933G","X1246G","X0162G","X1157G","X1965G","X0381G","X1111G","X1174G","X1168G","X0111G","X1458G","X1256G","X1139G","X1227G"))
        },
        "_periDsal"={
            IDs_to_keep=which(!colnames(betaThis)%in%c("X198036","X198443","X198384","X198394","X198029","X198226","X198509","X198454","X198271","X198124","X198437","X198347","X198456","X198339","X198087","X198088","X198120","X724984","X198279","X198349","X198145","X198216","X198341","X198056","X198203","X198467","X198163","X198002","X198310","X198174","X198172","X198428","X198189","X198042","X198135","X198001"))
        }
    )
    x=apply(betaThis,1,function(x) mean(!is.na(x)))
    CpG_to_keep <- which(ann$snp==0 & ann$CHR%in%1:22 & x>=0.03)
    betaThis <- betaThis[CpG_to_keep,IDs_to_keep]
    ann=ann[CpG_to_keep,]
    clinThis=clinThis[IDs_to_keep,]


    # # Removing CpG sites with >3% of NAs
    # # --------------------------------
    # # set 1
    # remove_NAs <- function(x){
    # if (table(is.na(x))/ length(x) < 0.97) {TRUE} else {FALSE}
    # }

    # res_set1 <- apply(set1, 1, remove_NAs)

    # set1_nona <- set1[(res_set1) == F,]; dim(set1_nona)

    # # set 2

    # res_set2 <- apply(set2, 1, remove_NAs)

    # set2_nona <- set2[(res_set2) == F,]; dim(set2_nona)

    # Imputing mean methylation levels on the other NAs
    # --------------------------------------------

    # set 2
    betaThis <- t(apply(betaThis, 1, impute_mean))


    # reFACTor calculations
    # -------------------------
    # Reduction of the number of CpG sites (acc. to Elior Rahmani)

    # x1 <- set1
    x2 <- betaThis

    # mean_x1 <- apply(x1, 1, mean)
    # names(mean_x1) <- rownames(x1)
    mean_x2 <- apply(x2, 1, mean)
    names(mean_x2) <-rownames(x2)

    x2_proc <- mean_x2[mean_x2 < 0.9 & mean_x2 > 0.1]

    betaThis <- x2[which(rownames(x2) %in% names(x2_proc)),]

    #save(betaThis, file = paste("beta_matrix",datType,".RData",sep=""))
    #load(file = paste("beta_matrix",datType,".RData",sep=""))

    # refactor from Elior''s matlab code

    # Running a standard PCA...
    # ---------------------------------

    O_prime <- (betaThis)

    # z score
    zscore <- function(x){
        (x - mean(x)) / sd(x)
    }

    z_norm_O_prime <- t(apply(O_prime, 1, zscore))

    K = 6

    t0 <- Sys.time()
    res_PCA <- princomp(z_norm_O_prime)
    Sys.time() - t0


    save(res_PCA, file = paste("res_PCA_for_reFACTOR",datType,".RData",sep=""))

    # Compute a low rank approximation of input data and rank sites...
    # ---------------------------------

    x <- res_PCA$scores[,c(1:K)] %*% t((res_PCA$loadings[,c(1:K)]))

    # An = bsxfun(@minus,O',mean(O',1));
    An <- z_norm_O_prime

    # Bn=bsxfun(@minus,x,mean(x,1));
    Bn <- t(apply(x, 1, zscore))

    # An=bsxfun(@times,An,1./sqrt(sum(An.^2,1)));

    low_rank <- function(data){
        data2 <- (data^2)
        data3 <- sum(data2)
        data4 <- sqrt(data3)
        data5 <- 1/data4
        data*data5
    }
    An2 <- t(apply(An, 1, low_rank))

    # Bn=bsxfun(@times,Bn,1./sqrt(sum(Bn.^2,1)));

    Bn2 <- t(apply(Bn, 1, low_rank))


    # Find the distance of each site from its low rank approximation
    # ---------------------------------

    res_AnBn <- An2-Bn2
    # save(res_AnBn, file = "For_cluser_res_AnBn.RData")

    #load("For_cluser_res_AnBn.RData")

    # distances = sum((An-Bn).^2,1).^0.5 ;
    # [~,ranked_list] = sort(distances);

    res_interm <- apply((res_AnBn)^2, 1, FUN = sum)
    distances <- sqrt(res_interm)
    ranked_list <- sort(distances)

    # Compute ReFACTor components...
    # ---------------------------------
    t = 500
    # sites = ranked_list(1:t);
    sites <- ranked_list[1:t]

    # [~,R_est] = pca(zscore(O(sites,:)''));

    t_O_sites <- t(O_prime[which(rownames(z_norm_O_prime) %in% names(sites)),])

    t_O_sites_zs <- (apply(t_O_sites, 2, zscore))

    R_est <- PCA(t_O_sites_zs, graph = F, ncp =6)

    Refactor_dat <- R_est$ind$coord[,c(1:6)]

    save(Refactor_dat, file = paste("Refactor_dat",datType,".RData",sep=""))

    #load(file = paste("Refactor_dat",datType,".RData",sep=""))

}
