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
