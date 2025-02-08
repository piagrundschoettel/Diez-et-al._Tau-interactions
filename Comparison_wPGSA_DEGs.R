library(dplyr)

setwd("/Users/koki/Documents/Projects/Lisa/20241110 New corrected data/Analysis_ver_20241110/wPGSA/Comparison_wPGSA_DEGs")

KDvsBase_wPGSA_nonTE <- read.csv("../list_tscore4_p10Emin6/TauKDvsBase_nonTE.csv")
OEvsBase_wPGSA_nonTE <- read.csv("../list_tscore4_p10Emin6/TauOEvsBase_nonTE.csv")
OEvsKD_wPGSA_nonTE <- read.csv("../list_tscore4_p10Emin6/TauOEvsKD_nonTE.csv")

KDvsBase_DEGs <- read.csv("../DEG_Lists/KD_vs_Base.csv")
OEvsBase_DEGs <- read.csv("../DEG_Lists/OE_vs_Base.csv")
OEvsKD_DEGs <- read.csv("../DEG_Lists/OE_vs_KD.csv")

wPGSA_ID_conv<- read.csv("wPGSA_to_GeneSyn.csv")


KDvsBase_TF_Conv <- inner_join(wPGSA_ID_conv[,1:2], KDvsBase_wPGSA_nonTE, by = c("wPGSA_IDs" =  "TF"))
OEvsBase_TF_Conv <- inner_join(wPGSA_ID_conv[,1:2], OEvsBase_wPGSA_nonTE, by = c("wPGSA_IDs" =  "TF"))
OEvsKD_TF_Conv <- inner_join(wPGSA_ID_conv[,1:2], OEvsKD_wPGSA_nonTE, by = c("wPGSA_IDs" =  "TF"))

KDvsBase_DEGs$GENEID <- toupper(KDvsBase_DEGs$GENEID)
OEvsBase_DEGs$GENEID <- toupper(OEvsBase_DEGs$GENEID)
OEvsKD_DEGs$GENEID <- toupper(OEvsKD_DEGs$GENEID)


KDvsBase_DEG_TF_Match <- inner_join(KDvsBase_TF_Conv, KDvsBase_DEGs, by = c("Gene_Symbol_Mus" = "GENEID"))
OEvsBase_DEG_TF_Match <- inner_join(OEvsBase_TF_Conv, OEvsBase_DEGs, by = c("Gene_Symbol_Mus" = "GENEID"))
OEvsKD_DEG_TF_Match <- inner_join(OEvsKD_TF_Conv, OEvsKD_DEGs, by = c("Gene_Symbol_Mus" = "GENEID"))

testcounts<-  read.csv("test.csv")
testcounts <- testcounts[0:19213,]
test<- inner_join(wPGSA_ID_conv,testcounts, by =c( "Gene_Symbol_Mus" ="GENEID" ))

testantij <- anti_join(wPGSA_ID_conv, test)

write.csv(KDvsBase_DEG_TF_Match, file = "./KDvsBase_DEG_TF_Match.csv", row.names = F)
write.csv(OEvsBase_DEG_TF_Match, file = "./OEvsBase_DEG_TF_Match.csv", row.names = F)
write.csv(OEvsKD_DEG_TF_Match, file = "./OEvsKD_DEG_TF_Match.csv", row.names = F)

