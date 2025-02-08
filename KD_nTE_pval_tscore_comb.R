library(dplyr)
setwd("/Users/koki/Documents/Projects/Lisa/20241110 New corrected data/Analysis_ver_20241110/wPGSA/KDvsBase_nonTE")

tscore <- read.table("KDvsBase_nonTE_TF_wPGSA_t_score.txt",header = T)
pval <- read.table("KDvsBase_nonTE_TF_wPGSA_p_value.txt",header = T)

tp_join <- inner_join(pval[c(1,3)],tscore[c(1,3)], by="TF")

tp_join[c("FDR_p_val")] <-p.adjust(tp_join$p_value_mean, method="fdr")  

tp_exp <- tp_join

write.csv(tp_exp, file = "./Lisa_TauKDvsBase_nonTE_exp.csv", row.names = F)




#volcano plot
library(ggrepel)
library(EnhancedVolcano)


pdf("Lisa_TauKDvsBase_nonTE_volc.pdf", width=10, height=6)
EnhancedVolcano(tp_exp,
                lab = tp_exp[,1],
                x = "t_score_mean",
                y = "FDR_p_val",
                xlab = bquote('wPGSA t score'),
                ylab = bquote("-"~Log[10]~ 'Adj. PVal.'),
                pCutoff = 1*10^(-5),
                FCcutoff = 4,
                pointSize = 2,
                labSize = 5,
                colAlpha = 0.8,
                col= c("black","black","black","red"),
                legendLabSize = 12,
                legendIconSize = 4.0,
                legendPosition = "none",
                drawConnectors = T,
                widthConnectors = 0.5,
                colConnectors = 'black',
                xlim = c(-12,12),
                ylim = c(0,50),
                subtitle = NULL,
                title = NULL,
                max.overlaps = 10
)
dev.off()  

  tiff("Lisa_TauKDvsWT_nonTE_volc.tiff", units="cm", width=25, height=15, res=300)
EnhancedVolcano(tp_exp,
                lab = tp_exp[,1],
                x = "t_score_mean",
                y = "FDR_p_val",
                xlab = bquote('wPGSA t score'),
                ylab = bquote("-"~Log[10]~ 'Adj. PVal.'),
                pCutoff = 1*10^(-5),
                FCcutoff = 4,
                pointSize = 2,
                labSize = 5,
                colAlpha = 0.8,
                col= c("black","black","black","red"),
                legendLabSize = 12,
                legendIconSize = 4.0,
                legendPosition = "none",
                drawConnectors = T,
                widthConnectors = 0.5,
                colConnectors = 'black',
                xlim = c(-12,12),
                ylim = c(0,50),
                subtitle = NULL,
                title = NULL,
                max.overlaps = 10
                )
dev.off()  




