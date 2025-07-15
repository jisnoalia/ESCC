##Extended Data Fig. 2-----
#Extended Data Fig. 2 a 
Idents(sc)<-sc$L3_C;my<-subset(sc,ident=c('CD4'));Idents(my)<-my$L4_C
T_plotG <- c('CD55','CCR7','ABLIM1','LEF1','KLF2','TOX2','PDCD1','TOX','IL6ST','CXCR5','IL22','IL17A','IL17F','CCR6','RORA',
             'CCL5','GZMK','NKG7','PRF1','CTSW','ISG15','IFI6','IFI44L','IFIT1','IRF7','TCF7','TXNIP','MYC','NOSIP','TLE5','ANXA1','FTH1','MYADM','IL7R','RGCC',
             'CD69','PTGER4','SLC2A3','BTG2','GPR183','SESN3','FCRL3','RTKN2','IKZF2','TGIF1','HLA-DRB1','FOXP3','HLA-DRA',
             'GBP5','IL2RA','TNFRSF4','TNFRSF18','LAIR2','CTLA4','TNFRSF9','ITGB1','AQP3','TRADD','CDC25B','SH3BP5','CXCL13',
             'IFNG','IL21','BHLHE40','GADD45G')
celltype <- c('CD4_C1_CCR7', 'CD4_C10_CXCR5', 'CD4_C11_IL17A','CD4_C12_NKG7', 'CD4_C13_ISG15','CD4_C2_TCF7',  
              'CD4_C3_ANXA1', 'CD4_C4_GPR183','CD4_C5_RTKN2','CD4_C6_CD25','CD4_C7_OX40','CD4_C8_ITGB1','CD4_C9_IFNG')
ctidx <- c(5,5,5,5,5,5,5,5,5,5,5,5,5)
getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))
dotplot(my,T_plotG,celltype,ctidx,color = rev(getPalette(10)))

#Extended Data Fig. 2 b 
Idents(sc)<-sc$L3_C;my<-subset(sc,ident=c('CD8'));Idents(my)<-my$L4_C
T_plotG <- c('CCR7','LEF1','SELL','RPS13','RPL32','CCL4','FCGR3A','CCL3','SPON2','BHLHE40','ISG15','IFI44L','IFIT1','IFIT3','OAS1','KLRB1','NCR3','ZBTB16',
             'SLC4A10','TRAV1-2','IL7R','RGCC','SLC2A3','ZNF331','PTGER4','LTB','TCF7','TXNIP','LDLRAP1','FLT3LG','LMNA','ANXA1','CRIP1','CAPG','PLP2',
             'GZMK','TNFSF9','SH2D1A','CRTAM','TRAT1','CXCL13','ENTPD1','TIGIT','HAVCR2','CTLA4','CX3CR1','ADGRG1','PLEK','S1PR5',
             'FCRL6','KIR3DL2','KIR2DL4','ZNF683','KIR2DL3','XCL1','TYROBP','KLRC2','KLRF1','KLRC3','CD160')
celltype <- c('CD8_C1_CCR7',  'CD8_C10_CD16', 'CD8_C11_ISG15',  'CD8_C12_CD161','CD8_C2_IL7R',   'CD8_C3_TCF7',  'CD8_C4_ANXA1',   'CD8_C5_GZMK',
              'CD8_C6_CD39','CD8_C7_CX3CR1','CD8_C8_KIR',  'CD8_C9_KLRC2')
ctidx <- c(5,5,5,5,5,5,5,5,5,5,5,5)
getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))
dotplot(my,T_plotG,celltype,ctidx,color = rev(getPalette(10)))

#Extended Data Fig. 2 c
Idents(sc)<-sc$L1_C;my<-subset(sc,ident=c('B','Plasma'));Idents(my)<-my$L4_C
T_plotG <- c('NME1','EIF4A1','CCR7','EIF5A','CRIP1','DUSP4','ITGB1','ANXA2','HSPA1A','GPR183','FOSB','KLF6','TCL1A','FCER2',
             'IGHD','FCRL1','IFITM1','ISG15','IFI6','IFIT2','BCL6','BCL7A','CD38','SERPINA9','STMN1','HMGB2','UBE2C','MKI67',
             'JCHAIN','IGHG1','MZB1','SSR4')
celltype <- c('B_C1_CCR7','B_C2_DUSP4','B_C3_GPR183','B_C4_TCL1A','B_C5_ISG15', 'B_C6_BCL6','B_C7_MKI67','Plasma')
ctidx <- c(4,4,4,4,4,4,4,4)
getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))
dotplot(my,T_plotG,celltype,ctidx,color = rev(getPalette(10)))

#Extended Data Fig. 2 d
my<-subset(sc,L1_C=='Endothelium');Idents(my)<-my$L4_C
T_plotG <- c('ACKR1','SELP','SELE','CCL23','SEMA3G','FBLN5','GJA5','SERPINE2','PLVAP','RGCC','APLNR','IGFBP5',
             'CCL21','PROX1','PDPN','FLT4')
celltype <- c('Endo_C1_ACKR1','Endo_C2_FBLN5','Endo_C3_RGCC','Endo_C4_CCL21')
ctidx <- c(4,4,4,4)
getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))
dotplot(my,T_plotG,celltype,ctidx,color = rev(getPalette(10)))

#Extended Data Fig. 2 e
Idents(sc)<- sc$L2_C
my<-subset(sc,ident=c('Myeloid'));Idents(my)<-my$L4_C
T_plotG <- c('CLEC9A','XCR1','CPNE3','SNX3','CD1C','FCER1A','CLEC10A','HLA-DQA1','CD1A','CD1E','LTB','CD207','LAMP3','CCR7',
             'CCL19','EBI3','IL3RA','GZMB','JCHAIN','LILRA4','NLRP3','EREG','IL1B','AREG','SPP1','MMP12','MMP9','CSTB','TREM2',
             'APOC1','C1QC','LIPA','LYVE1','CCL18','FOLR2','CCL2','CXCL10','CXCL9','ISG15','GBP1','MKI67','STMN1','TUBA1B','HMGN2',
             'S100A9','S100A8','FCN1','VCAN','RHOC','LST1','MS4A7','LILRB2','CXCL8','G0S2','FCGR3B','CXCR2')
celltype <- c('DC_C1_CLEC9A','DC_C2_CD1C','DC_C3_CD1A', 'DC_C4_LAMP3','DC_C5_IL3RA','Mac_C1_NLRP3', 'Mac_C2_SPP1','Mac_C3_C1QC',
              'Mac_C4_LYVE1', 'Mac_C5_CXCL10','Mac_C6_MKI67', 'Mono_C1_CD14', 'Mono_C2_CD16', 'Neutrophil')
ctidx <- c(4,4,4,4,4,4,4,4,4,4,4,4,4,4)
getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))
dotplot(my,T_plotG,celltype,ctidx,color = rev(getPalette(10)))
#Extended Data Fig. 2 f
my<-subset(sc,L1_C=='Fibroblast');Idents(my)<-my$L4_C;my<-subset(my,ident='Smooth muscle cell',invert=T)
T_plotG <- c("CFD","IGFBP6","GPX3","PI16","IGF1","DPT","CXCL12","COL4A4","MMP11","POSTN","COL1A1","CTHRC1","APOD","COL13A1",
             "CXCL14","CXCL1","HIGD1B","COX4I2","KCNJ8","PDGFRB","MUSTN1","ADIRF","ACTA2","ADAMTS9")
celltype <- c('FB_C1_CFD','FB_C2_IGF1','FB_C3_COL1A1','FB_C4_APOE','FB_C5_PDGFRB','FB_C6_ACTA2')
ctidx <- c(4,4,4,4,4,4)
getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))
dotplot(my,T_plotG,celltype,ctidx,color = rev(getPalette(10)))