#!/usr/bin/Rscript



#(head -n 1 SNP_replace_BBID_table.tsv; tail -n +2 SNP_replace_BBID_table.tsv | perl -lane '$F[0]=~s/\.\d+$//;for ($i=1;$i<scalar(@F);$i++){$F[$i]=~s/:\d+$//;@arr=split(/\//,$F[$i]);print STDERR "Error: heterzygous site" unless ($arr[0]==$arr[1]);$F[$i]=$arr[0];} print join("\t",@F);') > SNP_replace_BBID_table.clean.tsv







filename="20251215.BB_unique_replace01.csv"
filename="SNP_replace_BBID_table.clean.tsv"



### load packages
install.packages("eulerr")
library(eulerr)


### read data
data<-read.table(filename, header=T)
dim(data)
head(data)
data<-data[which(data$BB==0 & data$SS<2 & data$St<2 & data$EE<2 & data$RR<2 & data$HH<2 &data$TT <2),c("CHROM","BB","SS","St","EE","RR","HH", "TT")]
dim(data)
#Total SNPs: 983718
head(data)
nrow(data[which(data$HH==0 | data$St==0 | data$RR==0 | data$EE==0 | data$TT==0 | data$SS==0),])
#Total SNPs: 924993 (94.03%)


unique(data$BB)
unique(data$SS)
unique(data$St)
unique(data$EE)
unique(data$RR)
unique(data$HH)
unique(data$TT)

### number of sequences
length(unique(data$CHROM))
#10439






table1<-data.frame(Species=factor(c("HH","St","RR","EE","TT")),SNPs=c(
	nrow(data[which(data$BB==0 & data$SS==1 & data$HH==0),]),
	nrow(data[which(data$BB==0 & data$SS==1 & data$St==0),]),
	nrow(data[which(data$BB==0 & data$SS==1 & data$RR==0),]),
	nrow(data[which(data$BB==0 & data$SS==1 & data$EE==0),]),
	nrow(data[which(data$BB==0 & data$SS==1 & data$TT==0),]))
)
#HH 82171
#St 82822
#RR 82005
#EE 82260
#TT 81035
#table1<-data.frame(Species=c("HH","St","RR","EE","TT"),SNPs=c(82171,82822,82005,82260,81035))
table1$Species<-factor(table1$Species, levels = c("HH","St","RR","EE","TT"))
plot_noSS<- ggplot(data = table1, aes(x = Species, y = SNPs)) +
  geom_col(colour = "#5C5C5C", size = 0.25, width=0.5,
           fill = "#5C5C5C", alpha = 1) +
  ylab("BB SNPs") +
  coord_cartesian(ylim=c(80000,83500)) +
  theme(
    axis.title = element_text(size = 10, family = "ArialMT", color = "black"),
    axis.text = element_text(size = 8, family = "ArialMT", color = "black")
)
plot_noSS
ggsave ("St_explains_more_BB_allele.svg", width=8,height=8,units="cm",dpi=600)



table2<-data.frame(Species=c("HH","St","RR","EE","TT","SS"),
	SNPs=c(nrow(data[which(data$BB == 0 & data$HH == 0), ]),
		   nrow(data[which(data$BB == 0 & data$St == 0), ]),
		   nrow(data[which(data$BB == 0 & data$RR == 0), ]),
		   nrow(data[which(data$BB == 0 & data$EE == 0), ]),
		   nrow(data[which(data$BB == 0 & data$TT == 0), ]),
		   nrow(data[which(data$BB == 0 & data$SS == 0), ])
		)
	)
table2$Species<-factor(table2$Species, levels = c("HH","St","RR","EE","TT","SS"))
#517127 646601 643345 708397 750014 813948
plot_SS<- ggplot(data = table2, aes(x = Species, y = SNPs)) +
  geom_col(colour = "#5C5C5C", size = 0.75, width=0.5,
           fill = "#5C5C5C", alpha = 1) +
  ylab("BB SNPs") +
  coord_cartesian(ylim=c(500000,820000)) +
  theme(
    axis.title = element_text(size = 10, family = "ArialMT", color = "black"),
    axis.text = element_text(size = 8, family = "ArialMT", color = "black")
)

library(patchwork)
patchwork <- plot_SS + plot_noSS
patchwork <- patchwork + plot_annotation(tag_levels = 'A')
ggsave ("St_explains_more_BB_allele.17cm.svg", plot=patchwork, width=17,height=8,units="cm",dpi=600)
ggsave ("St_explains_more_BB_allele.8cm.svg", plot=patchwork, width=8,height=8,units="cm",dpi=600)



### plot1: BB St SS

data1<-data[which((data$SS + data$St)==1 & data$BB==0), c("CHROM","BB", "SS", "St")]
dim(data1)
head(data1)
length(unique(data1$CHROM))

BBSSSt1<-nrow(data1[which(data1$BB==0 & data1$SS==0 & data1$St==0),])
BBSS1<-nrow(data1[which(data1$BB==0 & data1$SS==0 & data1$St==1),])
BBSt1<-nrow(data1[which(data1$BB==0 & data1$SS==1 & data1$St==0),])
SS1<-nrow(data1[which(data1$BB==0 & data1$SS==1),])
St1<-nrow(data1[which(data1$BB==0 & data1$St==1),])
BB1<-nrow(data1[which(data1$BB==0),])-BBSSSt1-BBSS1-BBSt1


svg("euler_plot.BBspcific.BEHRT.svg", width = 8.5, height = 8.5, family = "ArialMT")
num1<-c("BB"=BB1, "SS"=SS1, "St"=St1, "BB&SS"=BBSS1, "BB&St"=BBSt1, "BB&SS&St"=BBSSSt1)
#plot(euler(num1, shape = "ellipse"), quantities = c(BB1,SS1,St1, BBSS1,BBSt1,BBSSSt1))
plot(euler(num1, shape = "ellipse"), quantities = TRUE)



### plot2: BB unique gene: BB EE HH RR TT
data2<-data[which(data$BB==0 & data$SS==1 & data$St==1), c("BB", "EE", "RR", "HH", "TT")]
dim(data2)
head(data2)


BEHRT<-nrow(data2[which(data2$BB==0 & data2$EE==0 & data2$HH==0 & data2$RR==0 & data2$TT==0),])
BEHR<-nrow(data2[which(data2$BB==0 & data2$EE==0 & data2$HH==0 & data2$RR==0),])-BEHRT
BEHT<-nrow(data2[which(data2$BB==0 & data2$EE==0 & data2$HH==0 & data2$TT==0),])-BEHRT
BERT<-nrow(data2[which(data2$BB==0 & data2$EE==0 & data2$RR==0 & data2$TT==0),])-BEHRT
BHRT<-nrow(data2[which(data2$BB==0 & data2$HH==0 & data2$RR==0 & data2$TT==0),])-BEHRT
EHRT<-nrow(data2[which(data2$EE==0 & data2$HH==0 & data2$RR==0 & data2$TT==0),])-BEHRT
BEH<-nrow(data2[which(data2$BB==0 & data2$EE==0 & data2$HH==0),])-BEHR-BEHT-BEHRT
BER<-nrow(data2[which(data2$BB==0 & data2$EE==0 & data2$RR==0),])-BEHR-BERT-BEHRT
BET<-nrow(data2[which(data2$BB==0 & data2$EE==0 & data2$TT==0),])-BERT-BEHT-BEHRT
BHR<-nrow(data2[which(data2$BB==0 & data2$HH==0 & data2$RR==0),])-BEHR-BHRT-BEHRT
BHT<-nrow(data2[which(data2$BB==0 & data2$HH==0 & data2$TT==0),])-BEHT-BHRT-BEHRT
BRT<-nrow(data2[which(data2$BB==0 & data2$RR==0 & data2$TT==0),])-BERT-BHRT-BEHRT
EHR<-nrow(data2[which(data2$EE==0 & data2$HH==0 & data2$RR==0),])-BEHR-EHRT-BEHRT
EHT<-nrow(data2[which(data2$EE==0 & data2$HH==0 & data2$TT==0),])-BEHT-EHRT-BEHRT
ERT<-nrow(data2[which(data2$EE==0 & data2$RR==0 & data2$TT==0),])-BERT-EHRT-BEHRT
HRT<-nrow(data2[which(data2$HH==0 & data2$RR==0 & data2$TT==0),])-BHRT-EHRT-BEHRT
BE<-nrow(data2[which(data2$BB==0 & data2$EE==0),])-BEH-BER-BET-BEHR-BEHT-BERT-BEHRT
BH<-nrow(data2[which(data2$BB==0 & data2$HH==0),])-BEH-BHR-BHT-BEHR-BEHT-BHRT-BEHRT
BR<-nrow(data2[which(data2$BB==0 & data2$RR==0),])-BER-BHR-BRT-BEHR-BERT-BHRT-BEHRT
BT<-nrow(data2[which(data2$BB==0 & data2$TT==0),])-BET-BHT-BRT-BEHT-BERT-BHRT-BEHRT
EH<-nrow(data2[which(data2$EE==0 & data2$HH==0),])-BEH-EHR-EHT-BEHR-BEHT-EHRT-BEHRT
ER<-nrow(data2[which(data2$EE==0 & data2$RR==0),])-BER-EHR-ERT-BEHR-BERT-EHRT-BEHRT
ET<-nrow(data2[which(data2$EE==0 & data2$TT==0),])-BET-EHT-ERT-BEHT-BERT-EHRT-BEHRT
HR<-nrow(data2[which(data2$HH==0 & data2$RR==0),])-BHR-EHR-HRT-BEHR-BHRT-EHRT-BEHRT
HT<-nrow(data2[which(data2$HH==0 & data2$TT==0),])-BHT-EHT-HRT-BEHT-BHRT-EHRT-BEHRT
RT<-nrow(data2[which(data2$RR==0 & data2$TT==0),])-BRT-ERT-HRT-BERT-BHRT-EHRT-BEHRT
B<-nrow(data2[which(data2$BB==0),])-BE-BH-BR-BT-BEH-BER-BET-BHR-BHT-BRT-BEHR-BEHT-BERT-BHRT-BEHRT
E<-nrow(data2[which(data2$EE==0),])-BE-EH-ER-ET-BEH-BER-BET-EHR-EHT-ERT-BEHR-BEHT-BERT-EHRT-BEHRT
H<-nrow(data2[which(data2$HH==0),])-BH-EH-HR-HT-BEH-BHR-BHT-EHR-EHT-HRT-BEHR-BEHT-BHRT-EHRT-BEHRT
R<-nrow(data2[which(data2$RR==0),])-BR-ER-HR-RT-BER-BHR-BRT-EHR-ERT-HRT-BEHR-BERT-BHRT-EHRT-BEHRT
TT<-nrow(data2[which(data2$TT==0),])-BT-ET-HT-RT-BET-BHT-BRT-EHT-ERT-HRT-BEHT-BERT-BHRT-EHRT-BEHRT


#c(B, E, H, R, TT, BE, BH, BR, BT, BEH, BER, BET, BHR, BHT, BRT, EHR, EHT, HRT, BEHR, BEHT, BHRT, EHRT, BEHRT)

num<-c("BB"=B, "EE"=E, "HH"=H, "RR"=R, "TT"=TT)
num<-c("BB"=B, "EE"=E, "HH"=H, "RR"=R, "SS"=S, "St"=St, "TT"=TT)
num<-c("BB"=B, "EE"=E, "HH"=H, "RR"=R, "SS"=S, "St"=St, "TT"=TT, "BB&EE"=BBEE, "BB&HH"=BBHH, "BB&RR"=BBRR, "BB&SS"=BBSS, "BB&St"=BBSt, "BB&TT"=BBTT)

num1<-c("BB"=B, "EE"=E, "HH"=H, "RR"=R, "TT"=TT, "BB&EE"=BE, "BB&HH"=BH, "BB&RR"=BR, "BB&TT"=BT, "BB&EE&HH"=BEH, "BB&EE&RR"=BER, "BB&EE&TT"=BET, "BB&HH&RR"=BHR, "BB&HH&TT"=BHT, "BB&RR&TT"=BRT, "EE&HH&RR"=EHR, "EE&HH&TT"=EHT, "HH&RR&TT"=HRT, "BB&EE&HH&RR"=BEHR, "BB&EE&HH&TT"= BEHT, "BB&HH&RR&TT"=BHRT, "EE&HH&RR&TT"=EHRT,"BB&EE&HH&RR&TT"=BEHRT)
plot(euler(num1, shape = "ellipse"), quantities = TRUE)

fit1<-plot(euler(num1, shape = "ellipse"), quantities = TRUE)

svg("euler_plot.BBspcific.BEHRT.svg", width = 8.5, height = 8.5, units="cm", family = "ArialMT")
plot(fit1, quantities = T)
dev.off()



### plot3: BBEEHHRR



#data2<-data[which(data$BB==0 & data$SS==1 & data$St==1), c("BB", "EE", "RR", "HH", "TT")]
#dim(data2)
#head(data2)


BEHR<-nrow(data2[which(data2$BB==0 & data2$EE==0 & data2$HH==0 & data2$RR==0),])
BEH<-nrow(data2[which(data2$BB==0 & data2$EE==0 & data2$HH==0),])-BEHR
BER<-nrow(data2[which(data2$BB==0 & data2$EE==0 & data2$RR==0),])-BEHR
BHR<-nrow(data2[which(data2$BB==0 & data2$HH==0 & data2$RR==0),])-BEHR
EHR<-nrow(data2[which(data2$EE==0 & data2$HH==0 & data2$RR==0),])-BEHR
BE<-nrow(data2[which(data2$BB==0 & data2$EE==0),])-BEH-BER-BEHR
BH<-nrow(data2[which(data2$BB==0 & data2$HH==0),])-BEH-BHR-BEHR
BR<-nrow(data2[which(data2$BB==0 & data2$RR==0),])-BER-BHR-BEHR
EH<-nrow(data2[which(data2$EE==0 & data2$HH==0),])-BEH-EHR-BEHR
ER<-nrow(data2[which(data2$EE==0 & data2$RR==0),])-BER-EHR-BEHR
HR<-nrow(data2[which(data2$HH==0 & data2$RR==0),])-BHR-EHR-BEHR
B<-nrow(data2[which(data2$BB==0),])-BE-BH-BR-BEH-BER-BHR-BEHR
E<-nrow(data2[which(data2$EE==0),])-BE-EH-ER-BEH-BER-EHR-BEHR
H<-nrow(data2[which(data2$HH==0),])-BH-EH-HR-BEH-BHR-EHR-BEHR
R<-nrow(data2[which(data2$RR==0),])-BR-ER-HR-BER-BHR-EHR-BEHR

num<-c("BB"=B, "EE"=E, "HH"=H, "RR"=R, "BB&EE"=BE, "BB&HH"=BH, "BB&RR"=BR, "BB&EE&HH"=BEH, "BB&EE&RR"=BER, "BB&HH&RR"=BHR, "EE&HH&RR"=EHR, "BB&EE&HH&RR"=BEHR)

fit2<-plot(euler(num, shape = "ellipse"), quantities = TRUE)

svg("euler_plot.BBspcific.BEHR.svg", width = 8.5, height = 8.5, units="cm", family = "ArialMT")
plot(fit2, quantities = TRUE)
dev.off()
