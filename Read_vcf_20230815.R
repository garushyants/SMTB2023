#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# args = c("ZI_allchrom_clean_bi_correctANN.ann.pick.100000.vcf",
#          "test2",
#          "/home/garushyants/Documents/Teaching/SMTB2023/20230815_day3")

library(tidyr)
library(dplyr)
library(ggplot2)

if (length(args)==0) {
  stop("At least two arguments must be supplied (input file).n",
       call.=FALSE)
} else if(length(args)==1){
  args[2]="tmp"
  args[3]="./"
} else if(length(args)==2){
  args[3]="./"
}

setwd(args[3])

#create folder
if (!dir.exists(args[2])){
  dir.create(args[2])
}else{
  print("dir exists")
}


TestVCF<-read.csv(args[1], comment.char = "#", sep="\t", header = F)

# VCFSep<-separate(TestVCF,V8, into=c("AC","AN","ANN"), sep=';')
#This was modified in order to better preserve all the information in original VCF file
VCFSep<-TestVCF %>% 
  separate(
    V8, into = c("AC","AN","ANN"), sep = ";", extra = "merge"
  )

VCFSep$ACmod<-gsub("AC=","",VCFSep$AC)
VCFSep$ANmod<-gsub("AN=","",VCFSep$AN)
dim(VCFSep)
VCFSepFilt<-VCFSep %>% filter(!grepl(",",V5))
dim(VCFSepFilt)
VCFSepF<-separate(VCFSepFilt,ACmod,sep=",",into=c("ACm","N"))
VCFSepF[is.na(VCFSepF)]<-0
VCFSepF$NSeqFlies<-as.numeric(VCFSepF$ANmod) - as.numeric(VCFSepF$N)


###Filtering procedure

#Filter by the number of sequenced flies in each position
PlotSeqFlies<-ggplot(VCFSepF)+
  geom_bar(aes(x=NSeqFlies))+
  theme_bw()+
  xlab("Number of flies with sequenced position")
PlotSeqFlies
#look at quantiles
# PossibleCutoffs<-quantile(VCFSepF$NSeqFlies, probs = c(0.05, 0.1, 0.25), na.rm = FALSE,
#          names = TRUE)
FixedCutOff<-180

PlotSeqFliesToSave<-PlotSeqFlies+
  geom_vline(xintercept=FixedCutOff,
             color="red")
PlotSeqFliesToSave

#
VCFSepFilterApplied<-subset(VCFSepF,VCFSepF$NSeqFlies>FixedCutOff)
##################
#Get only singletons
plot2<-ggplot(VCFSepFilterApplied)+
  geom_bar(aes(x=as.numeric(ACm)),
           width =1)+
  xlab("# of samples with AF")+
  theme_bw()
plot2
#Do we need other allele frequencies??

print("Number of singletons in data")
dim(subset(VCFSepFilterApplied,as.numeric(VCFSepFilterApplied$ACm)==1))[1]
#
VCFDataToSaveSingletons<-subset(VCFSepFilterApplied,
                                VCFSepFilterApplied$ACm==1)[,c(1:10)]

##########################
#Save plots
ggsave("Plot1.png",PlotSeqFliesToSave,path=args[2])
ggsave("Plot2.png",plot2,path=args[2])

#Save VCF without header

outfilename<-paste(gsub(".vcf","",args[1]),"_singletons_",FixedCutOff,".vcf",
                               sep="")

colnames(VCFDataToSaveSingletons)<-c('##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=.,Type=String,Description="Genotype">
##INFO=<ID=AC,Number=2,Type=Integer,Description="Allele count in genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##contig=<ID=2L,length=23011544>
##contig=<ID=2R,length=21146708>
##contig=<ID=3L,length=24543557>
##contig=<ID=3R,length=27905053>
##contig=<ID=X,length=22422827>
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: \'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO\'">
##INFO=<ID=LOF,Number=.,Type=String,Description="Predicted loss of function effects for this variant. Format: \'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected\'">
##INFO=<ID=NMD,Number=.,Type=String,Description="Predicted nonsense mediated decay effects for this variant. Format: \'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected\'">
##bcftools_viewVersion=1.18+htslib-1.18
##bcftools_viewCommand=view -h ZI_allchrom_clean_bi_correctANN.ann.vcf; Date=Mon Aug  7 21:02:02 2023
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ZI103',rep("",9))

write.table(VCFDataToSaveSingletons,
            file=outfilename,
            sep="\t",
            quote=F,
            col.names = T,
            row.names = F)
