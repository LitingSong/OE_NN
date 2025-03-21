# conda activate Py39_R43_Ju10
# run step1-2 using Rscript

.libPaths(c('/sc/arion/projects/CommonMind/liting/conda/envs/Py39_R43_Ju10/lib/R/library',.libPaths()))
.libPaths(c("/sc/arion/projects/roussp01a/jaro/programs/R_libs_4_2", .libPaths()))
source('/sc/arion/projects/roussp01a/jaro/atacseq_ad/NYGC_AD_R01/downstream_magma.R')
magmaGeneSetResultReader=function(testID,outFile){
  print(testID)
  z=data.frame(testID=testID,read.table(outFile,comment.char="#",quote="",header=T,stringsAsFactors=F),stringsAsFactors=F)
  if("FULL_NAME" %in% colnames(z)){  #it is only in there if one ore more gene set names are truncated
    z$name=z$FULL_NAME
    z$FULL_NAME=NULL
  }else{
    if("SET" %in% colnames(z)){
      z$name=z$SET #magma v1.06
    }else{
      z$name=z$VARIABLE #magma v1.07
    }
  }
  z$SET=NULL
  z=z[order(z$P),]
  return(z)
}

doFdr=function(x){

  #library(fdrtool)
  #x$global_FDR_AdjP=fdrtool(x$P, statistic="pvalue", plot=F,cutoff.method="fndr",verbose=F)$qval
  x$global_Bonf_AdjP=p.adjust(x$P, method="bonferroni")
  x$global_BH_AdjP=p.adjust(x$P, method="BH")
  x$plotLabel=""
  x$plotLabel[x$P<0.05]="*"
  x$plotLabel[x$global_BH_AdjP<0.05]="#"

  return(x)
}


# gene cluster or gene trend


clusterx <- c('magma_graph_trend','magma_graph_cluster',"magma_tf_targets",'magma_de_overlap','magma_demerge_overlap')[3]

testInfo <- read.delim(paste0('/sc/arion/projects/CommonMind/liting/ENT/',clusterx,'/meta-files/testInfo.tsv'))

geneSetResults=do.call(rbind,mapply(magmaGeneSetResultReader,testInfo$testID,testInfo$outFile,SIMPLIFY=F,USE.NAMES=F))

#add most relvant cols from testInfo
geneSetResults=cMerge(geneSetResults,testInfo[,colnames(testInfo) %in% c("geneMetaSets","geneTypeName","gwasAcronym","gwasTrait", "testID"),drop=F],"testID")

#convert gwas acronyms to ordered factors for proper plotting
geneSetResults$gwasAcronym=ordered(geneSetResults$gwasAcronym,levels=unique(geneSetResults$gwasAcronym))
geneSetResults$gwasTrait=ordered(geneSetResults$gwasTrait,levels=unique(geneSetResults$gwasTrait))


geneSetGroupAggOverview=unique(geneSetResults[,c("geneTypeName","geneMetaSets")])
geneSetGroupAggOverview$aggID=make.names(apply(geneSetGroupAggOverview,1,paste,collapse=" "))
if(any(duplicated(geneSetGroupAggOverview$aggID))) stop("unfortunate combination of geneTypeName and geneMetaSets led to duplicated aggID")
rownames(geneSetGroupAggOverview)=geneSetGroupAggOverview$aggID



# file filter
# merged dataset

doGeneSetGroupAgg=function(aggID,geneTypeName,geneMetaSets){
  z=geneSetResults[geneSetResults$geneTypeName==geneTypeName & geneSetResults$geneMetaSets==geneMetaSets   ,] # & grepl('inner|outer',geneSetResults$VARIABLE)
  z$aggID=aggID
  z=doFdr(z)
  return(z)
}
geneSetGroupAgg=mapply(doGeneSetGroupAgg,geneSetGroupAggOverview$aggID,geneSetGroupAggOverview$geneTypeName,geneSetGroupAggOverview$geneMetaSets,SIMPLIFY=F)

geneSetGroupAgg_rs <- geneSetGroupAgg$ensemblProtCodGenes35kb10kbAutosomesNoBmhc.myGeneSets
sig_gwasTrait <- subset(geneSetGroupAgg_rs,P < 0.05 )[,'gwasTrait']
#sig_gwasTrait <- subset(geneSetGroupAgg_rs,global_BH_AdjP < 0.1)[,'gwasTrait']
#geneSetGroupAgg_rs_sig <- subset(geneSetGroupAgg_rs,gwasAcronym%in%unlist(mytraits))
geneSetGroupAgg_rs_sig <- subset(geneSetGroupAgg_rs, gwasTrait%in%sig_gwasTrait)

library(RColorBrewer)
myPalette = colorRampPalette(brewer.pal(9, "Greens"), space="Lab")
#myPalette2way=colorRampPalette(rev(c("#67001F","#B2182B","#D6604D","#F4A582","#FDDBC7","#F7F7F7","#D1E5F0","#92C5DE","#4393C3","#2166AC","#053061")), space="Lab")


ggplot(geneSetGroupAgg_rs_sig, aes(x=gwasTrait, y=VARIABLE))+
  geom_tile(aes(fill=-log10(P))) +
  scale_fill_gradientn(colours = myPalette(100))+
  geom_text(aes(label=plotLabel),alpha=0.7)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+ylab('Trend pattern')+
  ggtitle('graph')

# #plotAndAggMagma(myInput='/sc/arion/projects/CommonMind/liting/ENT/magma',outDir='/sc/arion/projects/CommonMind/liting/ENT/magma', calcFdr=T )
#
#
# dev.print(pdf, file='/sc/arion/projects/roussp01a/liting/Olf/figures/magma_trend.pdf')
#


# # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#

# mytraits=list(`Psychiatric`=c("sz3",'asd',"ptsd_1",'adhd_ipsych',"mdd2",'asd_ipsych', "bip",'insomn2', "bdAndSczShared", 'mdd_without_23andMe',"adhd", "neu2"),
#               Behavioral=c("drinking", "smoking",'eduAttainment','alcoholism_2019','child_intel','gaming'),
#               Neurological=c('alz2','als','ms','pd_without_23andMe','stroke'),
#               `Non-brain related`=c('bmi','height',"covid19_A1", "covid19_A2", "covid19_B1", "covid19_B2")
# )


# mytraits=list(`Psychiatric`=c('scz_mvp4_pgc3_eur_afr','adhd_ipsych',"mdd2",'asd_ipsych', "bip_koromina_2024",'insomn2', "bdAndSczShared"),
#               Behavioral=c("drinking", "smoking",'eduAttainment','alcoholism_2019',"intel",'gaming'),
#               Neurological=c('alz2','als','ms','pd_without_23andMe','stroke')
#               #`Non-brain related`=c('bmi','height',"covid19_A1", "covid19_A2", "covid19_B1", "covid19_B2")
# )


mytraits=list(`Psychiatric`=c('scz_mvp4_pgc3_eur','adhd_ipsych',"mdd_ipsych",'asd_ipsych', "bip_koromina_2024",'insomn2', "bdAndSczShared"),
              Behavioral=c("drinking", "smoking",'eduAttainment','alcoholism_2019',"intel",'gaming'),
              Neurological=c('alzBellenguez','als2021','ms','pd_without_23andMe','stroke')#,
              #`Non-brain related`=c('bmi','height',"covid19_A1", "covid19_A2", "covid19_B1", "covid19_B2")
)

mytraits_df <- data.frame()
for (term in names(mytraits)){
  disease <- mytraits[[term]]
  df <- cbind(rep(term, length(disease)),disease)
  mytraits_df <- rbind( mytraits_df,df)
}

rownames(mytraits_df) <- mytraits_df$disease


geneSetGroupAgg_rs_sig <- subset(geneSetGroupAgg_rs,gwasAcronym%in%unlist(mytraits))

library(RColorBrewer)
myPalette = colorRampPalette(brewer.pal(9, "Greens"), space="Lab")
#myPalette2way=colorRampPalette(rev(c("#67001F","#B2182B","#D6604D","#F4A582","#FDDBC7","#F7F7F7","#D1E5F0","#92C5DE","#4393C3","#2166AC","#053061")), space="Lab")

geneSetGroupAgg_rs_sig$category <- factor(mytraits_df[as.character(geneSetGroupAgg_rs_sig$gwasAcronym),"V1"],levels = names(mytraits))
geneSetGroupAgg_rs_sig <- doFdr(geneSetGroupAgg_rs_sig)


geneSetGroupAgg_rs_sig$VARIABLE <- gsub('X','',geneSetGroupAgg_rs_sig$VARIABLE)
#save(geneSetGroupAgg_rs_sig, file='/sc/arion/projects/roussp01a/liting/Olf/data/magma_TF_tragets.RData')


ggplot(geneSetGroupAgg_rs_sig, aes(x=str_split(gwasAcronym,'_',simplify = T)[,1] , 
                                   y=VARIABLE
                                   #y=factor( str_to_title(VARIABLE), levels =str_to_title(c('down','trans_up','up','trans_down')))
))+
  geom_tile(aes(fill=-log10(P))) +theme_bw()+
  scale_fill_gradientn(colours = myPalette(100))+
  geom_text(aes(label=plotLabel),alpha=0.7)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position='top')+
  ylab('Gene clusters')+xlab('')+
  facet_grid(.~category,scales = 'free_x',space='free')

ggplot(geneSetGroupAgg_rs_sig, aes(x=str_split(gwasAcronym,'_',simplify = T)[,1] , 
                                   y=VARIABLE
                                  #y=factor( str_to_title(VARIABLE), levels =str_to_title(c('down','trans_up','up','trans_down')))
))+
  geom_tile(aes(fill=BETA)) +theme_bw()+
  scale_fill_gradient2()+
  geom_text(aes(label=plotLabel),alpha=0.7)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position='top')+
  ylab('Gene clusters')+xlab('')+
  facet_grid(.~category,scales = 'free_x',space='free')

dev.print(pdf, file='/sc/arion/projects/roussp01a/liting/Olf/figures/magma_trend_bytraitcat.pdf', width=5, height=3.5)

