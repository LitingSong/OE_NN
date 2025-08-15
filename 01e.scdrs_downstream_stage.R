
setwd('/sc/arion/projects/roussp01a/liting/Olf/figures/downstream_iOSN/')
scdrs_outps <-  list.files('./', pattern = 'scdrs_group.N_types_stage')

cor_traits <- c()
for (trait in scdrs_outps){
  cor_trait <- read.table(trait, header = T)
  cor_trait$trait <- trait
  cor_traits <- rbind(cor_traits, cor_trait )
}


mytraits=list(`Psychiatric`=c('scz_mvp4_pgc3_eur','adhd_ipsych',"mdd_ipsych",'asd_ipsych', "bip_koromina_2024"),
              Neurological=c('alzBellenguez','als2021','ms','pd_without_23andMe','stroke')

)

cor_traits$Trait <- gsub(".scdrs_group.N_types_stage", '', cor_traits$trait)
cor_traits <- subset(cor_traits, Trait%in% unlist(mytraits))

cor_traits$assoc_mcfdr <- p.adjust(cor_traits$assoc_mcp,method = 'fdr')
cor_traits$hetero_mcfdr<- p.adjust(cor_traits$hetero_mcp,method = 'fdr')

cor_traits$plotLabel <- ifelse(cor_traits$assoc_mcfdr < 0.1,'#',ifelse(cor_traits$assoc_mcp < 0.1,'x','')) 
cor_traits$trait <- gsub('.scdrs_group.N_types_stage','',cor_traits$trait)
cor_traits_sig <- subset(cor_traits, trait%in%subset(cor_traits, assoc_mcp < 0.1)$trait)
cor_traits_sig$plotLabel2 <- ifelse(cor_traits_sig$assoc_mcp < 0.1 & cor_traits_sig$hetero_mcp < 0.1, 'x','' ) 


mytraits_df <- data.frame()
for (term in names(mytraits)){
  disease <- mytraits[[term]]
  df <- cbind(rep(term, length(disease)),disease)
  mytraits_df <- rbind( mytraits_df,df)
}

rownames(mytraits_df) <- mytraits_df$disease

select_trait <- subset(cor_traits,trait%in%mytraits_df$disease)
select_trait$category <- factor(mytraits_df[select_trait$trait,'V1'], levels = names(mytraits))
  
library(RColorBrewer)
myPalette = colorRampPalette(brewer.pal(9, 'RdYlBu'), space="Lab")
library(stringr)
select_trait$Traits <- gsub('\\d','',str_split(select_trait$trait,'_',simplify = T)[,1])
select_trait$Traits <- str_to_title(select_trait$Traits)
select_trait$Traits[select_trait$Traits=='Alzbellenguez'] <- 'AD'


select_trait$Traits <- ifelse(select_trait$Traits%in%str_to_title(c('adhd','asd','bip','mdd','scz','als','ms','pd')),toupper(select_trait$Traits),select_trait$Traits)
ggplot(select_trait, aes(x= Traits, y= factor(group, levels = c("GBC",    "e_iOSN", "l_iOSN" ,"mOSN"  ))))+
#ggplot(subset(cor_traits_sig), aes(x=trait, y=group))+  
  geom_tile(aes(fill=assoc_mcz)) + 
  #scale_fill_distiller(palette = "RdBu")+
  scale_fill_gradient2(low='#2765A2',high='#A40F14')+
  geom_text(aes(label=plotLabel),alpha=0.7)+
  theme_bw()+ylab('Cell types')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),legend.position ='bottom',
        #plot.caption = element_text(hjust = 1,vjust = 1)
        )+ facet_grid(~category,scales = 'free_x',space='free')+
  ylab('Cell types')+xlab('')+labs(fill='z-score') #caption = "*: p < 0.05",

#dev.print(pdf, file='/sc/arion/projects/roussp01a/liting/Olf/figures/scdrs_bycelltypestage.pdf', width=6, height=4.5)
dev.print(pdf, file='/sc/arion/projects/roussp01a/liting/Olf/figures/scdrs_bycelltypestage_fdr.pdf', width=5, height=4)


