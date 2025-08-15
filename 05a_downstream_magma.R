##################################################################
#  MAGMA pipeline 
#
#
#    NOTE: 
#      Prequisites: prepared sumstats (see prepare_magma folder)
#      MAGMA doesn't output gene sets for gene sets where there are no included genes
#      CAVEAT: calcMagmaGeneSetPvals will probably not work if downstream_greatR.R has been sourced 
#      (as that will cause the aggMagmaPlotter to also source that file and overwrite fncts w. shared names)
#
#
#     annotateGenesMagma (second step)
#       annotate snps to genes
#       requires files generated under prepare_magma/ scripts
#
#
#     getGenePvals (third step)
#       calculate pvals for genes
#
#
#     calcMagmaGeneSetPvals: (fourth step)
#      returns jobarray info in case you wan't to link it to dependencies
#      this can take some time to run, so maybe run in screen on interactive node
#     
#      
#     plotAndAggMagma: (fifth step)
#       called by previous step to agg results
#       
#
#
#  #################################
#  #EXAMPLE:
#  
#  #2nd and 3rd step
#  source("/sc/arion/projects/roussp01a/jaro/atacseq_ad/NYGC_AD_R01/downstream_magma.R");
#  annotateGenesMagma();
#  getGenePvals()
#  
#  #4th and 5th step (jobs in steps above must have completed first)
#  load("/sc/arion/projects/roussp01b/resources/databases/gene-info/gene-sets-misc/gseafunctions_noChrY_noChrM_v3.Rdata");  #This is an output from offlineFisherGsea in downstream_greatR.R and has var standardGeneSets
#  source("bb/atacseq_ad/NYGC_AD_R01/downstream_magma.R")
#  myAnalysis=calcMagmaGeneSetPvals(  outDir="/sc/arion/projects/roussp01b/resources/databases/magma/v004/test/standardGeneSets_A",  geneMetaSets=standardGeneSets,  pvalInfo=grabGenePvalInfo(selectedGeneTypeName="ensemblProtCodGenes35kb10kbAutosomesNoBmhc")  )
#  
#
#  
#  #################################
#  #EXAMPLE: using hg38 gene sets but hg19 gene and variant positions #FIXME:
#      
#  #4th and 5th step (jobs in steps above must have completed first)
#  z=local({load("/sc/arion/projects/roussp01b/resources/databases/gene-info/gene-sets-misc/gseafunctions_noChrY_noChrM_v4_hg38.Rdata"); environment()})
#  standardGeneSets=z$standardFisher$standardGeneSets
#  source("bb/atacseq_ad/NYGC_AD_R01/downstream_magma.R")
#  myAnalysis=calcMagmaGeneSetPvals(  outDir="/sc/arion/projects/roussp01b/resources/databases/magma/v004/test/standardGeneSets_hg38hg19fusion",  geneMetaSets=standardGeneSets,  pvalInfo=grabGenePvalInfo(selectedGeneTypeName="ensemblProtCodGenes35kb10kbAutosomesNoBmhc")  )
#  
#  #As above but with updated gene sets
#  z=local({load("/sc/arion/projects/roussp01b/resources/databases/gene-info/gene-sets-misc/gseafunctions_noChrY_noChrM_v5_hg38.Rdata"); environment()})
#  standardGeneSets=z$standardFisher$standardGeneSets
#  source("bb/atacseq_ad/NYGC_AD_R01/downstream_magma.R")
#  myAnalysis=calcMagmaGeneSetPvals(  outDir="/sc/arion/projects/roussp01b/resources/databases/magma/v004/test/standardGeneSets_hg38hg19fusion_t002",  geneMetaSets=standardGeneSets,  pvalInfo=grabGenePvalInfo(selectedGeneTypeName="ensemblProtCodGenes35kb10kbAutosomesNoBmhc")  )
#  
#  #As above but with even further updated gene sets
#  z=local({load("/sc/arion/projects/roussp01b/resources/databases/gene-info/gene-sets-misc/gseafunctions_noChrY_noChrM_v6_hg38.Rdata"); environment()})
#  standardGeneSets=z$standardFisher$standardGeneSets
#  source("bb/atacseq_ad/NYGC_AD_R01/downstream_magma.R")
#  myAnalysis=calcMagmaGeneSetPvals(  outDir="/sc/arion/projects/roussp01b/resources/databases/magma/v004/test/standardGeneSets_hg38hg19fusion_t003",  geneMetaSets=standardGeneSets,  pvalInfo=grabGenePvalInfo(selectedGeneTypeName="ensemblProtCodGenes35kb10kbAutosomesNoBmhc")  )
#
#
#  #As above but with even further updated gene sets
#  z=local({load("/sc/arion/projects/roussp01b/resources/databases/gene-info/gene-sets-misc/gseafunctions_noChrY_noChrM_v8_hg38.Rdata"); environment()})
#  standardGeneSets=z$standardFisher$standardGeneSets
#  source("bb/atacseq_ad/NYGC_AD_R01/downstream_magma.R")
#  myAnalysis=calcMagmaGeneSetPvals(  outDir="/sc/arion/projects/roussp01b/resources/databases/magma/v004/test/standardGeneSets_hg38hg19fusion_t004",  geneMetaSets=standardGeneSets,  pvalInfo=grabGenePvalInfo(selectedGeneTypeName="ensemblProtCodGenes35kb10kbAutosomesNoBmhc")  )
#
#  
#  #As above but with even further updated gene sets including syngo pruned
#  z=local({load("/sc/arion/projects/roussp01b/resources/databases/gene-info/gene-sets-misc/gseafunctions_noChrY_noChrM_v9_hg38.Rdata"); environment()})
#  standardGeneSets=z$standardFisher$standardGeneSets
#  source("bb/atacseq_ad/NYGC_AD_R01/downstream_magma.R")
#  myAnalysis=calcMagmaGeneSetPvals(  outDir="/sc/arion/projects/roussp01b/resources/databases/magma/v004b/test/standardGeneSets_hg38hg19fusion_t004b",  geneMetaSets=standardGeneSets,  pvalInfo=grabGenePvalInfo(selectedGeneTypeName="ensemblProtCodGenes35kb10kbAutosomesNoBmhc")  )
#  
#
#  
#  #################################
#  #################################
#  ### OLD EXAMPLES NOW MOVED TO V001 using code prior to 2019-february commits:
#  ###    
#  ####source("bb/atacseq_ad/NYGC_AD_R01/downstream_magma.R");annotateGenesMagma();getGenePvals() #pre-computed
#  ###load("data/gene-info/gene-sets-misc/gseafunctions_noChrY_noChrM.Rdata");  #if you wan't to load the standard gene sets. This is an output from offlineFisherGsea downstream_greatR.R and has var standardGeneSets. You can have a look at how the data structure is for the gene sets 
#  ###source("bb/atacseq_ad/NYGC_AD_R01/downstream_magma.R")
#  ###myAnalysis=calcMagmaGeneSetPvals(  outDir="test/magma2/",  geneMetaSets=standardGeneSets,  pvalInfo=grabGenePvalInfo(selectedGeneTypeName="ensemblProtCodGenes35kb10kbAutosomesNoBmhc")  )
#  ###
#  ###load("data/gene-info/gene-sets-misc/gseafunctions_noChrY_noChrM_v3.Rdata"); 
#  ###source("bb/atacseq_ad/NYGC_AD_R01/downstream_magma.R")
#  ###myAnalysis3=calcMagmaGeneSetPvals(  outDir="test/magma3/",  geneMetaSets=standardGeneSets,  pvalInfo=grabGenePvalInfo(selectedGeneTypeName="ensemblProtCodGenes35kb10kbAutosomesNoBmhc")  )
#  
#  
#
#

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
#######
#######     GENERIC HELPER FUNCTIONS
#######
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

#Info for copying source to analysis dir:
if(!exists("sourcedFiles")){   primarySourceDir=dirname(normalizePath(sys.frame(1)$ofile));   sourcedFiles=normalizePath(sys.frame(1)$ofile) }

createScript=function(scriptString,
  settings,
  outDir="", 
  scriptName=NA,
  mem=20000,
  cores=1,
  queue="premium",
  walltime="3:00",
  clusterSection="bode", #this is currently disabled, as I have removed this line below: #BSUB -m ", gs("clusterSection"), "
  accountForJobs="acc_roussp01b",
  arrayLength=1,
  submit=T,
  dependencies=list(),
  isRscript=F,
  Rmodule="module load R/4.0.3",
  sourcePrevSourcedRfiles=T
){
  z=list()
  #settings is a list (potentially empty) that can overwrite settings for more concise coding
  gs=function(x){if(is.null(settings[[x]])){get(x)}else{settings[[x]]}}

  #if scriptName is not present, just take variable name
  if(is.na(scriptName))
    z$scriptName=deparse(substitute(scriptString))
  
  z$dryRun=gs("dryRun")

  z$scriptFile=paste0(gs("outDir"),"/scripts/",z$scriptName,".sh")

  z$scriptString=paste0("#!/bin/bash
    #BSUB -q ", queue,"
    #BSUB -P ", gs("accountForJobs"), " 
    #BSUB -n ", cores, " 
    #BSUB -R 'span[hosts=1]'
    #BSUB -R 'rusage[mem=", mem, "]'
    #BSUB -W ", walltime, "
    #BSUB -oo ", gs("outDir"), "/log/", z$scriptName, ".%I.out
    #BSUB -eo ", gs("outDir"), "/log/", z$scriptName, ".%I.err

    if [ -z \"$LSB_JOBINDEX\" ]; then LSB_JOBINDEX=1; fi #for debugging

    module purge
    module load glibmm
    module load openssl
    ")

  if(!isRscript){
    #just append provided code
    z$scriptString=paste0(z$scriptString, scriptString) 
  }else{
    #if it is an R script we put that into a seperate file
    z$RscriptFile=paste0(gs("outDir"),"/scripts/",z$scriptName,"_Rcode.R")
    #Should we source the same files that we've already sourced?
    #per default the sourced files are copied to the output dir
    if(sourcePrevSourcedRfiles & exists("sourcedFiles")){ 
      z$RscriptString=paste0("source(\"",gs("outDir"),"/scripts/",basename(sourcedFiles),"\")\n",collapse="")
    }else{
      z$RscriptString=""
    }

    #now add actual script code
    z$RscriptString=paste0(z$RscriptString,scriptString)
    write(z$RscriptString,z$RscriptFile)
    z$scriptString=paste0(z$scriptString,"\n", Rmodule, "\nRscript ", z$RscriptFile) 
  }

  write(z$scriptString,z$scriptFile)

  myDependencies=if(length(dependencies)){
    paste0(   " -w '",    paste(sapply(1:length(dependencies),function(x){paste0("numdone(", dependencies[[x]]$jobId,",*)")}),collapse=" && "),     "'"     ) #sorry for illegibility
  }else{""}

  z$submitString=paste0("cat ",z$scriptFile, "| bsub", myDependencies, " -J '", z$scriptName, "_", settings$runID ,"[1-",arrayLength ,"]'")

  if(submit){return(submitJob(z))}else{return(z)}
}

######################
# These are not currently used
bashPaste=function(...){
  z=list()
  z$arrayLength=1
  z$text==""
  z$type="script"
  for(i in list(...)){
    if(is.list(i)){
      z$text=paste0(i$text,"\n", z$text,"$",i$name)
    }else if(is.data.frame(i)){
      for(j in colnames(i))
      z$text=paste0(z$text,j,"=$(sed \"${LSB_JOBINDEX}q;d\" <<< '", paste(i[,j],collapse="\n"), "')\n")
    }
  }
}
bsV=function(x){ paste(mapply(function(x,y)paste0(x,"=\"",y,"\""),sapply(match.call(expand.dots=TRUE)[-1], deparse),list(...)),collapse="\n") }
bV=function(x){list(name=deparse(substitute(x)),varString=paste0(deparse(substitute(x)),"=\"",x,"\""), type="bashVar")}



######################
submitJob=function(job,verbose=T){ #"job" must be a list containing "$submitString" and "$dryRun"
  z=job
  if(!z$dryRun)
    z$systemReturnString=system(z$submitString, intern=T)
  else
    z$systemReturnString="Job <00000000> is submitted to queue <SOMETHING>"

  z$jobId=gsub("^.+<","",gsub(">.+$","",z$systemReturnString))

  if(verbose){
    print(z$submitString)
    print(z$systemReturnString)
    print(z$jobId)
  }

  return(z)
}


###################
createOutputFileStructure=function(outDir,copySourcedFiles=T){
  for(i in c("files", "log", "meta-files", "results", "scripts"))
    dir.create(paste0(outDir, "/", i),showWarnings=F,recursive=T) 
  if(exists("sourcedFiles") & copySourcedFiles)
    mapply(  file.copy,  sourcedFiles,  paste0(outDir,"/scripts/",basename(sourcedFiles)),  MoreArgs=list(overwrite=T)  )
}

###################
fixRelativePath=function(file,check=F){
  if( ifelse(is.na(file),F,!grepl("^[\\w\\d\\-\\_\\./]+$",normalizePath(file),perl=T)) )
    myStop(paste0("Filenames cannot have weird characters or spaces in them. This variable has that ", deparse(substitute(file)), " and points to: ", file))

  if(check)
    if(!file.exists(file))
      myStop(paste0("Required variable ", deparse(substitute(file)), " points to non-existent file:", file))

  return( ifelse(is.na(file),file,normalizePath(file)) )
}


###################
myStop=function(...) eval.parent(substitute({
  debugEnv <<- as.environment(as.list(environment(), all.names=T))
  stop(paste0(...,". Script aborted. Variables saved to the 'debugEnv' environment for debugging purposes. To access the variables use the normal variable name preceeded by this and a dollar sign."),call.=F)
}))


###################
cMerge=function(df1,df2,col){cbind(df1,df2[match(df1[,col],df2[,col]),colnames(df2)!=col,drop=F])} #like built in merge fct but keeps order of rows and cols


###################
mtsv=function(x,outDir,myHeader=T,gz=F){
  myFile=paste0(outDir, "/meta-files/",  deparse(substitute(x)),".tsv",ifelse(gz,".gz",""))
  if(gz) myFile=gzfile(myFile, "w")
  write.table(x,file=myFile,sep="\t",quote=F,row.names=F,col.names=myHeader)
  if(gz) close(myFile)
}




####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
#######
#######     END OF GENERIC FUNCTIONS
#######
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

suppressMessages(library(GenomicRanges))
suppressMessages(library(ggplot2))

removeMHC=function(myGr){ #hg19, not currently used
  mhc=GRanges("chr6",IRanges(25E6,width=10E6),strand="*")
  myGr[!overlapsAny(myGr,mhc)]
}

##################################################################
# snp to gene annotation
##################################################################
#
# NOTE: remember to delete output folder before rerun
# NOTE: chr sizes are for hg19
#
# genes, and gene sets will all come in three flavors:
# - "All"
# - "AutosomesNoBmhc"
# - "NoChryNoMhc"
#
# reference panels will further come in one category
# - "Magma1kgp3"
#
# snp to gene annotation will further come in two sub-sub-categories
# - "NoPadding"
# - "Padding35us10ds"
#
# That is, snpToGene{"Magma1kgp3"} x {"All"|"AutosomesNoBmhc"|"NoChryNoMhc"} x {"NoPadding"|"Padding35us10ds"}.genes.annot

#hard links!!! #hg19!!!!
annotateGenesMagma=function(
  baseDir="/sc/arion/projects/roussp01b/resources/databases/magma/v004",
  myMagma="/sc/arion/projects/roussp01b/resources/software/magma/magma_1_10/magma",
  standardGeneFile="/sc/arion/projects/roussp01b/resources/software/magma/magma_1_08/2019-02-25/NCBI37.3.gene.loc",
  lincGeneFile="/sc/arion/projects/roussp01b/resources/databases/cage/linc_analyzed/lv3_robust.tsv.gz",
  ENSEMBL_INFO="/sc/arion/projects/roussp01b/resources/databases/gene-info/ensembl/my-downloads/muchEnsemblInfo_hg19.tsv.gz", #Ensembl genes as described in "STEP3.R"
  refPanel="/sc/arion/projects/roussp01b/resources/software/magma/magma_1_08/2019-02-25/g1000_eur.bim",
  dryRun=F
){

#hg19!!!!!
chromSizes=read.csv(stringsAsFactors=F, header=F,text="1,249250621
2,243199373
3,198022430
4,191154276
5,180915260
6,171115067
7,159138663
8,146364022
9,141213431
10,135534747
11,135006516
12,133851895
13,115169878
14,107349540
15,102531392
16,90354753
17,81195210
18,78077248
20,63025520
19,59128983
22,51304566
21,48129895
X,155270560
Y,59373566")
  
  colnames(chromSizes)=c("chr","end")
  
  
  #Derived vars
  outDir=paste0(baseDir,"/genes")
  dir.create(outDir, showWarnings=F, recursive=T)
  myLogFile=paste0(outDir,"/log.txt")
  write(unclass(Sys.time()),file=myLogFile)
  
  
  #Functions
  mySystem=function(...){
    #Compared to the standard R "system"-function, we create a dry-run option, and make "outDir" available in all calls
    myCommand=paste0("outDir=",outDir,";", ...)
    print(myCommand)
    write(myCommand,file=myLogFile,append=T)
    if(!dryRun)
      system(myCommand)
  }
  
  
  rmMHC=function(z){
    z[!  (  z$chr==6  &  ((z$start>25000000 & z$start<35000000)|(z$end>25000000 & z$end<35000000))  )  ,] #hg19
  }
  
  
  addPadding=function(myGenes){
    #Now for padded genes (we also ensure it does stay within chromosomal boundaries
    if(any(! myGenes$chr %in% chromSizes$chr)) stop("there's one or more genes on chromosomes for which the size is not given in chromSizes table.")
    if(any(!grepl("^(\\+|-)$",myGenes$strand))) stop("something's wrong with the strand information")
    myGenes$start = myGenes$start - ifelse(myGenes$strand=="+",35E3,10E3)
    myGenes$end = myGenes$end + ifelse(myGenes$strand=="+",10E3,35E3)
    myGenes$start = pmax(  myGenes$start,  1  ) 
    myGenes$end = pmin(  myGenes$end,  chromSizes[match(myGenes$chr,chromSizes$chr),"end"]  )
    if(any(myGenes$start>=myGenes$end)) stop("This is weird")
    return(myGenes)
  }
  
  
  writeAndAdd=function(x,prefix=""){
    xName=paste0(prefix,deparse(substitute(x)))
    xFile=paste0(outDir,"/",xName, ".tsv")
    write.table(x,file=xFile,sep="\t",quote=F,row.names=F,col.names=F)
    annotatedGeneFile=paste0("$outDir/", "snpToGene_", xName)
    mySystem(myMagma, " --annotate --snp-loc ", refPanel, " --gene-loc ", xFile, " --out ", annotatedGeneFile )
  }
  
  
  genStandardSet=function(GenesAll){
    myPrefix=deparse(substitute(GenesAll))
  
    writeAndAdd(GenesAll,myPrefix)
    
    GenesNoChryNoBmhc = rmMHC(  GenesAll[GenesAll$chr %in% c(1:22,"X"),]  )
    writeAndAdd(GenesNoChryNoBmhc,myPrefix)
    
    GenesAutosomesNoBmhc = rmMHC(  GenesAll[GenesAll$chr %in% 1:22,]  )
    writeAndAdd(GenesAutosomesNoBmhc,myPrefix)
    
    Genes35kb10kbAll = addPadding(GenesAll)
    writeAndAdd(Genes35kb10kbAll,myPrefix)
    
    Genes35kb10kbNoChryNoBmhc = rmMHC(  Genes35kb10kbAll[Genes35kb10kbAll$chr %in% c(1:22,"X"),]  )
    writeAndAdd(Genes35kb10kbNoChryNoBmhc,myPrefix)
    
    Genes35kb10kbAutosomesNoBmhc = rmMHC(  Genes35kb10kbAll[Genes35kb10kbAll$chr %in% 1:22,]  )
    writeAndAdd(Genes35kb10kbAutosomesNoBmhc,myPrefix)
  }
  
  
  #################################
  # Creating reference genes
  
  #standard genes (entrez)
  magma=read.table(standardGeneFile,col.names=c("geneID","chr","start","end","strand","geneName"),stringsAsFactors=F)
  genStandardSet(magma)

  #linc genes
  linc=read.delim(lincGeneFile,stringsAsFactors=F)
  linc$geneID=linc$ID
  linc$chr=sub("^chr","",linc$chr)
  linc=linc[,c("geneID","chr","start","end","strand","geneName")]
  genStandardSet(linc)

  #Ensembl genes
  library(plyr)
  ensemblInfo=read.delim(ENSEMBL_INFO,stringsAsFactors=F)
  ensemblProtCod=unique(ensemblInfo[ensemblInfo$Gene.type=="protein_coding", c("Ensembl.Gene.ID", "Strand", "Gene.End..bp.", "Gene.Start..bp.", "Associated.Gene.Name", "Chromosome.Name")])
  ensemblProtCod=ensemblProtCod[grepl("^(\\d|X|Y)",ensemblProtCod$Chromosome.Name,perl=T),]
  ensemblProtCod$Strand=ifelse(ensemblProtCod$Strand==1,"+","-")
  ensemblProtCod=unique(ensemblProtCod)

  ##NOTE: currently not addressing duplicated gene names, this should be an untested fix. NOTE: one could choose to filtering with ensemblInfo$Transcript.type=="protein_coding"
  #ensemblProtCod=ensemblProtCod[order(ensemblProtCod$ID),]
  #ensemblProtCodWithDuplicatedNames=unique(ensemblProtCod$Associated.Gene.Name[duplicated(ensemblProtCod$Associated.Gene.Name)])
  #ensemblProtCod$geneName=ensemblProtCod$Associated.Gene.Name
  #z=setNames(rep(0,length(ensemblProtCodWithDuplicatedNames)),ensemblProtCodWithDuplicatedNames)
  #for(i in 1:nrow(ensemblProtCod)){
  #  j=ensemblProtCod$geneName[i]
  #  if(j %in% ensemblProtCodWithDuplicatedNames){
  #    z[[j]]=z[[j]]+1
  #    ensemblProtCod$geneName[i]=paste0(j,"_v",z[[j]])
  #  }
  #}
  #ensemblProtCod$Associated.Gene.Name_orig=ensemblProtCod$Associated.Gene.Name
  #ensemblProtCod$Associated.Gene.Name=ensemblProtCod$geneName

  myCols=c(
    "Chromosome.Name"="chr",
    "Gene.Start..bp."="start",
    "Gene.End..bp."="end",
    "Strand"="strand",
    "Ensembl.Gene.ID"="geneID",
    "Associated.Gene.Name"="geneName"
  )
  colnames(ensemblProtCod)=revalue(colnames(ensemblProtCod),myCols)
  ensemblProtCod=ensemblProtCod[,c("geneID","chr","start","end","strand","geneName")]
  rm(myCols,ensemblInfo)
  genStandardSet(ensemblProtCod)
}
##################################################################


##################################################################
# calculating gene p-values
##################################################################
# NOTE: all dirs must be absolute

getGenePvals=function(
  baseDir="/sc/arion/projects/roussp01b/resources/databases/magma/v004",
  myMagma="/sc/arion/projects/roussp01b/resources/software/magma/magma_1_10/magma",
  refPanelPrefix="/sc/arion/projects/roussp01b/resources/software/magma/magma_1_08/2019-02-25/g1000_eur",
  snpAnnotatedGenes=c("ensemblProtCodGenes35kb10kbNoChryNoBmhc","lincGenes35kb10kbNoChryNoBmhc","magmaGenes35kb10kbNoChryNoBmhc","ensemblProtCodGenes35kb10kbAutosomesNoBmhc","lincGenes35kb10kbAutosomesNoBmhc","magmaGenes35kb10kbAutosomesNoBmhc"),
  mafCutOff=0.01,
  infoCutOff=0.8,
  runID="run1",
  dryRun=F,
  accountForJobs="acc_roussp01a",
  clusterSection="bode"
){
  outDir=baseDir

  settings=list(
    dryRun=dryRun,
    outDir=outDir,
    accountForJobs=accountForJobs,
    clusterSection=clusterSection,
    runID=runID
  )

  gwasDir=paste0(baseDir,"/filteredSumstats_maf_",mafCutOff,"_info_", infoCutOff)
  gwasPattern="_sumstats-RenamedAndFiltered.gz$"
  gwasFiles=list.files(gwasDir,pattern=gwasPattern,full.names=T)
  if(any(duplicated(sub(gwasPattern,"",basename(gwasFiles))))) stop("unfotunate choice of sumstat names")
  sumstatTempDir=paste0(baseDir,"/temp_unzipped_sumstats")

  dir.create(paste0(outDir,"/scripts"),recursive=T,showWarnings=F)
  dir.create(paste0(outDir,"/log"),recursive=T,showWarnings=F)
  dir.create(sumstatTempDir,recursive=T,showWarnings=F)

  analyses=expand.grid(gwasFiles=gwasFiles,snpAnnotatedGenes=snpAnnotatedGenes,stringsAsFactors=F)
  analyses$gwasNames=sub(gwasPattern,"",basename(analyses$gwasFiles))
  analyses$tempSumstats=paste0(sumstatTempDir,"/",analyses$gwasNames)
  analyses$tempSumstat_commands=paste0("zcat ", analyses$gwasFiles, " > ", analyses$tempSumstats)
  analyses$snpAnnotatedGenesFile=paste0(baseDir,"/genes/snpToGene_",analyses$snpAnnotatedGenes,".genes.annot")
  analyses$analysisNames=paste0(analyses$snpAnnotatedGenes,"__",analyses$gwasNames)
  analyses$outFile=paste0(baseDir,"/genePvals_maf_",mafCutOff,"_info_", infoCutOff,"/", analyses$analysisNames, "/", analyses$analysisNames)
  analyses$command=paste0("module purge;module load glibmm;", myMagma, " --bfile ", refPanelPrefix, " synonyms=0 --pval ", analyses$tempSumstats, " ncol=N --gene-annot ", analyses$snpAnnotatedGenesFile, " --out ",analyses$outFile) #we've already fixed the vars so no synonyms

  lapply(unique(analyses$tempSumstat_commands),function(x){message(x);system(x)})


  pvalScript=paste0("
    myCommand=$(sed \"${LSB_JOBINDEX}q;d\" <<< '", paste(analyses$command,collapse="\n"), "')
    myDir=$(sed \"${LSB_JOBINDEX}q;d\" <<< '", paste(dirname(analyses$outFile),collapse="\n"), "')

    mkdir -p $myDir
    cd $myDir
    eval $myCommand
  ")
  calcPvalJob=createScript(pvalScript,settings,mem=5E3,walltime="1:00",arrayLength=nrow(analyses))


  #lapply(paste("rm",unique(analyses$tempSumstats)),system)

  return(analyses)
}
##################################################################


##################################################################
# calc magma gene set pvals
##################################################################

calcMagmaGeneSetPvals=function(
  outDir,
  geneMetaSets, #either a named list, or in the same format as in downstream_greatR.R potentially containing multiple gene set sets
  pvalInfo, #data.frame must contain cols rawGenePvals, stdGenePvals, pvalSetName, gwasAcronym, gwasTrait. Further, file names must be absolute paths
  geneInfo=NA, #info about the genes e.g. geneName. NOT CURRENTLY IMPLEMENTED
  runID="run1",
  dryRun=F,
  accountForJobs="acc_roussp01a",
  clusterSection="bode",
  checkIfOutDirExists=T,
  myMagma="/sc/arion/projects/roussp01b/resources/software/magma/magma_1_10/magma",
  magma1_07_or_newer=T,
  background=NULL,
  covariatesTable=NULL,
  covariatesName=NULL,
  nperm=NULL #for fwer in magma 1.06 or prior (eg 10000). Must be NULL for 1.07 or newer
){

  message("initializing")
  
  if(checkIfOutDirExists & file.exists(outDir)) myStop("output dir mustn't already exist")
  createOutputFileStructure(outDir)
  outDir=fixRelativePath(outDir)
  myMagma=fixRelativePath(myMagma,check=T)

  settings=list(
    dryRun=dryRun,
    outDir=outDir,
    accountForJobs=accountForJobs,
    clusterSection=clusterSection,
    runID=runID
  )


  #if the input is simply a named list convert it to the format we used in greatR (a somewhat uggly check)
  if(   if(length(names(geneMetaSets[[1]]))!=2){ T }else{ !all(sort(names(geneMetaSets[[1]]))==c("metadata","sets")) }    ){
    message("converting geneSet format to contain metadata")
    geneMetaSets=list(myGeneSets=list(
      sets=setNames(geneMetaSets,make.names(names(geneMetaSets))),
      metadata=data.frame(name=make.names(names(geneMetaSets)),name_full=names(geneMetaSets),stringsAsFactors=F)
    ))
    rownames(geneMetaSets$myGeneSets$metadata)=geneMetaSets$myGeneSets$metadata$name
  }

  #check that gene set meta set names are safe, otherwise fix them
  if(any(duplicated(make.names(names(geneMetaSets))))) myStop("unfortunate choice of gene set group names")

  #check that all names of gene sets are safe:
  for(geneSetName in names(geneMetaSets)){
    z=make.names(names(geneMetaSets[[geneSetName]]$sets))
    if(any(duplicated(z))) myStop("unfortunate choice of gene set names")
    names(geneMetaSets[[geneSetName]]$sets)=z
    z=make.names(geneMetaSets[[geneSetName]]$metadata$name)
    if(any(duplicated(z))) myStop("unfortunate choice of gene set names")
    geneMetaSets[[geneSetName]]$metadata$name=z
  }
  rm(geneSetName)

  #test pvalInfo
  if(any(duplicated(pvalInfo$pvalSetName))) myStop("duplicated pvalSetName in pvalInfo")

  #grab all gene pvals and save it
  grabGenePvals=function(pvalSetName,myFile) data.frame(pvalSetName=pvalSetName,read.table(myFile,stringsAsFactors=F,header=T),stringsAsFactors=F)
  genePvals=do.call(rbind,mapply(grabGenePvals,pvalInfo$pvalSetName,pvalInfo$stdGenePvals,SIMPLIFY=F))
  mtsv(genePvals,outDir,gz=T)
  
  if(!is.null(background)) {
    write.table(background, quote=F, col.names=F, row.names=F, file=file.path(outDir, "meta-files", "background.lst"))
  }
  if(!is.null(covariatesTable) & !is.null(covariatesName)) {
    write.table(covariatesTable, quote=F, col.names=F, row.names=F, file=file.path(outDir, "meta-files", "covar.tble"))
  }
  
  #make table of tests to conduct
  testInfo=expand.grid(geneMetaSets=names(geneMetaSets),pvalSetName=pvalInfo$pvalSetName,stringsAsFactors=F)
  testInfo$testID=make.names(apply(testInfo,1,paste,collapse=" "))
  testInfo$outPrefix=paste0(outDir,"/results/geneSetPvals/",testInfo$testID,"/",testInfo$testID)
  testInfo$geneSetFile=paste0(outDir,"/files/geneMetaSets/",make.names(testInfo$geneMetaSets))
  testInfo$outFile=paste0(testInfo$outPrefix,ifelse(magma1_07_or_newer,".gsa.out",".sets.out"))
  testInfo=cMerge(testInfo,pvalInfo,"pvalSetName")
  permutationString=if(!is.null(nperm)){paste0(" --model fwer=", sprintf("%i",nperm))} else{""} #only works with magma 1.06 or prior
  geneIncludeLimit=if(!is.null(background)){paste0(" --settings gene-include=", file.path(outDir, "meta-files", "background.lst"))} else{""}
  covarInclude=if(!is.null(covariatesTable)){paste0(" --gene-covar ", file.path(outDir, "meta-files", "covar.tble"), " --model direction=pos condition='", covariatesName, "'")} else{""}
  
  testInfo$command=paste0(
    myMagma,
    " --gene-results ", testInfo$rawGenePvals, 
    " --set-annot ", testInfo$geneSetFile, " col=2,1",
    permutationString,
    geneIncludeLimit,
    covarInclude,
    " --out ", testInfo$outPrefix
  )
  
  rm(permutationString)

  #create output dirs
  sapply(dirname(testInfo$outPrefix),dir.create,showWarnings=F,recursive=T)
  sapply(unique(dirname(testInfo$geneSetFile)),dir.create,showWarnings=F,recursive=T)


  #save gene sets to files for magma
  writeGeneSetFile=function(geneSetName,myFile){
    df=do.call(rbind,lapply(names(geneMetaSets[[geneSetName]]$sets),function(setName)data.frame(gene.set=setName,gene=geneMetaSets[[geneSetName]]$sets[[setName]])))
    write.table(df,file=myFile,quote=F,sep="\t",col.names=F,row.names=F)
  }
  mapply(writeGeneSetFile,unique(testInfo$geneMetaSets),unique(testInfo$geneSetFile))

  #create script to calc gene setPvals
  genePvalScript=paste0("
    myCommand=$(sed \"${LSB_JOBINDEX}q;d\" <<< '", paste(testInfo$command,collapse="\n"), "')
    myDir=$(sed \"${LSB_JOBINDEX}q;d\" <<< '", paste(dirname(testInfo$outFile),collapse="\n"), "')

    mkdir -p $myDir
    cd $myDir
    eval $myCommand
  ")
  calcGenePvalJob=createScript(genePvalScript,settings,mem=5E3,walltime="1:00",arrayLength=nrow(testInfo))


  #FIXME: add aggegation script and copy sourced files to test dir
  gwasAggMagma=paste0( "plotAndAggMagma(\"", outDir, "\")" )
  gwasAggMagmaJob=createScript(gwasAggMagma,settings,walltime="1:00",isRscript=T,dependencies=list(calcGenePvalJob))
  rm(gwasAggMagma)


  mtsv(testInfo,outDir)

  save(list=ls(all=T), file=paste0(outDir, "/meta-files/calcMagmaGeneSetPvals.Rdata"), envir=environment())

  return(list(calcGenePvalJob=calcGenePvalJob,testInfo=testInfo,geneMetaSets=geneMetaSets))
}
##################################################################


##################################################################
# Aggregation function
##################################################################

plotAndAggMagma=function( #FIXME: not redone
  myInput, #currently needs to be output dir of previous step, but I've also thought about implementing it using just a data frame
  outDir=NA, #if NA and myInput is output dir of prev step, we output aggregated results to same dir
  calcFdr=T
  ){

  library(ggplot2)
  library(RColorBrewer)

  if(!is.data.frame(myInput)){
    inputWasDf=F
    load(paste0(myInput,"/meta-files/calcMagmaGeneSetPvals.Rdata"))
    if(is.na(outDir)) outDir=myInput
  }else{
    stop("not currently implemented")
    testInfo=myInfo #NOT CURENTLY IMPLEMENTED
    inputWasDf=T
  }

  #initialize
  if(!is.na(outDir))
    dir.create(paste0(outDir,"/results/"),recursive=T,showWarnings=F)

  #grab results
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
  geneSetResults=do.call(rbind,mapply(magmaGeneSetResultReader,testInfo$testID,testInfo$outFile,SIMPLIFY=F,USE.NAMES=F))

  #add most relvant cols from testInfo
  geneSetResults=cMerge(geneSetResults,testInfo[,colnames(testInfo) %in% c("geneMetaSets","geneTypeName","gwasAcronym","gwasTrait", "testID"),drop=F],"testID")
  
  #convert gwas acronyms to ordered factors for proper plotting
  geneSetResults$gwasAcronym=ordered(geneSetResults$gwasAcronym,levels=unique(geneSetResults$gwasAcronym))
  geneSetResults$gwasTrait=ordered(geneSetResults$gwasTrait,levels=unique(geneSetResults$gwasTrait))


  #add gene names #FIXME: make check in earlier step that all metadata contains name and name_full cols
  if(exists("geneMetaSets") & "geneMetaSets" %in% colnames(geneSetResults)){
    z=do.call(rbind,lapply(names(geneMetaSets),function(x)data.frame(geneMetaSets=x,geneMetaSets[[x]]$metadata[,c("name","name_full")],stringsAsFactors=F)))
    z$metaName=paste(z$name,z$geneMetaSets)
    geneSetResults$name_full=z$name_full[match(paste(geneSetResults$name,geneSetResults$geneMetaSets),z$metaName)]
  }else{
    geneSetResults$name_full=geneSetResults$name
  }

  #Sort by pval
  geneSetResults$LogP=-log10(geneSetResults$P)

  #save that as a table
  mtsv(geneSetResults,outDir,gz=T)

  #IDEA: we could also in addtion save the results by themselves

  #PLOTS:  always per geneTypeName. both pval and BETA
  # barplots all gs (if not too many)
  #  of top gs 
  # barplots all gs (if not too many)
  # barplots of top gs 
  # heatmaps all gs (if not too many)
  # heatmaps of top gs 
  # barplot fnct should be able to also plot one gene set in all traits
  # w/wo clustering

  #################################
  #heatmaps
  

  #################################
  #FDR functions:
  doFdr=function(x){
    if(calcFdr){
      #library(fdrtool)
      #x$global_FDR_AdjP=fdrtool(x$P, statistic="pvalue", plot=F,cutoff.method="fndr",verbose=F)$qval
      x$global_Bonf_AdjP=p.adjust(x$P, method="bonferroni")
      x$global_BH_AdjP=p.adjust(x$P, method="BH")
      x$plotLabel=""
      x$plotLabel[x$P<0.05]="·"
      x$plotLabel[x$global_BH_AdjP<0.05]="#"
    }else{
      x$plotLabel="-"
    }
    return(x)
  }

  #################################
  # aggregate by geneMetaSets (and geneTypeName)

  geneSetGroupAggOverview=unique(geneSetResults[,c("geneTypeName","geneMetaSets")])
  geneSetGroupAggOverview$aggID=make.names(apply(geneSetGroupAggOverview,1,paste,collapse=" "))
  if(any(duplicated(geneSetGroupAggOverview$aggID))) stop("unfortunate combination of geneTypeName and geneMetaSets led to duplicated aggID")
  rownames(geneSetGroupAggOverview)=geneSetGroupAggOverview$aggID

  doGeneSetGroupAgg=function(aggID,geneTypeName,geneMetaSets){
    z=geneSetResults[geneSetResults$geneTypeName==geneTypeName & geneSetResults$geneMetaSets==geneMetaSets,]
    z$aggID=aggID
    z=doFdr(z)
    return(z)
  }
  geneSetGroupAgg=mapply(doGeneSetGroupAgg,geneSetGroupAggOverview$aggID,geneSetGroupAggOverview$geneTypeName,geneSetGroupAggOverview$geneMetaSets,SIMPLIFY=F)


  #################################
  # Heatmap plots

  if(!is.na(outDir)){
    message("GSEA: plotting")

    #plot all pathways
    dir.create(outDir,showWarnings=F,recursive=T)
    lapply(names(geneSetGroupAgg),function(i){
      if(length(unique(geneSetGroupAgg[[i]]$name))<101){ #we cannot plot like 4k pathways
        mySetName=geneSetGroupAgg[[i]]$geneMetaSets[1]
        subDir=paste0(outDir,"/results/plots/",mySetName)
        dir.create(subDir,showWarnings=F,recursive=T)
        aggMagmaPlotter(df=geneSetGroupAgg[[i]],annoName=i,outDir=subDir)
        aggMagmaPlotter(df=geneSetGroupAgg[[i]],annoName=i,outDir=subDir,doCluster=T,geneSets=geneMetaSets[[mySetName]]$sets)
      }
    })
      
    #Top something best pathways from each peakSet
    topGeneSetGroupAgg=list()
    for(i in names(geneSetGroupAgg)){
      w=geneSetGroupAgg[[i]]
      w=w[order(w$P),]
      for(myTop in c(5,10,15,20,25)){
        myName=paste0(i,"_top",myTop)
        if(length(unique(w$name))>=myTop){
          myPathways=unique(unlist(lapply(unique(w$gwasTrait),function(x){w$name_full[w$gwasTrait==x][1:myTop]})))
          topGeneSetGroupAgg[[myName]]=w[w$name_full %in% myPathways,]
          mySetName=topGeneSetGroupAgg[[myName]]$geneMetaSets[1]
          subDir=paste0(outDir,"/results/plots/",mySetName)
          dir.create(subDir,showWarnings=F,recursive=T)
          aggMagmaPlotter(df=topGeneSetGroupAgg[[myName]], annoName=myName, outDir=subDir)
          aggMagmaPlotter(df=topGeneSetGroupAgg[[myName]], annoName=myName, outDir=subDir, doCluster=T,  geneSets=geneMetaSets[[mySetName]]$sets)
          #if(length(unique(topGeneSetGroupAgg[[myName]]$gwasTrait))>1) #we can only bicluster with 2+ samples
          #  plotGseaBiclust(df=topGeneSetGroupAgg[[myName]], annoName=myName, outDir=subDir)
        }
      }
    }
    rm(i,myTop,w,myPathways,myName)
  }

}
##################################################################


##################################################################
# Plotting function
##################################################################
  
aggMagmaPlotter=function(df,annoName,outDir,doCluster=F,geneSets=NA){
  message(paste0("plotting ", annoName))
  
  #init
  library(ggplot2)
  library(RColorBrewer)
  myPalette = colorRampPalette(brewer.pal(9, "Greens"), space="Lab")
  myPalette2way=colorRampPalette(rev(c("#67001F","#B2182B","#D6604D","#F4A582","#FDDBC7","#F7F7F7","#D1E5F0","#92C5DE","#4393C3","#2166AC","#053061")), space="Lab")


  if(doCluster){
    if(!is.list(geneSets)) stop("if you wan't to cluster, you must provide the original gene sets as a list")
    df$name_full=as.character(df$name_full)
    myClustering=geneSetClusteringMagma(metadata=df,sets=geneSets,plotName=annoName,outDir=outDir)
    df$name_full=ordered(df$name_full,levels=df$name_full[match(myClustering$orderedGeneSets_Ward.D2,df$name)])
  }

  #plots

  #proper symmetric scaling for BETA
  if(!all(is.infinite(df$BETA))){
    myLims=max(abs(df$BETA[!is.infinite(df$BETA)]))
    if(myLims==0) myLims=1
    myLims=c(-1,1)*myLims
  }else{ myLims=c(-1,1) }
  
  for(flipped in c(F,T)){
    if(flipped){
      df$xCol=df$gwasTrait
      df$yCol=df$name_full
    }else{
      df$xCol=df$name_full
      df$yCol=df$gwasTrait
    }

    zz=ggplot(df,aes(xCol,yCol)) +
        ggtitle(annoName) +
        scale_y_discrete(expand = c(0, 0)) +
        scale_x_discrete(expand = c(0, 0)) +
        ggtitle(annoName) +
        theme_classic() +
        coord_fixed() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(axis.title = element_blank())
  
    filename=paste0(outDir,"/plot_gsea_",annoName,"_clustered_",doCluster,ifelse(flipped,"_flipped.pdf",".pdf"))
    #if(flipped) zz= zz  + coord_flip() + scale_x_discrete(limits = rev(levels(df$name_full)), expand = c(0, 0))
    #how much space the names take up is rather unpredictable without extensive coding, so this will have to do. Also coord_fixed breaks with coord_flip, which I have not addressed.
    geneSetWidth=7+0.2*length(unique(df$name_full)) 
    peakWidth=7+0.2*length(unique(df$peaksFullName))
    myHeight=ifelse(flipped,geneSetWidth+0,peakWidth)
    myWidth=ifelse(flipped,peakWidth-0,geneSetWidth)

    pdf(filename,height=myHeight,width=myWidth);
      print( zz + geom_tile(aes(fill=LogP)) + scale_fill_gradientn(colours = myPalette(100)))
      print( zz + geom_tile(aes(fill=LogP)) + scale_fill_gradientn(colours = myPalette(100))+theme(legend.position="none"))
      print( zz + geom_tile(aes(fill=LogP)) + scale_fill_gradientn(colours = myPalette(100))+theme(legend.position="none")+geom_text(aes(label=plotLabel),alpha=0.7))
      print( zz + geom_tile(aes(fill=LogP)) + scale_fill_gradientn(colours = myPalette(100),trans="sqrt"))
      print( zz + geom_tile(aes(fill=LogP)) + scale_fill_gradientn(colours = myPalette(100),trans="sqrt")+theme(legend.position="none"))
      print( zz + geom_tile(aes(fill=LogP)) + scale_fill_gradientn(colours = myPalette(100),trans="sqrt")+theme(legend.position="none")+geom_text(aes(label=plotLabel),alpha=0.7))
      print( zz + geom_tile(aes(fill=BETA)) + scale_fill_gradientn(colours = myPalette2way(100),limits=myLims))
      print( zz + geom_tile(aes(fill=BETA)) + scale_fill_gradientn(colours = myPalette2way(100),limits=myLims)+theme(legend.position="none"))
      print( zz + geom_tile(aes(fill=BETA)) + scale_fill_gradientn(colours = myPalette2way(100),limits=myLims)+theme(legend.position="none")+geom_text(aes(label=plotLabel),alpha=0.7))
    dev.off()
  }
}
##################################################################


##################################################################
# Cluster gene sets based on Jaccard
##################################################################

geneSetClusteringMagma=function(metadata,sets,plotName=NA,outDir=NA){

  library(vegan)

  #lets restrict analysis to sets with metadata
  sets=sets[names(sets) %in% metadata$name]

  #create a matrix with gene membership
  myGenes=unique(unlist(sets))
  gm=matrix(data=F,nrow=length(myGenes),ncol=length(sets),dimnames=list(myGenes=myGenes,sets=names(sets)))
  for(i in names(sets)) gm[,i]=myGenes %in% sets[[i]]; rm(i)

  dist.mat=vegdist(data.frame(t(gm)),method="jaccard",binary=T)
  myClust_Complete=hclust(dist.mat)
  myClust_Ward.D2=hclust(dist.mat,method="ward.D2")

  #Plot if output dir is provided
  if(!is.na(outDir)){
    pdf(paste0(outDir,"/plot_gsea_",plotName,"_clustering.pdf"),height=18,width=5+0.2*length(sets));
      print(plot(myClust_Complete))
      print(plot(myClust_Complete,hang=-1))
      print(plot(myClust_Ward.D2))
      print(plot(myClust_Ward.D2,hang=-1))
    dev.off()
  }

  #version without labs to put together with heatmap in post
  if(!is.na(outDir)){
    pdf(paste0(outDir,"/plot_gsea_",plotName,"_clustering_no_lab.pdf"),height=3.5+0.2*length(sets)/4,width=3+0.2*length(sets));
      print(plot(myClust_Complete,labels = F))
      print(plot(myClust_Complete,labels = F,hang=-1))
      print(plot(myClust_Ward.D2,labels = F))
      print(plot(myClust_Ward.D2,labels = F,hang=-1))
    dev.off()
  }
  
  z=list()
  z$orderedGeneSets_Complete=myClust_Complete$labels[myClust_Complete$order]
  z$clustering_Complete=myClust_Complete
  z$orderedGeneSets_Ward.D2=myClust_Ward.D2$labels[myClust_Ward.D2$order]
  z$clustering_Ward.D2=myClust_Ward.D2
  return(z)
}
##################################################################


##################################################################
# Misc small helper functions that you can call yourself
##################################################################

#################################
# Grab some sample sumstats

grabGenePvalInfo=function(
  selectedGeneTypeName, #e.g. selectedGeneTypeName="ensemblProtCodGenes35kb10kbNoChryNoBmhc" # gwas pval files must be named in manner [geneTypeName]__[gwasAcronym] e.g. ensemblProtCodGenes35kb10kbNoChryNoBmhc__alz
  myPath="/sc/arion/projects/roussp01b/resources/databases/magma/v004/genePvals_maf_0.01_info_0.8/",
  sumstatMetadataFile="/sc/arion/projects/roussp01b/resources/databases/ldsc/metadata.tsv", #needs to have cols "gwasAcronym" and "gwasTrait"
  useOnlyTraitsFromMetadata=F
){
  #grab gene pvals
  myInfo=data.frame(rawGenePvals=list.files(myPath,full.names=T,recursive=T,pattern="\\.genes\\.raw$"),stringsAsFactors=F)
  myInfo$stdGenePvals=sub("\\.raw$",".out",myInfo$rawGenePvals)
  myInfo$pvalSetName=sub("\\.genes\\.raw$","",basename(myInfo$rawGenePvals))
  if(any(unlist(lapply(strsplit(myInfo$pvalSetName,split="__"),length))!=2)) myStop("unfortunate choice of naming causes error in intifying traits")
  myInfo$geneTypeName=do.call(rbind,strsplit(myInfo$pvalSetName,split="__"))[,1]
  myInfo$gwasAcronym=do.call(rbind,strsplit(myInfo$pvalSetName,split="__"))[,2]
  myInfo=myInfo[myInfo$geneTypeName %in% selectedGeneTypeName,]
  if(any(!file.exists(myInfo$stdGenePvals))) myStop("one or more genes.out file(s) are missing")

  #add gwas metadata
  ssmd=read.delim(sumstatMetadataFile,stringsAsFactors=F, quote="")
  ssmd=ssmd[ssmd$gwasAcronym %in% myInfo$gwasAcronym,]
  w=c(ssmd$gwasAcronym,myInfo$gwasAcronym[!myInfo$gwasAcronym %in% ssmd$gwasAcronym])
  myInfo=cbind(   myInfo[match(w,myInfo$gwasAcronym),],   ssmd[match(w,ssmd$gwasAcronym),colnames(ssmd)!="gwasAcronym"]   )
  myInfo$myName=myInfo$gwasTrait
  myInfo$myName[is.na(myInfo$myName)]=myInfo$gwasAcronym[is.na(myInfo$myName)] #if we don't have metadata just use the acronym as a name
  myInfo$gwasTrait=myInfo$myName
  myInfo$myName=NULL
  
  if(useOnlyTraitsFromMetadata) {
    myInfo=myInfo[myInfo$gwasAcronym %in% ssmd$gwasAcronym,]
  }
  
  return(myInfo)
}

#################################
# Test set-up with classic gene sets
# NOT YET DONE
