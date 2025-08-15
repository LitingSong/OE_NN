# Single-Nucleus Analysis of the Adult Human Olfactory Epithelium Uncovers Shared Neurogenesis Programs with the Brain

## Figure 1: The cell taxonomy of olfactory epithelium

s1.1 merge MSSM-OE data: 

&emsp;*01a_merge_outer_entdata.ipynb*

s1.2 QC and data processing: 

&emsp;*01b.ent_process_outerMerge_0807.ipynb & 01c.NN_process-0807.ipynb*

s1.3 integrate OSNs from MSSM-OE dataset and Durante et al’s data 

&emsp;*01c.integrate_Neuron_nn_ent-0807.ipynb*

s1.4 Cell stage composition and nasal site origin for OSNs; statistics of cell cycle score: 

&emsp;*01d.statstitic.ipynb*

s1.5 scDRS association between brain-related disease and OSNs 

&emsp;*01e.scdrs_downstream_stage.R*


## Figure 2: Dynamics across OSN genesis and development. 

s2 trajectory inference & identify trajDEGs & identify clusters & functional enrichment

&emsp;*02a.traj_integrate_600-6-40-50-graph_k7_FDR_BG-Cov.ipynb*

## Figure 3: Similarities and differences in gene expression between OSNs and CENs during development

s3.1 expression correlation between OSNs and brain cells 

&emsp;*03a_OSN_vs_BRAIN_lister-corr.ipynb & 03a_OSN_vs_BRAIN_liwang-corr.ipynb*

s3.2 gene module score

&emsp;*03b_module_score.ipynb & 03b_module_score_heatmap.ipynb*

s3.3 overlap between trajDEGs between OSNs and CENs

&emsp;*03c.OSN_CEN_DEoverlap.ipynb*

s3.4 scatter plot of enriched GOBPs

&emsp;*03d.trend_pattern_OlfvsBrain-gprofiler2-ExN2cux2-k7_FDR_BG-graph.ipynb*

## Figure 4: Expression trajectory alignments between OSNs and CENs for ASD risk genes

s4.1 Expression trajectory alignments and identification of alignment patterns

&emsp;*04a_genes2genes-psyad-sfari.ipynb*

s4.2 Functional enrichment of genes in different alignment patterns

&emsp;*04b_genes2genes_pathway-psyad-sfari.ipynb*

## Fig. 5: Association of OSN stage-specific TFs with neuropsychiatric disorders
s5.1.1 GRN construction by SCENIC

&emsp;*05a.ENTNN_N_SCENIC-allG.py*

s5.1.2 MAGMA enrichment

&emsp;*05a_downstream_magma.R & 05a_run_magma.R*


S5.2 Visualization of TF regulon activity and enrichment in psychiatric disorders.

&emsp;*05b.TF_heatmap_tcf4_network.ipynb*

## R sessionInfo()
R version 4.3.1 (2023-06-16)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: Rocky Linux 9.4 (Blue Onyx)

Matrix products: default
BLAS/LAPACK: /sc/arion/projects/roussp01a/liting/software/conda/Monocle3/lib/libopenblasp-r0.3.21.so;  LAPACK version 3.9.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: America/New_York
tzcode source: system (glibc)

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] stringr_1.5.1               uwot_0.1.16                
 [3] Matrix_1.6-5                reticulate_1.40.0          
 [5] scuttle_1.12.0              BiocManager_1.30.25        
 [7] patchwork_1.3.0             zellkonverter_1.13.3       
 [9] dplyr_1.1.4                 circlize_0.4.16            
[11] ComplexHeatmap_2.18.0       slingshot_2.10.0           
[13] TrajectoryUtils_1.10.1      SingleCellExperiment_1.24.0
[15] SummarizedExperiment_1.32.0 Biobase_2.62.0             
[17] GenomicRanges_1.54.1        GenomeInfoDb_1.38.8        
[19] IRanges_2.36.0              S4Vectors_0.40.2           
[21] BiocGenerics_0.48.1         MatrixGenerics_1.14.0      
[23] matrixStats_1.3.0           princurve_2.1.6            
[25] Seurat_5.0.3                SeuratObject_5.0.2         
[27] sp_2.1-3                    ggVennDiagram_1.5.2        
[29] SeuratWrappers_0.3.5        schard_0.0.1               
[31] scales_1.4.0                RColorBrewer_1.1-3         
[33] scattermore_1.2             ggplot2_3.5.2              

loaded via a namespace (and not attached):
  [1] RcppAnnoy_0.0.22          splines_4.3.1            
  [3] later_1.3.1               pbdZMQ_0.3-11            
  [5] filelock_1.0.3            bitops_1.0-7             
  [7] tibble_3.2.1              R.oo_1.26.0              
  [9] polyclip_1.10-6           basilisk.utils_1.14.1    
 [11] fastDummies_1.7.4         lifecycle_1.0.4          
 [13] doParallel_1.0.17         globals_0.18.0           
 [15] lattice_0.22-6            MASS_7.3-60              
 [17] magrittr_2.0.3            plotly_4.10.4            
 [19] remotes_2.4.2.1           httpuv_1.6.11            
 [21] sctransform_0.4.1         spam_2.10-0              
 [23] spatstat.sparse_3.1-0     cowplot_1.1.3            
 [25] pbapply_1.7-2             abind_1.4-8              
 [27] zlibbioc_1.48.0           Rtsne_0.17               
 [29] purrr_1.0.2               R.utils_2.12.3           
 [31] RCurl_1.98-1.14           GenomeInfoDbData_1.2.11  
 [33] ggrepel_0.9.5             irlba_2.3.5.1            
 [35] listenv_0.9.1             spatstat.utils_3.1-0     
 [37] goftest_1.2-3             RSpectra_0.16-1          
 [39] spatstat.random_3.3-2     fitdistrplus_1.2-1       
 [41] parallelly_1.37.1         DelayedMatrixStats_1.24.0
 [43] leiden_0.4.3.1            codetools_0.2-19         
 [45] DelayedArray_0.28.0       shape_1.4.6.1            
 [47] tidyselect_1.2.1          farver_2.1.1             
 [49] base64enc_0.1-3           spatstat.explore_3.3-2   
 [51] jsonlite_1.8.7            GetoptLong_1.0.5         
 [53] progressr_0.14.0          ggridges_0.5.6           
 [55] survival_3.5-8            iterators_1.0.14         
 [57] foreach_1.5.2             tools_4.3.1              
 [59] ica_1.0-3                 Rcpp_1.0.11.6            
 [61] glue_1.6.2                gridExtra_2.3            
 [63] SparseArray_1.2.2         IRdisplay_1.1            
 [65] withr_3.0.2               fastmap_1.1.1            
 [67] basilisk_1.14.3           digest_0.6.33            
 [69] rsvd_1.0.5                R6_2.6.1                 
 [71] mime_0.12                 colorspace_2.1-0         
 [73] tensor_1.5                dichromat_2.0-0.1        
 [75] spatstat.data_3.1-2       R.methodsS3_1.8.2        
 [77] tidyr_1.3.1               generics_0.1.3           
 [79] data.table_1.15.2         httr_1.4.7               
 [81] htmlwidgets_1.6.2         S4Arrays_1.2.0           
 [83] pkgconfig_2.0.3           gtable_0.3.6             
 [85] lmtest_0.9-40             XVector_0.42.0           
 [87] htmltools_0.5.6.1         dotCall64_1.1-1          
 [89] clue_0.3-65               png_0.1-8                
 [91] spatstat.univar_3.0-1     rjson_0.2.21             
 [93] reshape2_1.4.4            uuid_1.2-0               
 [95] nlme_3.1-164              GlobalOptions_0.1.2      
 [97] repr_1.1.7                zoo_1.8-12               
 [99] KernSmooth_2.23-22        parallel_4.3.1           
[101] miniUI_0.1.1.1            pillar_1.10.2            
[103] vctrs_0.6.4               RANN_2.6.1               
[105] promises_1.2.1            beachmat_2.18.0          
[107] xtable_1.8-4              cluster_2.1.6            
[109] evaluate_1.0.3            cli_3.6.1                
[111] compiler_4.3.1            rlang_1.1.1              
[113] crayon_1.5.3              future.apply_1.11.2      
[115] plyr_1.8.9                stringi_1.7.12           
[117] BiocParallel_1.36.0       viridisLite_0.4.2        
[119] deldir_2.0-4              lazyeval_0.2.2           
[121] spatstat.geom_3.3-3       dir.expiry_1.10.0        
[123] IRkernel_1.3.2            RcppHNSW_0.6.0           
[125] sparseMatrixStats_1.14.0  future_1.33.2            
[127] shiny_1.9.1               ROCR_1.0-11              
[129] igraph_2.0.3       

## python version
3.9.19 | packaged by conda-forge | (main, Mar 20 2024, 12:50:21) [GCC 12.3.0]

## python libraries
typing_extensions==4.12.2
traitlets==5.14.3
jupyter_core==5.7.2
tornado==6.4.1
six==1.16.0
platformdirs==4.3.6
importlib_metadata==8.5.0
jupyter_client==8.6.3
ipykernel==6.29.5
executing==2.1.0
asttokens==2.4.1
pure_eval==0.2.3
stack_data==0.6.3
pygments==2.18.0
importlib.metadata==8.5.0
ptyprocess==0.7.0
pexpect==4.9.0
pickleshare==0.7.5
backcall==0.2.0
decorator==5.1.1
wcwidth==0.2.13
prompt_toolkit==3.0.48
parso==0.8.4
jedi==0.19.1
IPython==8.12.3
comm==0.2.2
psutil==6.1.0
packaging==24.1
debugpy==1.8.7
zipp==3.20.2
numpy==1.22.4
h5py==3.11.0
natsort==8.4.0
pytz==2024.2
pandas==2.0.3
scipy==1.10.1
tqdm==4.66.5
torch==2.4.1
anndata==0.9.2
pyparsing==3.1.4
cycler==0.12.1
kiwisolver==1.4.7
matplotlib==3.7.5
importlib_resources==6.4.5
statistics==1.0.3.5
patsy==0.5.6
statsmodels==0.14.1
seaborn==0.13.2
llvmlite==0.41.1
importlib.resources==6.4.5
numba==0.58.1
joblib==1.4.2
threadpoolctl==3.5.0
scanpy==1.9.8
regex==2024.9.11
tabulate==0.9.0
leven==1.0.4
urllib3==2.2.3
charset_normalizer==3.4.0
idna==3.10
certifi==2024.8.30
requests==2.32.3
gseapy==1.1.3
mpmath==1.3.0
blitzgsea==1.3.47
genes2genes==0.2.0

## Contact
Dr. Liting Song (liting.song@mssm.edu)


