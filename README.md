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

## Contact
Dr. Liting Song (liting.song@mssm.edu)


