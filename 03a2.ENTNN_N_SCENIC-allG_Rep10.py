
# create basic row and column attributes for the loom file:
import loompy as lp
import os
import numpy as np
import pandas as pd
import scanpy as sc
import glob
from cytoolz import compose
import numpy as np
from pyscenic.transform import df2regulons
np.bool = np.bool_
import os, glob, re, pickle
from pyscenic.aucell import aucell
import matplotlib.pyplot as plt
import seaborn as sns
from pyscenic.utils import modules_from_adjacencies, load_motifs
import subprocess
from pyscenic.binarization import binarize

sc.settings.verbosity = 0 # verbosity: errors (0), warnings (1), info (2), hints (3)
#sc.logging.print_versions()
sc.set_figure_params(dpi=150, fontsize=10, dpi_save=600)


# ranking databases
f_db_glob = "/sc/arion/projects/roussp01a/pengfei/tools/SCENIC/gene_based/*feather"
f_db_names = ' '.join( glob.glob(f_db_glob) )

# motif databases
f_motif_path = "/sc/arion/projects/roussp01a/pengfei/tools/SCENIC/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
f_tfs = "/sc/arion/projects/roussp01a/pengfei/tools/SCENIC/allTFs_hg38.txt" # human

cmc_dir='/sc/arion/projects/CommonMind/roussp01a/ENT/snRNAseq/'
ADJACENCIES_FNAME=os.path.join(cmc_dir, 'qc_scanpy/{}_allG.adj.csv'.format('ENTNN'))
REGULONS_FNAME=os.path.join(cmc_dir, 'qc_scanpy/{}_allG.reg.csv'.format('ENTNN'))
REGULONS_DAT_FNAME = os.path.join(cmc_dir, 'qc_scanpy/{}_allG.regulons.dat'.format('ENTNN'))
AUCELL_MTX_FNAME = os.path.join(cmc_dir, 'qc_scanpy/{}_allG.auc.csv'.format('ENTNN'))


N_adata_f = sc.read(cmc_dir+'qc_scanpy/ent_nn_merge_rawcount.h5ad')

meta_data = pd.read_csv('/sc/arion/projects/roussp01a/liting/Olf/data/ent_nn_merge_cca.metadata.csv', index_col=0)

N_adata_f = N_adata_f[meta_data.index]
N_adata_f.obs["N_types_stage"]= meta_data.loc[N_adata_f.obs.index,'cca_N_types_stage'] 



sc.pp.filter_genes(N_adata_f, min_cells=3)

df_exp=pd.DataFrame(N_adata_f.X.toarray(), index = N_adata_f.obs_names, columns = N_adata_f.var_names)

row_attrs = {
    "Gene": np.array(N_adata_f.var_names) ,
}
col_attrs = {
    "CellID": np.array(N_adata_f.obs_names) ,
    "nGene": np.array( np.sum(N_adata_f.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(N_adata_f.X.transpose() , axis=0)).flatten(),
}
lp.create( cmc_dir+'ENTNN_n_filtered_scenic_allG.loom', N_adata_f.X.transpose(), row_attrs, col_attrs)
f_loom_path_scenic=cmc_dir+'ENTNN_n_filtered_scenic_allG.loom'


def derive_regulons(motifs,):
    motifs.columns = motifs.columns.droplevel(0)
    
    def contains(*elems):
        def f(context):
            return any(elem in context for elem in elems)
        return f

    # For the creation of regulons we only keep the 10-species databases and the activating modules. We also remove the
    # enriched motifs for the modules that were created using the method 'weight>50.0%' (because these modules are not part
    # of the default settings of modules_from_adjacencies anymore.
    motifs = motifs[np.fromiter(map(contains('activating'), motifs.Context), dtype=np.bool)]

    # We build regulons only using enriched motifs with a NES of 3.0 or higher; we take only directly annotated TFs or TF annotated
    # for an orthologous gene into account; and we only keep regulons with at least 10 genes.
    regulons = list(filter(lambda r: len(r) >= 5, df2regulons(motifs[(motifs['NES'] >= 0.1) 
                                                                      & ((motifs['Annotation'] == 'gene is directly annotated')
                                                                        | (motifs['Annotation'].str.startswith('gene is orthologous to')
                                                                           & motifs['Annotation'].str.endswith('which is directly annotated for motif')))
                                                                     ])))
    
    # Rename regulons, i.e. remove suffix.
    return list(map(lambda r: r.rename(r.transcription_factor), regulons))



for repi in range(52,53):
    ## STEP 1: Network inference based on GRNBoost2 from CLI
    #!/hpc/users/songl05/.conda/envs/Py39_R43_Ju10/bin/pyscenic grn {f_loom_path_scenic} {f_tfs} -o {ADJACENCIES_FNAME} --num_workers 8
    command1=f"/hpc/users/songl05/.conda/envs/Py39_R43_Ju10/bin/pyscenic grn {f_loom_path_scenic} {f_tfs} -o {ADJACENCIES_FNAME} --num_workers 12"
    subprocess.run(command1, shell=True, capture_output=True, text=True)

    ## STEP 2-3: Regulon prediction aka cisTarget from CLI
    
    command2=f"/hpc/users/songl05/.conda/envs/Py39_R43_Ju10/bin/pyscenic ctx {ADJACENCIES_FNAME} \
        {f_db_names} \
        --annotations_fname {f_motif_path} \
        --expression_mtx_fname {f_loom_path_scenic} \
        --output {REGULONS_FNAME} \
        --mask_dropouts \
        --num_workers 12 --nes_threshold 2  "
    
    subprocess.run(command2, shell=True, capture_output=True, text=True)
    from pyscenic.utils import modules_from_adjacencies, load_motifs

    ### STEP 4: Cellular enrichment aka AUCell

    df_motifs = load_motifs(REGULONS_FNAME)
    regulons = derive_regulons(df_motifs)

    import warnings
    warnings.filterwarnings("ignore")
    dfs = []

    for i in range(df_motifs.shape[0]):
        df=pd.DataFrame(df_motifs["TargetGenes"][i],)
        df['tf']=df_motifs["TargetGenes"].index[i][0]
        df['MotifID']=df_motifs["TargetGenes"].index[i][1]
        dfs.append(df)

    tf_target = pd.concat(dfs, ignore_index=True)
    tf_target = pd.DataFrame(tf_target)
    tf_target = tf_target.rename(columns={'0': 'target','1':'weight'})


    #tf_target.to_csv(cmc_dir+'entnn_allG_tf_target.csv')
    tf_target.to_csv('/sc/arion/projects/CommonMind/liting/ENT/scenic/entnn_allG_tf_target'+ str(repi) + '.csv')


    #### AUCell

    auc_mtx = aucell(df_exp, regulons, num_workers=26) 
    #auc_mtx.to_csv(AUCELL_MTX_FNAME) 

    ### CELL TYPE SPECIFIC REGULATORS - Z-SCORE


    signature_column_names = list(auc_mtx.columns)

    import pandas as pd
 

    df_scores=pd.concat([auc_mtx, N_adata_f.obs['N_types_stage']],sort=False, axis=1, join='outer' )

    df_results = ((df_scores.groupby(by='N_types_stage').mean() - df_scores[signature_column_names].mean())/ df_scores[signature_column_names].std())#.stack().reset_index().rename(columns={'level_1': 'regulon', 0:'Z'})

    df_results.to_csv('/sc/arion/projects/CommonMind/liting/ENT/scenic/tf_aucscore_rep'+ str(repi) + '.csv')
    
    bin_mtx, thresholds = binarize(auc_mtx)
    bin_mtx.to_csv('/sc/arion/projects/CommonMind/liting/ENT/scenic/entnn_allG_tf_target'+ str(repi) + '.csv')
    thresholds.to_frame().rename(columns={0:'threshold'}).to_csv('/sc/arion/projects/CommonMind/liting/ENT/scenic/threshold'+ str(repi) + '.csv')

