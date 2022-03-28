#Import libraies needed 
import os
import scib
import pandas as pd 
import anndata
import scanpy as sc
#library(anndata, warn.conflicts = FALSE, quietly = TRUE)

#Change working directory
os.getcwd()
os.chdir('/datos_2/phD_raw_data_2/OUTPUTS_2_ae_peaks_test_7_100pcent/q0/Summer_phase_1_data/starter_kit-joint_embedding-python/')

#Metrics for integration by:: https://github.com/theislab/scib 
data = pd.read_csv("shared_space_01.csv") 
adata_latent_space = sc.read("openproblems_bmmc_multiome_phase1v2.method.output.h5ad")


#pred
from sklearn.metrics import mean_squared_error

os.getcwd()
os.chdir('/datos_2/phD_raw_data_2/OUTPUTS_2_ae_peaks_test_7_100pcent/q0/Summer_phase_1_data/starter_kit-predict_modality-python/')
pred = sc.read("openproblems_bmmc_multiome_phase1v2_mod2.method.output.h5ad")
pred = pred.X.toarray()

os.chdir('/datos_2/phD_raw_data_2/OUTPUTS_2_ae_peaks_test_7_100pcent/q0/Summer_phase_1_data/starter_kit-predict_modality-python/output/datasets_phase1v2/predict_modality/openproblems_bmmc_multiome_phase1v2_mod2/')
test = sc.read("openproblems_bmmc_multiome_phase1v2_mod2.censor_dataset.output_test_mod2.h5ad")
test = test.X.toarray()

rms = mean_squared_error(test, pred, squared=False)


#################
#annData for shared space
adata_latent_space = anndata.AnnData(X=data.iloc[0:,1:],obs= data.iloc[0:,0:1],var=data.iloc[0:0,1:].T, uns = {'dataset_id': 'organism', 'AE': 'AE', },)

#adata_latent_space.obs["cell"] = data.iloc[:,0].str.slice(stop=16).array
#adata_latent_space.obs["batch"] = data.iloc[:,0].str.slice(start=-4).array
#adata_latent_space.obs["cell_type"] = data.iloc[:,0].str.slice(start=17, stop=18).array

#################
#Solution data
os.chdir('/datos_2/phD_raw_data_2/OUTPUTS_2_ae_peaks_test_7_100pcent/q0/Summer_phase_1_data/starter_kit-joint_embedding-python/output/datasets_phase1v2/joint_embedding/openproblems_bmmc_multiome_phase1v2/')
solution = sc.read("openproblems_bmmc_multiome_phase1v2.censor_dataset.output_solution.h5ad")

os.chdir('/datos_2/phD_raw_data_2/OUTPUTS_2_ae_peaks_test_7_100pcent/q0/Summer_phase_1_data/starter_kit-joint_embedding-python/output/datasets_phase1v2/joint_embedding/openproblems_bmmc_cite_phase1v2/')
solution = sc.read("openproblems_bmmc_cite_phase1v2.censor_dataset.output_solution.h5ad")


#
#del adata_latent_space.obsm['Unnamed: 0']
#del adata_orig_rna.obsm['Unnamed: 0']
#scib.metrics.metrics(adata_latent_space, adata_orig_rna, batch_key='batch', label_key='label')



my_list = list()
###ALL IN ONE
#########################################################################################################
#########################################################################################################
## VIASH START
par = dict(
    input_prediction=adata_latent_space,
    input_solution=solution,
    output="openproblems_bmmc_multiome_starter.ari.had",
    debug=True
)

## VIASH END
print('Importing libraries')
import pprint
import scanpy as sc
import anndata
from scib.metrics import silhouette_batch
from scib.metrics import silhouette
from scib.metrics import cell_cycle
from scib.metrics import graph_connectivity
from scib.metrics.clustering import opt_louvain
from scib.metrics import nmi
from scib.metrics import trajectory_conservation

if par['debug']:
    pprint.pprint(par)

OUTPUT_TYPE = 'graph'
METRIC = 'asw_batch'

input_prediction = par['input_prediction']
input_solution = par['input_solution']
output = par['output']

print("Read prediction anndata")
adata = input_prediction
dataset_id = adata.uns['dataset_id']

print("Read solution anndata")
adata_solution = solution

print('Transfer obs annotations')
adata.obs['batch'] = adata_solution.obs['batch'][adata.obs_names]
adata.obs['cell_type'] = adata_solution.obs['cell_type'][adata.obs_names]

print('Preprocessing')
adata.obsm['X_emb'] = adata.X

print('Compute score')
sil_clus = silhouette_batch(
    adata,
    batch_key='batch',
    group_key='cell_type',
    embed='X_emb',
    verbose=False
)
#score = sil_clus['silhouette_score'].mean()
score = sil_clus
my_list.append(score)

##2
print('Preprocessing')
adata.obsm['X_emb'] = adata.X

print('Compute score')
score = silhouette(adata, group_key='cell_type', embed='X_emb')
my_list.append(score)

##3
organism = adata_solution.uns['organism']

recompute_cc = 'S_score' not in adata_solution.obs_keys() or \
               'G2M_score' not in adata_solution.obs_keys()

print('Preprocessing')
adata.obsm['X_emb'] = adata.X

print('Compute score')
score = cell_cycle(
    adata_pre=adata_solution,
    adata_post=adata,
    batch_key='batch',
    embed='X_emb',
    recompute_cc=recompute_cc,
    organism=organism
)
my_list.append(score)

##4
print('Preprocessing')
adata.obsm['X_emb'] = adata.X
sc.pp.neighbors(adata, use_rep='X_emb')

print('Compute score')
score = graph_connectivity(adata, label_key='cell_type')
my_list.append(score)

#5
print('Preprocessing')
adata.obsm['X_emb'] = adata.X
sc.pp.neighbors(adata, use_rep='X_emb')

print('Clustering')
opt_louvain(
    adata,
    label_key='cell_type',
    cluster_key='cluster',
    plot=False,
    inplace=True,
    force=True
)

print('Compute score')
score = nmi(adata, group1='cluster', group2='cell_type')
my_list.append(score)

##6
print('Preprocessing')
adt_atac_trajectory = 'pseudotime_order_ATAC' if 'pseudotime_order_ATAC' in adata_solution.obs else 'pseudotime_order_ADT'

adata.obsm['X_emb'] = adata.X
sc.pp.neighbors(adata, use_rep='X_emb')

print('Compute scores')
obs_keys = adata_solution.obs_keys()

if 'pseudotime_order_GEX' in obs_keys:
    score_rna = trajectory_conservation(
        adata_pre=adata_solution,
        adata_post=adata,
        label_key='cell_type',
        pseudotime_key='pseudotime_order_GEX'
    )
else:
    score_rna = np.nan

if adt_atac_trajectory in obs_keys:
    score_adt_atac = trajectory_conservation(
        adata_pre=adata_solution,
        adata_post=adata,
        label_key='cell_type',
        pseudotime_key=adt_atac_trajectory
    )
else:
    score_adt_atac = np.nan

score_mean = (score_rna + score_adt_atac) / 2
my_list.append(score_mean)
my_list









#########################################################################################################
#########################################################################################################
## VIASH START
par = dict(
    input_prediction=adata_latent_space,
    input_solution=solution,
    output="openproblems_bmmc_multiome_starter.ari.had",
    debug=True
)

## VIASH END

print('Importing libraries')
import pprint
import scanpy as sc
import anndata
from scib.metrics import silhouette_batch

if par['debug']:
    pprint.pprint(par)

OUTPUT_TYPE = 'graph'
METRIC = 'asw_batch'

input_prediction = par['input_prediction']
input_solution = par['input_solution']
output = par['output']

print("Read prediction anndata")
adata = input_prediction
dataset_id = adata.uns['dataset_id']

print("Read solution anndata")
adata_solution = solution

print('Transfer obs annotations')
#adata.obs['batch'] = adata_solution.obs['batch'][adata.obs['Unnamed: 0']].array
#adata.obs['cell_type'] = adata_solution.obs['cell_type'][adata.obs['Unnamed: 0']].array

print('Preprocessing')
adata.obsm['X_emb'] = adata.X

print('Compute score')
sil_clus = silhouette_batch(
    adata,
    batch_key='batch',
    group_key='cell_type',
    embed='X_emb',
    verbose=False
)
#score = sil_clus['silhouette_score'].mean()
score = sil_clus


#########################################################################################################
#########################################################################################################
# VIASH START
par = dict(
    input_prediction=adata_latent_space,
    input_solution=solution,
    output="openproblems_bmmc_multiome_starter.ari.had",
    debug=True
)

## VIASH END

print('Importing libraries')
import pprint
import scanpy as sc
import anndata
from scib.metrics import silhouette

if par['debug']:
    pprint.pprint(par)

OUTPUT_TYPE = 'graph'
METRIC = 'asw_label'

input_prediction = par['input_prediction']
input_solution = par['input_solution']
output = par['output']

print("Read prediction anndata")
adata = input_prediction
dataset_id = adata.uns['dataset_id']

print("Read solution anndata")
adata_solution = solution

print('Transfer obs annotations')
adata.obs['batch'] = adata_solution.obs['batch'][adata.obs['Unnamed: 0']].array
adata.obs['cell_type'] = adata_solution.obs['cell_type'][adata.obs['Unnamed: 0']].array

print('Preprocessing')
adata.obsm['X_emb'] = adata.X

print('Compute score')
score = silhouette(adata, group_key='cell_type', embed='X_emb')


#########################################################################################################
#########################################################################################################
## VIASH START
par = dict(
    input_prediction=adata_latent_space,
    input_solution=solution,
    output="openproblems_bmmc_multiome_starter.ari.had",
    debug=True
)
## VIASH END

debug = par['debug']

print('Importing libraries')
import pprint
import scanpy as sc
import anndata
import numpy as np
from scib.metrics import cell_cycle

if debug:
    pprint.pprint(par)

OUTPUT_TYPE = 'graph'
METRIC = 'cc_cons'

input_prediction = par['input_prediction']
input_solution = par['input_solution']
output = par['output']

print("Read prediction anndata")
adata = input_prediction
dataset_id = adata.uns['dataset_id']

print("Read solution anndata")
adata_solution = solution
organism = adata_solution.uns['organism']

print('Transfer obs annotations')
adata.obs['batch'] = adata_solution.obs['batch'][adata.obs['Unnamed: 0']].array
adata.obs['cell_type'] = adata_solution.obs['cell_type'][adata.obs['Unnamed: 0']].array
recompute_cc = 'S_score' not in adata_solution.obs_keys() or \
               'G2M_score' not in adata_solution.obs_keys()

print('Preprocessing')
adata.obsm['X_emb'] = adata.X

print('Compute score')
score = cell_cycle(
    adata_pre=adata_solution,
    adata_post=adata,
    batch_key='batch',
    embed='X_emb',
    recompute_cc=recompute_cc,
    organism=organism
)


#########################################################################################################
#########################################################################################################
## VIASH START
par = dict(
    input_prediction=adata_latent_space,
    input_solution=solution,
    output="openproblems_bmmc_multiome_starter.ari.had",
    debug=True
)

## VIASH END

print('Importing libraries')
import pprint
import scanpy as sc
import anndata
from scib.metrics import graph_connectivity

if par['debug']:
    pprint.pprint(par)

OUTPUT_TYPE = 'graph'
METRIC = 'graph_conn'

input_prediction = par['input_prediction']
input_solution = par['input_solution']
output = par['output']

print("Read prediction anndata")
adata = input_prediction
dataset_id = adata.uns['dataset_id']

print("Read solution anndata")
adata_solution = solution

print('Transfer obs annotations')
adata.obs['batch'] = adata_solution.obs['batch'][adata.obs['Unnamed: 0']].array
adata.obs['cell_type'] = adata_solution.obs['cell_type'][adata.obs['Unnamed: 0']].array

print('Preprocessing')
adata.obsm['X_emb'] = adata.X
sc.pp.neighbors(adata, use_rep='X_emb')

print('Compute score')
score = graph_connectivity(adata, label_key='cell_type')


#########################################################################################################
#########################################################################################################
## VIASH START
par = dict(
    input_prediction=adata_latent_space,
    input_solution=solution,
    output="openproblems_bmmc_multiome_starter.ari.had",
    debug=True
)

## VIASH END

print('Importing libraries')
import pprint
import scanpy as sc
import anndata
from scib.metrics.clustering import opt_louvain
from scib.metrics import nmi

if par['debug']:
    pprint.pprint(par)

OUTPUT_TYPE = 'graph'
METRIC = 'nmi'

input_prediction = par['input_prediction']
input_solution = par['input_solution']
output = par['output']

print("Read prediction anndata")
adata = input_prediction
dataset_id = adata.uns['dataset_id']

print("Read solution anndata")
adata_solution = solution

print('Transfer obs annotations')
adata.obs['batch'] = adata_solution.obs['batch'][adata.obs['Unnamed: 0']].array
adata.obs['cell_type'] = adata_solution.obs['cell_type'][adata.obs['Unnamed: 0']].array

print('Preprocessing')
adata.obsm['X_emb'] = adata.X
sc.pp.neighbors(adata, use_rep='X_emb')

print('Clustering')
opt_louvain(
    adata,
    label_key='cell_type',
    cluster_key='cluster',
    plot=False,
    inplace=True,
    force=True
)

print('Compute score')
score = nmi(adata, group1='cluster', group2='cell_type')


#########################################################################################################
#########################################################################################################
## VIASH START
par = dict(
    input_prediction=adata_latent_space,
    input_solution=solution,
    output="openproblems_bmmc_multiome_starter.ari.had",
    debug=True
)
## VIASH END

debug = par['debug']

print('Importing libraries')
import pprint
import numpy as np
import scanpy as sc
import anndata
from scib.metrics import trajectory_conservation

if debug:
    pprint.pprint(par)

OUTPUT_TYPE = 'graph'
METRIC = 'ti_cons'

input_prediction = par['input_prediction']
input_solution = par['input_solution']
output = par['output']

print("Read prediction anndata")
adata = input_prediction
dataset_id = adata.uns['dataset_id']

print("Read solution anndata")
adata_solution = solution

print('Transfer obs annotations')
adata.obs['batch'] = adata_solution.obs['batch'][adata.obs['Unnamed: 0']].array
adata.obs['cell_type'] = adata_solution.obs['cell_type'][adata.obs['Unnamed: 0']].array
adt_atac_trajectory = 'pseudotime_order_ATAC' if 'pseudotime_order_ATAC' in adata_solution.obs else 'pseudotime_order_ADT'

print('Preprocessing')
adata.obsm['X_emb'] = adata.X
sc.pp.neighbors(adata, use_rep='X_emb')

print('Compute scores')
obs_keys = adata_solution.obs_keys()

if 'pseudotime_order_GEX' in obs_keys:
    score_rna = trajectory_conservation(
        adata_pre=adata_solution,
        adata_post=adata,
        label_key='cell_type',
        pseudotime_key='pseudotime_order_GEX'
    )
else:
    score_rna = np.nan

if adt_atac_trajectory in obs_keys:
    score_adt_atac = trajectory_conservation(
        adata_pre=adata_solution,
        adata_post=adata,
        label_key='cell_type',
        pseudotime_key=adt_atac_trajectory
    )
else:
    score_adt_atac = np.nan

score_mean = (score_rna + score_adt_atac) / 2




