import logging
import anndata as ad
import numpy as np
import pandas as pd
import leidenalg
import pkg_resources
import tensorflow as tf 
import keras
import scanpy as scanpy
import scipy
import os
import multiprocessing as mp
import itertools
import gc
from sklearn.model_selection import GridSearchCV

from keras.models import Model
from keras import regularizers, activations, initializers, constraints
from keras.layers import InputSpec
from keras.layers import LeakyReLU
from keras.layers import Dense, Dropout
from keras.layers import Layer
from keras.callbacks import EarlyStopping, ReduceLROnPlateau
from keras.constraints import UnitNorm
from keras.constraints import Constraint
from keras import backend as K
from keras.wrappers.scikit_learn import KerasClassifier
from keras.wrappers.scikit_learn import KerasRegressor
from tensorflow.keras.layers import BatchNormalization
from tensorflow.keras.optimizers import Adam

logging.basicConfig(level=logging.INFO)

## IMPORT DATA
dataset_path = "your_path_to_h5ad_data"

par = {
    'input_mod1': dataset_path + 'mod1.h5ad',
    'input_mod2': dataset_path + 'mod2.h5ad',
}

## READ DATA
logging.info('Reading `h5ad` files...')
ad_mod1 = ad.read_h5ad(par['input_mod1'])
ad_mod2 = ad.read_h5ad(par['input_mod2'])

## GENERATING VARIABLES ENVIROMENT
mod1_obs = ad_mod1.obs
mod1_uns = ad_mod1.uns

#CHANGE DATA FORMATS
ad_mod1_raw = ad_mod1.X.toarray()
ad_mod1_raw_matrix = pd.DataFrame(data=ad_mod1_raw, index=ad_mod1.obs_names, columns=ad_mod1.var_names)

ad_mod2_raw = ad_mod2.X.toarray()
ad_mod2_raw_matrix = pd.DataFrame(data=ad_mod2_raw, index=ad_mod2.obs_names, columns=ad_mod2.var_names)

#scanpy for MVG obtain
adata = scanpy.AnnData(ad_mod1_raw_matrix)
scanpy.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var.highly_variable]
np.array(adata.var_names)

#filtered matrix for RNA base on 2000 MVG
ad_mod1_raw_matrix= ad_mod1_raw_matrix[np.array(adata.var_names)]

#GET DIMENSIONS FROM REMAINING DATA
ad_mod1_raw_matrix_x = ad_mod1_raw_matrix.values[:,0:(ad_mod1_raw_matrix.shape[1])]
ad_mod2_raw_matrix_x = ad_mod2_raw_matrix.values[:,0:(ad_mod2_raw_matrix.shape[1])]

ncol_mod_1 = ad_mod1_raw_matrix_x.shape[1]
ncol_mod_2 = ad_mod2_raw_matrix_x.shape[1]


###############
###############
###############
#AUTO_fine_tune

#MODEL dynamic generator base on hyperparameters
def createmodel(n_layers, n_nodes, alpha, dropout, batch_size, mid_layer):
  model = keras.Sequential()
  for i in range(0, n_layers):
    if i==0:
    #Input_layer
      model.add(Dense(n_nodes, input_dim=ncol_mod_2))
      model.add(LeakyReLU(alpha=alpha))
      model.add(Dropout(dropout))
    #intermediate_encoder_layers
    else:
      model.add(Dense(n_nodes/(2*i)))
      model.add(LeakyReLU(alpha=alpha))
      model.add(Dropout(dropout))
  #middle_layer
  model.add(Dense(mid_layer, name = 'Bottleneck'))
  model.add(LeakyReLU(alpha=alpha))
  model.add(Dropout(dropout))
    #intermediate_decoder_layers    
  for i in range(n_layers-1, -1, -1):
    if i!=0:
      model.add(Dense(n_nodes/(2*i)))
      model.add(LeakyReLU(alpha=alpha))
      model.add(Dropout(dropout))   
    else:
      model.add(Dense(n_nodes))
      model.add(LeakyReLU(alpha=alpha))
      #model.add(Dropout(dropout))          
  #Output_layer
  model.add(Dense(ncol_mod_1, activation='relu'))
  model.compile(optimizer = Adam(lr=0.001), loss='mse', metrics=['mse'])
  return model

########
#FIX THE GRID TO BE USED BY THE MODEL
n_layers = [1,2,3,4,5,6] 
n_nodes = [128,256,512,1024,2048] 
alpha = [0.05,0.1,0.3,0.5] 
dropout = [0.1,0.2,0.3,0.4]  
batch_size = [7000] 
mid_layer = [10,30,50,70] 

########
#TRAIN MODELS IN PARALLEL
def fast_tune(n_layers,n_nodes,alpha,dropout,batch_size,mid_layer):
  model = createmodel(n_layers,n_nodes,alpha,dropout,batch_size,mid_layer)
  print(model.summary())
  #compile
  model.compile(optimizer = Adam(lr=0.001), loss='mean_squared_error')
  #callbacks
  early_stop = EarlyStopping(monitor='loss', mode='min', verbose=1 ,patience=20)
  lr_plateau_callback = ReduceLROnPlateau(monitor="loss", factor=0.1, patience=15, min_lr=0.00001, verbose = 1)
  #fit
  trained_model= model.fit(ad_mod2_raw_matrix_x, ad_mod1_raw_matrix_x  , epochs=1500, batch_size=7000, validation_split = 0, verbose = 1, shuffle = True, callbacks=[early_stop, lr_plateau_callback]) 
              
  #Retrieve middle_layer from trained_model
  layer_name = 'Bottleneck'
  intermediate_layer_model = Model(inputs=model.input,
                                   outputs=model.get_layer(layer_name).output)
  intermediate_output = intermediate_layer_model.predict(ad_mod2_raw_matrix_x)
  
  #Esport to csv the embedding "LIBRA_components"
  middle_layer =pd.DataFrame(intermediate_output)
  middle_layer.to_csv(dataset_path
  + 'n_layers' + str(n_layers) 
  + '_' + 'n_nodes' + str(n_nodes) 
  + '_' + 'alpha' + str(alpha) 
  + '_' + 'dropout' + str(dropout) 
  + '_' + 'batch_size' + str(batch_size) 
  + '_' + 'mid_layer' + str(mid_layer) 
  + '.csv', index=True)
  del(model)
  del(trained_model)
  del(intermediate_layer_model)
  del(intermediate_output)
  del(middle_layer)
  tf.keras.backend.clear_session()
  gc.collect()


a = [n_layers,n_nodes,alpha,dropout,batch_size,mid_layer]
input = list(itertools.product(*a))

n_jobs  = 5 #Base on datset size and computational power range it from 1(serie) to (>1) to train n models in parallel
pool = multiprocessing.Pool(processes=n_jobs)
pool.starmap(fast_tune, input)

