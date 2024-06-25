import os
import sys 
sys.path.append("Path/") #The parent directory of the current folder
import sICTA
import scanpy as sc
import numpy as np
import warnings 
import time
import pandas as pd
from sklearn.preprocessing import LabelEncoder
warnings.filterwarnings ("ignore")
import torch
print(torch.__version__)  
print(torch.cuda.get_device_capability(device=None),  torch.cuda.get_device_name(device=None))
import collections

our_query_adata = sc.read('Path/Muraro.h5ad')
current_path = os.getcwd()
print(our_query_adata)
sICTA.train(our_query_adata, gmt_path='human_gobp', label_name='broad_cell_type',epochs=50,current_path = current_path, project='hGOBP_demo',batch_size = 128)


# model_weight_path = './hGOBP_demo/model-49.pth_fin'
# adata, pre_list = sICTA.pre(our_query_adata, model_weight_path = model_weight_path,project='hGOBP_demo')
# new_adata = sc.AnnData(adata.X)
# new_adata.obs_names = adata.obs_names 
# new_adata.var_names = adata.var_names
# new_adata.obs['broad_cell_type_2'] = adata.obs['broad_cell_type_2']
# new_adata.write("/path/Muraro_att.h5ad")