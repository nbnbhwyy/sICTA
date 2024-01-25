import sys 
sys.path.append("/mnt/rao/home/chenhg/R/Relate_method/TOSICA/") 
import sICTA
import scanpy as sc
import numpy as np
import warnings 
import pandas as pd
from sklearn.preprocessing import LabelEncoder
warnings.filterwarnings ("ignore")
import torch
print(torch.__version__)  

print(torch.cuda.get_device_capability(device=None),  torch.cuda.get_device_name(device=None))

import collections
our_query_adata = sc.read('/mnt/rao/home/chenhg/R/Relate_method/TOSICA/sICTA/data/Muraro.h5ad')

sICTA.train(our_query_adata, gmt_path='human_gobp', label_name='broad_cell_type',epochs=50,project='hGOBP_demo',batch_size = 128)

# model_weight_path = './hGOBP_demo/model-49.pth_fin'

# new_adata, pre_list = sICTA.pre(our_query_adata, model_weight_path = model_weight_path,project='hGOBP_demo')

# adata = sc.AnnData(new_adata.X)

# adata.obs_names = new_adata.obs_names
# adata.var_names = new_adata.var_names
# adata.obs['broad_cell_type_2'] = new_adata.obs['broad_cell_type_2']
# adata.write("/path/Muraro_att.h5ad")