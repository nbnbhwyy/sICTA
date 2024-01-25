# sICTA: Interpretable Cell Type Annotation based on self-training

The full description of sICTA and its application on published single cell RNA-seq datasets are available.

Download archive with preprocessed data at: https://drive.google.com/drive/folders/1jbqSxacL_IDIZ4uPjq220C9Kv024m9eL.

The repository includes detailed installation instructions and requirements, scripts and demos.


## 1 The workflow of sICTA.

![](Flow.jpg)

**(a)** Combining cell expression and marker gene specificity to generate pseudo-labels. **(b)**  The downstream Transformer classifiers are first pre-trained based on cell type probability distributions (pseudo-labels), followed by iterative refinement of the classifiers through a self-training framework until convergence. The sICTA takes the a priori knowledge from the biological domain and uses masked learnable embeddings to transform the input data ($G$ genes) into $k$ input tokens representing each gene set (GS) and a class token (CLS).
## 2 Requirements

+ Linux/UNIX/Windows system
+ Python == 3.8.6
+ torch == 1.12.1
+ scanpy == 1.9.1

Topic_gene_embedding

## 3 Usage

### Data format

sICTA requires cell-by-cell-gene matrix and cell type information to be entered in csv object format.

### Training

```bash
python main.py
```

We provide default data for users to understand and debug sICTA code.


## Reference

If you use `sICTA` in your work, please cite