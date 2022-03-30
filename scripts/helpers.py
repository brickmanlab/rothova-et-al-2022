#!/usr/bin/env python3
import anndata
import os
import glob
import itertools
import numpy as np
import pandas as pd
import scanpy as sc
import networkx as nx

from typing import Set, List
from scipy import sparse
from natsort import natsorted

########################################################################################################################
# Variables
########################################################################################################################

random_seed = 12345

########################################################################################################################
## Trajectory inference
########################################################################################################################

def run_paga(adata, by='CellType'):
    sc.tl.paga(adata, groups=by)
    sc.pl.paga(adata, color=[by])
    sc.pl.paga(adata, threshold=0.1, show=False)
    sc.tl.draw_graph(adata, init_pos='paga')
    sc.pl.draw_graph(adata, color=by)

    return adata

########################################################################################################################
## Velocity
########################################################################################################################

def run_scvelo(adata):
    import scvelo as scv
    """Run basic workflow for computing velocities."""
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    sc.tl.pca(adata, random_state=random_seed)
    sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30, random_state=random_seed)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    scv.tl.recover_dynamics(adata, n_jobs=8)
    scv.tl.velocity(adata, mode='dynamical')
    scv.tl.velocity_graph(adata)
    
    return adata


# taken from https://github.com/alexdobin/STAR/issues/774
def buildAnndataFromStar(path):
    """Generate an anndata object from the STAR aligner output folder"""
    path=path
    # Load Read Counts
    X = sc.read_mtx(path+'Gene/raw/matrix.mtx')
    
    # Transpose counts matrix to have Cells as rows and Genes as cols as expected by AnnData objects
    X = X.X.transpose()
    
    # This matrix is organized as a sparse matrix with Row, Cols and 3 values columns for 
    # Spliced, Unspliced and Ambigous reads
    mtx = np.loadtxt(path+'Velocyto/raw/matrix.mtx', skiprows=3, delimiter=' ')
    # Extract sparse matrix shape informations from the third row
    shape = np.loadtxt(path+'Velocyto/raw/matrix.mtx', skiprows=2, max_rows = 1 ,delimiter=' ')[0:2].astype(int)

    # Read the sparse matrix with csr_matrix((data, (row_ind, col_ind)), shape=(M, N))
    # Subract -1 to rows and cols index because csr_matrix expects a 0 based index
    # Traspose counts matrix to have Cells as rows and Genes as cols as expected by AnnData objects
    spliced = sparse.csr_matrix((mtx[:,2], (mtx[:,0]-1, mtx[:,1]-1)), shape = shape).transpose()
    unspliced = sparse.csr_matrix((mtx[:,3], (mtx[:,0]-1, mtx[:,1]-1)), shape = shape).transpose()
    ambiguous = sparse.csr_matrix((mtx[:,4], (mtx[:,0]-1, mtx[:,1]-1)), shape = shape).transpose()
    
    # Load Genes and Cells identifiers
    obs = pd.read_csv(path+'Velocyto/raw/barcodes.tsv',
                  header = None, index_col = 0)
    # Remove index column name to make it compliant with the anndata format
    obs.index.name = None
    
    var = pd.read_csv(path+'Velocyto/raw/features.tsv', sep='\t',
                  names = ('gene_ids', 'feature_types'), index_col = 1)
    
    # Build AnnData object to be used with ScanPy and ScVelo
    adata = anndata.AnnData(X = X, obs = obs, var = var,
                        layers = {'spliced': spliced, 'unspliced': unspliced, 'ambiguous': ambiguous})
    adata.var_names_make_unique()
    
    # Subset Cells based on STAR filtering
    selected_barcodes = pd.read_csv(path+'Gene/filtered/barcodes.tsv', header = None)
    adata = adata[selected_barcodes[0]]
    
    return adata.copy()
