#!/usr/bin/env python3
import os
import glob
import itertools
import numpy as np
import pandas as pd
import scanpy as sc
# import scvelo as scv
import networkx as nx

from typing import Set, List
from scipy import sparse
from natsort import natsorted

########################################################################################################################
# Variables
########################################################################################################################

random_seed = 12345

spatial_genes = [
    "Tmsb10",
    "Trh",
    "Cpm",
    "Ifitm1",
    "Nid2",
    "Afp",
    #     "Selenop",
    "Sepp1",
    "Ctsh",
    "Trap1a",
    "Fmr1nb",
    "Myl6b",
    "Gsn",
    "Krt19",
    "Peg10",
    "Wfdc2",
    "Atp1b1",
    "Tmem37",
    "Tmem120a",
    "Ino80c",
    "Sat1",
    "Lgmn",
    "Slc39a8",
    "Apoe",
    "S100a10",
    "Tagln2",
]

########################################################################################################################
## Velocity
########################################################################################################################

def create_db(metadata):
    db = {}
    for pb in metadata.Pool_barcode.unique():
        tmp = metadata.loc[metadata.Pool_barcode == pb].reset_index()
        db[pb] = dict(zip(tmp.Cell_barcode, tmp.Well_ID))
    return db

def hamming_distance(string1, string2):
    return np.sum([1 for idx, char in enumerate(string1) if string2[idx] != char ])

def process_velocity_batch(path: str, batch: str, metadata):
    """
    Function for processing velocity output from Batch.
    
    Arguments:
        path [str]: paht to velocity folder
        batch [str]: Batch
        metadata [DataFrame]
    """
    print(f'Processing {batch}')
    print('Reading matrices ...')
    mtx = np.loadtxt(f'{path}/Solo.out/Velocyto/raw/matrix.mtx', skiprows=3, delimiter=' ')
    mtx_shape = np.loadtxt(f'{path}/Solo.out/Velocyto/raw/matrix.mtx', skiprows=2, max_rows = 1, delimiter=' ')[0:2].astype(int)
    
    spliced = sparse.csr_matrix((mtx[:,2], (mtx[:,0]-1, mtx[:,1]-1)), shape = mtx_shape)
    unspliced = sparse.csr_matrix((mtx[:,3], (mtx[:,0]-1, mtx[:,1]-1)), shape = mtx_shape)
    ambiguous = sparse.csr_matrix((mtx[:,4], (mtx[:,0]-1, mtx[:,1]-1)), shape = mtx_shape)
    
    # load genes
    mtx_genes = pd.read_table(f'{path}/Solo.out/Velocyto/raw/features.tsv', header=None)
    mtx_genes.columns = ['ENSEMBL_ID', 'Symbol']
    mtx_genes = mtx_genes.set_index("Symbol")

    # load cells
    mtx_cells = pd.read_table(f'{path}/Solo.out/Velocyto/raw/barcodes.tsv', header=None, index_col=0)
    mtx_cells.index = mtx_cells.index.str.strip()
    
    print('Matching cells ...')
    db = create_db(metadata.query('Batch == @batch'))
    
    mtx_cells['Well_ID'] = 'EMPTY_CELL'
    # Step #1
    # Match only barcodes that are identical
    for barcode in mtx_cells.index:
        pb = barcode[:4].strip()
        cb = barcode[4:].strip()

        pb_cells = db.get(pb, 'EMPTY_POOL')
        if pb_cells != 'EMPTY_POOL':
            # does the cell barcode exist in db
            mtx_cells.loc[barcode, 'Well_ID'] = pb_cells.pop(cb, 'EMPTY_CELL')
        else:
            mtx_cells.loc[barcode, 'Well_ID'] = 'EMPTY_POOL'

    # Step #2
    correction_th = [1]
    for correction in correction_th:
        for barcode in mtx_cells[mtx_cells.Well_ID == 'EMPTY_CELL'].index:
            pb = barcode[:4].strip()
            cb = barcode[4:].strip()

            pb_cells = db.get(pb, 'EMPTY_POOL')
            if pb_cells != 'EMPTY_POOL':
                dist = np.array([hamming_distance(cb, x) for x in pb_cells])
                dist_idx = np.where(dist == correction)[0]
                if dist_idx.size > 1:
                    # found correction
                    key = list(pb_cells.keys())[dist_idx[0]]
                    mtx_cells.loc[barcode, 'Well_ID'] = pb_cells.pop(key)
            else:
                mtx_cells.loc[barcode, 'Well_ID'] = 'EMPTY_POOL'

    print(f'EMPTY_POOL: {sum(mtx_cells.Well_ID == "EMPTY_POOL")}')
    print(f'EMPTY_CELL: {sum(mtx_cells.Well_ID == "EMPTY_CELL")}')
    
    # match proper index to cell index from metadata
    # basically reorder the dataset but remember the index
#     mtx_cells = mtx_cells.reset_index()
#     mtx_cells = mtx_cells.iloc[pd.Index(mtx_cells['Well_ID']).get_indexer(metadata.index)]
#     mtx_cells = mtx_cells.reset_index().set_index('Well_ID')
    print('-------------------------')
    print('Summary')
    
    mtx_cells = mtx_cells.reset_index().set_index('Well_ID')
    
    spliced_ann = sc.AnnData(spliced.T, obs=mtx_cells, var=mtx_genes)
    spliced_ann.var_names_make_unique()
    spliced_ann.obs_names_make_unique()

    unspliced_ann = sc.AnnData(unspliced.T, obs=mtx_cells, var=mtx_genes)
    unspliced_ann.var_names_make_unique()
    unspliced_ann.obs_names_make_unique()

    ambiguous_ann = sc.AnnData(ambiguous.T, obs=mtx_cells, var=mtx_genes)
    ambiguous_ann.var_names_make_unique()
    ambiguous_ann.obs_names_make_unique()
    
    return [spliced_ann, unspliced_ann, ambiguous_ann]

def merge_velocity(adata, velocity):
    
    # find common cells between velocities matrices
    spliced, unspliced, ambiguous = [], [], []
    for idx, batch in enumerate(adata.obs.Batch.unique()):
        cells = adata.obs.query("Batch == @batch").index
        spliced.append(velocity[idx][0][cells])
        unspliced.append(velocity[idx][1][cells])
        ambiguous.append(velocity[idx][2][cells])
    
    # merge velocity batches
    spliced = sc.AnnData.concatenate(*spliced, uns_merge='unique', index_unique=None)
    unspliced = sc.AnnData.concatenate(*unspliced, uns_merge='unique', index_unique=None)
    ambiguous = sc.AnnData.concatenate(*ambiguous, uns_merge='unique', index_unique=None)
    
    common_genes = natsorted(np.intersect1d(adata.var_names, spliced.var_names))
    
    print(f'Adata: {adata.shape}')
    print(f'Common genes: {len(common_genes)}')
    
    adata_sub = adata[:, common_genes]
    print(f'Updated Adata: {adata_sub.shape}')
    
    adata_sub.layers['spliced'] = spliced[:, common_genes].X
    adata_sub.layers['unspliced'] = unspliced[:, common_genes].X
    adata_sub.layers['ambiguous'] = ambiguous[:, common_genes].X

    return adata_sub

def run_scvelo(adata):
    """Run basic workflow for computing velocities."""
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2_000)
    sc.tl.pca(adata, random_state=random_seed)
    sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30, random_state=random_seed)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    scv.tl.recover_dynamics(adata, n_jobs=8)
    scv.tl.velocity(adata, mode='dynamical')
    scv.tl.velocity_graph(adata)
    
    return adata

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

########################################################################################################################
# Spatial
########################################################################################################################
from read_roi import read_roi_zip
from shapely import geometry

def get_centroids_coord(rois):

    points = []
    for _, values in rois.items():
        poly = geometry.Polygon(list(zip(values["x"], values["y"])))
        points.append([poly.centroid.x, poly.centroid.y])

    return np.array(points)


def load_object(path: str) -> sc.AnnData:
    if not os.path.exists(path):
        print(f"Provided {path} does not exists!")

    counts = glob.glob(f"{path}/*.csv")
    if len(counts) != 1:
        print("Can't find a csv count matrix")

    counts = counts[0]

    adata = sc.AnnData(X=pd.read_table(counts, skiprows=1, index_col=0).T)
    adata = adata[adata.obs_names[:-1], :]
    rois = read_roi_zip(f"{path}/RoiSet.zip")
    adata.obsm["spatial"] = get_centroids_coord(rois)

    return adata

def nzm(adata):
    normalized = adata.copy()
    nzm = np.apply_along_axis(lambda v: np.median(v[np.nonzero(v)]), 0, adata.X)
    normalized.X /= nzm
    normalized.X = np.log1p(normalized.X)

    return normalized