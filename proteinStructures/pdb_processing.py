#!/usr/bin/env python
# coding: utf-8
## Created: 10/30/2024
## Author: Daniel Zeiberg
## Email: zeiberg.d@northeastern.edu

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1
import gzip
from scipy.spatial.distance import pdist,cdist
from sklearn.metrics import pairwise_distances
from multiprocessing import Pool,Array

parser = PDBParser()

def get_residue_coords(res):
    # Get the coordinates of the atoms of a residue
    coords = np.stack([np.array(list(a.get_vector())) for a in res.get_atoms()])
    return coords

def get_residue_pair_dist(a,b,atom_distances,corresp_res):
    # Distance between two residues
    a_idxs = np.where(corresp_res == a)[0]
    b_idxs = np.where(corresp_res == b)[0]
    return a,b,atom_distances[np.ix_(a_idxs,b_idxs)].min()


def read_pdb_residues(pdb_file_path : str,chain:str):
    """Get residues from pdb file
    Args:
        pdb_file_path (str): path to pdb file
        chain (str): chain to read [A,B,C,...] (usually A)
    """
    if pdb_file_path[-2:] == "gz":
        with gzip.open(pdb_file_path, "rt") as file_handle:
            structure = parser.get_structure("protein",file_handle)
    else:
        structure = parser.get_structure("protein",pdb_file_path)
    chains = {c.id: c for c in structure.get_chains()}
    return list(chains[chain].get_residues())


def get_residue_distance_mat(pdb_file_path,chain):
    """Get residue distance matrix from pdb file
    
    Args:
        pdb_file_path (str): path to pdb file
        chain (str): chain to read [A,B,C,...] (usually A)

    Returns:
        residueDistances (np.array): residue distance matrix
        residues (list): list of residues
    """
    residues = read_pdb_residues(pdb_file_path,chain)
    n_residues = len(residues)
    coord_arrays = [get_residue_coords(r) for r in residues]
    coords = np.concatenate(coord_arrays)
    corresp_res = np.concatenate([[i] * len(arr) for i,arr in enumerate(coord_arrays)])
    mm_corr = np.memmap('corresp_res.array', dtype=int, mode='write', shape=corresp_res.shape)
    mm_corr[:] = corresp_res[:]
    atom_distances = pairwise_distances(coords,n_jobs=-1)
    mm_ad = np.memmap('atom_distances.array',dtype='float32',mode="write",shape=atom_distances.shape)
    mm_ad[:] = atom_distances[:]
    mm_corr = np.memmap('corresp_res.array', dtype=int, mode='readonly', shape=corresp_res.shape)
    mm_ad = np.memmap('atom_distances.array',dtype='float32',mode="readonly",shape=atom_distances.shape)
    with Pool(8) as p:
        pairDists = p.starmap(get_residue_pair_dist,
                              [(a,b,mm_ad,mm_corr) for a in range(n_residues) \
                               for b in range(a+1,n_residues)])
    

    residueDistances = np.zeros((n_residues,n_residues))
    for p in pairDists:
        residueDistances[p[0],p[1]] = p[2]
        residueDistances[p[1],p[0]] = p[2]
    return residueDistances,residues

def write_adj_file(distances,output_path,threshold=6):
    """
    Write adjacency file from distance matrix
    Args:
        distances (np.array): distance matrix
        output_path (str): path to write output
        threshold (float): distance threshold for adjacency (default 6)

    """
    with open(output_path,"w") as f:
        for i,dists_i in enumerate(distances):
            neighbors = list(set(np.where(dists_i < threshold)[0]) - set([i]))
            lineout = '\t'.join([str(j) for j in [i] + neighbors])
            f.write(f"{lineout}\n")

def write_label_files(residues,posfile,negfile,labels_file):
    """ Write labels for residues (specific for graphlet counting input)"""
    with open(posfile,"w") as f:
        f.write("0\n")
    with open(negfile,"w") as f:
        for i in range(1,len(residues)):
            f.write(f"{i}\n")
    with open(labels_file,"w") as f:
        for r in residues:
            try:
                f.write(protein_letters_3to1[r.resname])
            except KeyError:
                f.write("X")
