"""
.. codeauthor: Niklaus Johner <niklaus.johner@a3.epfl.ch>

This module contains functions for clustering

-ClusterOnPairwiseDistance(view,dist_cutoff=3,prop_name='cluster'):
This function clusters the atoms in the view based on the pairwise distances. 
A pair of atom with a distance lower than dist_cutoff will be part of the same cluster.
The clustering scales linearly with the number of atoms, unlike hierarchical clustering.

-HierarchichalClusteringOnPairwiseDistance(view,dist_cutoff,prop_name='cluster',method='single'):
This function uses the scipy hierarchical clustering to perform clustering of atoms based on their
pairwise distance. The same methods as in the scipy.cluster.hierarchy.linkage toolbox are available
('single', 'average', 'full').  

"""
from ost import *
import time
import numpy as npy
import os,math

def HierarchicalClusteringOnPairwiseDistance(view,dist_cutoff,prop_name='cluster',method='single'):
  import numpy as npy
  import scipy
  import scipy.cluster
  dist_mat=mol.alg.PariwiseDistanceMatrix(view,view)
  dm=[]
  for i,el in enumerate(dist_mat):
    for j,el2 in enumerate(el):
      if i>=j:continue
      dm.append(el2)
  dm=npy.array(dm)
  print 'starting clustering'   
  z=scipy.cluster.hierarchy.linkage(dm,method=method)
  m=scipy.cluster.hierarchy.fcluster(z,dist_cutoff,criterion='distance')
  c_list=[[] for i in range(max(m))]
  for i,el in enumerate(m):c_list[el-1].append(i)
  for a,c in zip(v.atoms,m):
    a.SetIntProp(prop_name,int(c))
  return (z,c_list)


def HoshenKopelman(neighbor_list):
  """
  This function performs a clustering using the Hoshen-Kopelman algorithm.
  It takes as input a list of neighbors, specifically for node i, 
  neighbor_list[i] is a list of all the neighbors of node i (a list of integers). So if
  neighbor_list[i]=[2,4,7] means that node i has nodes 2, 4, and 7 as neighbors.
  It returns a list of integers representing the cluster number for each node
  """
  node_list=[-1 for el in neighbor_list]
  node_list_labels=[]
  clusters=0
  for i,nl in enumerate(neighbor_list):
    #print i,node_list_labels,nl,node_list
    nl2=[node_list_labels[node_list[el]] for el in nl if node_list[el]!=-1]
    if len(nl2)==0:
      node_list_labels.append(clusters)
      node_list[i]=clusters
      clusters+=1
      continue
    nlp=min([node_list_labels[el] for el in nl2])
    node_list[i]=nlp
    #print 'nlp',nlp,node_list_labels
    for j in nl2:node_list_labels[j]=nlp
  #Correct the labels
  for i in range(len(node_list_labels)):
    N=i
    while (node_list_labels[N]<N):
      N=node_list_labels[N]
    node_list_labels[i]=N
  #Renumber labels
  ld={}
  clusters=1
  for el in node_list_labels:
    if not el in ld:
      ld[el]=clusters
      clusters+=1
  for i in range(len(node_list_labels)):
    node_list_labels[i]=ld[node_list_labels[i]]
  #Apply labels to array
  for i,el in enumerate(node_list):
    node_list[i]=node_list_labels[el]
  return node_list


def ClusterOnPairwiseDistance(view,dist_cutoff=3,prop_name='cluster'):
  an_l=[view.FindWithin(a.pos,dist_cutoff) for a in view.atoms]
  d={}
  for i,a in enumerate(view.atoms):d[a]=i
  an_l=[[d[a] for a in el] for el in an_l]
  cl=HoshenKopelman(an_l)
  for i in range(len(cl)):
    view.atoms[i].SetIntProp(prop_name,cl[i])
  return max(cl)

