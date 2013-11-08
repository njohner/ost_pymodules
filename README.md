This is a set of python modules that can be used to perform various
analyses on structures and dynamics of proteins and lipids.

The different modules and their functions are:

clustering.py:
=============

This module contains functions for performing distance based clustering.

-ClusterOnPairwiseDistance(view,dist_cutoff=3,prop_name='cluster'):
This function clusters the atoms in the view based on the pairwise distances.
A pair of atom with a distance lower than dist_cutoff will be part of the same cluster.
The clustering scales linearly with the number of atoms, unlike hierarchical clustering.

-HierarchichalClusteringOnPairwiseDistance(view,dist_cutoff,prop_name='cluster',method='single'):
This function uses the scipy hierarchical clustering to perform clustering of atoms based on their
pairwise distance. The same methids as in the scipy.cluster.hierarchy.linkage tollbox are available
('single', 'average', 'full').


