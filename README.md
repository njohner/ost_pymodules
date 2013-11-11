This is a set of python modules that can be used to perform various analyses 
on structures and dynamics of proteins and lipids in OpenStructure.

The main functionalities are presented below, but more documentation is available
in the modules themselves. For most functions python's *help(foo)* function will display the
documentation for *foo*. This is especially useful in OpenStructure's interactive python shell
which has auto completion, making it easy to see the list of functions available ina  certain module.

To import the modules in OpenStructure:

```python
import sys
sys.path.append(path_to_module)
import module
```python

The different modules and their functions are:

clustering.py:
=============

This module contains functions for performing distance based clustering.

- ***ClusterOnPairwiseDistance(view,dist_cutoff=3,prop_name='cluster')***:
This function clusters the atoms in the view based on the pairwise distances.
A pair of atom with a distance lower than dist_cutoff will be part of the same cluster.
The clustering scales linearly with the number of atoms, unlike hierarchical clustering.

- ***HierarchichalClusteringOnPairwiseDistance(view,dist_cutoff,prop_name='cluster',method='single')***:
This function uses the scipy hierarchical clustering to perform clustering of atoms based on their
pairwise distance. The same methids as in the scipy.cluster.hierarchy.linkage tollbox are available
('single', 'average', 'full').

hbond.py:
========

This module is a flexible rewrite of HBPlus, allowing to calculate hbond
conservation between different structures or over a trajectory. 
It uses customizable dictionaries to define donors and acceptors and can
account for equivalent HBonds such as involving for example either oxygen atom
of an Aspartic Acid.

Main available functions are:

- ***GetHbondListFromView(eh,hbond_donor_acceptor_dict={},verbose=True)***:
Returns a list of hydrogen bonds from an Entity or EntityView.

- ***GetHbondListFromTraj(t,eh,cutoff=0.7,stride=1,swap=False,...)***:
Returns a list of hydrogen bonds from an Entity or EntityView that are present
present in a fraction of the frames larger than *cutoff*.

- ***GetHbondListBetweenViews(eh1,eh2,...)***:
Returns the list of hydrogen bonds formed between two Entity or EntityView.

- ***CalculateHBondScore(ref_eh,eh2,ref_eh2=None,...)***:
Returns the fraction of H-bonds from ref_eh that are also present in eh2.
If ref_eh2 is specified, it uses as reference the Hbonds between ref_eh and ref_eh2. 
This allows to look at H-bonds between specific parts of proteins or so.
Alternatively ref_eh can be a list of H-bonds.

- ***AnalyzeHBondScore(ref_eh,t,eh2,ref_eh2=None,...)***:
Returns the same score as CalculateHBondScore, but for a trajectory.


entity_alg:
==========

This module contains basic functions to work with structures (Entities and EntityViews).
It notably includes several functions to work with periodic conditions in orthogonal 
and non-orthogonal unit cells.

- ***FindWithinWithPBC(eh,pos,radius,cell_center,cell_size)***:

- ***FindWithinWithNonOrthogonalPBC(eh,pos,radius,ucell_vecs,vecs_to_neighbor_ucells=None)***:

- ***ExtendEntityToNeighboringUnitCells(eh,vecs_to_neighbor_ucells)***:

- ***ExtendEntityWithPBC(eh,cell_center,ucell_vecs,extension_size=10)***:

- ***GenerateCrystalPackingFromPDB(pdb_filename,distance_cutoff=20,superpose_sele='aname=CA')***:



