This is a set of python modules that can be used to perform various analyses 
on structures and dynamics of proteins and lipids in OpenStructure. These modules should be used
with the latest stable release of OpenStructure (http://www.openstructure.org).

The main functionalities are presented below, but the complete documentation is included in docstrings in each module and function. It can be found online at http://njohner.github.io/ost_pymodules/ or it can be assembled into html files using sphinx from the source files. To build the documentation, simply go to the *doc* folder and :

```shell
make html
```
This should generate a set of html files in *doc/build/html*. To view the documentation, open the *index.html* file.
Some examples on how to use certain functions of the modules can be found in the *example* folder.

To import the modules in OpenStructure:

```python
import sys
sys.path.append(path_to_module)
import module
```
alternatively both OpenStructure (ost) and the modules can be imported in a python session (ost has to be first added to the PYTHONPATH for that):
```python
from ost import *
import sys
sys.path.append(path_to_module)
import module
```

The different modules and their functions are summarised below:

align_traj_on_density.py
====================

This module is used to align trajectories on a density. It was designed 
for complex lipidic phases such as the lipid cubic phase. It maximizes the overlap
of a selection of atoms with a reference density while minimizing the overlap
of a second selection with that same density. Typically the reference density is water.

lipid_analysis.py
====================

This module contains functions to determine lipid tilt and splay angles and
Calculate the elastic properties of the membrane from there.

If you use this code please cite ref. [1], where the method and its implementation 
is described in details or ref. [2], where the method is applied to complex lipidic phases.

To ensure proper treatment of the periodic boundary conditions, the trajectory
should first be extended to neighboring unit cells, and then aligned.
The tilts and splays are then calculated only for the central cell, the others 
being only used to ascertain correct treatment of the PBCs.

References
-------------
[1] Niklaus Johner, D. Harries and G. Khelashvili, 
       "Release note of the lipid tilt and splay code"

[2] Niklaus Johner, D. Harries and G. Khelashvili,
       "Curvature and lipid packing modulate the elastic properties of lipid assemblies: comparing the HII and lamellar phases."
       The Journal of Physical Chemistry Letters 5, no. 23 (2014), 4201-6.
 
lipid_order_parameters.py
======================
This module contains functions to calculate lipid order parameters form structures or trajectories.

density_alg.py
==============
This module is used to generate densities from trajectories and extract
surfaces from these densities in the form of sets of points.
Some algebra can be done on these surfaces, found in the surface_alg module

surface_alg.py
===============
This module is used to mainly to calculate curvatures of a point set surface
Such surfaces can be obtained from densities using functions form the density_alg module

entity_alg.py
=============
This module contains basic functions to work with structures (Entities and EntityViews).
It notably functions to generate the biounit or crystal packing from a pdb file

sequence_alg.py
=================
This module contains basic functions to work with sequences, notably to find motifs,
and build a position specific scoring matrix from an alignment.

hole.py
=========
This module contains functions to calculate the radius of a pore 
from a structure or a trajectory

principal_components.py
=======================
This module contains functions for principal components analysis. 
It can be used to calculate principal components from a trajectory,
graphically display the principal components, reconstruct a trajectory from
a chosen set of components and so on.

movies.py
==========
This module contains functions to render movies from OpenStructure using ffmpeg

secondary_structure.py
======================
This module contains functions to analyze the secondary structure of proteins (structures or trajectories).
It also has some functionalities to plot the secondary structure from a trajectory and to find secondary structure
elements in a protein.
It uses DSSP, which has to be installed and in the path for OpenStructure.

clustering.py:
=============
This module contains functions for performing distance-based clustering.
It can be used to perform hierarchical clustering or clustering using the
Hoshen-Kopelman algorithm.

pbc.py
=========
This module contains some functions to work with periodic boundary conditions. 

plot.py
===========
This module uses matplotlib and provides some high level plotting functions

angles.py
==========
This module contains functions to work with angles. It can notably be used
to correct for rotational symmetries of certain residues when comparing two structures
or different frames of a trajectory.

coarse_grained_conop.py
=========================
This module contains some functions to connect the atoms of coarse-grained martini structures.

hydrophobicity.py
===============
This module contains function to calculate hydrophobicity and hydrophobic
moments from different hydrophobic scales.

colvar.py
==========
This module is a collection of functions to work in conjunction with namd's colvar module.
It mainly contains functions to read som of the files output by namd.

file_utilities.py
=================
This module contains some helper functions for reading and writing files

