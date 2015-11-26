.. Python modules for OpenStructure documentation master file, created by
   sphinx-quickstart on Mon May 18 15:50:06 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Python modules for OpenStructure
============================================================

.. toctree::
   :hidden:
   :maxdepth: 1

   align_traj_on_density
   angles
   clustering
   coarse_grained_conop
   colvar
   density_alg
   entity_alg
   file_utilities
   hole
   hydrophobicity
   lipid_analysis
   lipid_order_parameters
   movies
   pbc_utilities
   plot
   principal_components
   secondary_structure
   sequence_alg
   surface_alg
   trajectory_utilities
   installation.rst
   references.rst

Introduction
--------------------------------------
This is a set of python modules based on the `OpenStructure <http://www.openstructure.org>`_ framework and mainly aimed at performing various analyses of molecular dynamics simulations but also offers some functions for sequences and molecular structures.

The functions are roughly clustered into modules, which are rather strongly interdependent. These modules provide functions for classical analyses of MD trajectories such as principal components analysis, correlations (covariance, mutual information), helper functions to make nice movies, but also more advanced  analysis tools notably for simulations of lipidic systems. For example there are functions to :doc:`align frames <align_traj_on_density>` from a simulation of complex lipidic assemblies and calculate :doc:`the local curvatures <surface_alg>` of the membrane (see ref. [2]_ and [3]_), or to extract :doc:`the elastic moduli <lipid_analysis>` of the membrane from the simulation (ref. [4]_, [5]_ and [6]_).

Several modules from this set have already been moved into *OpenStructure*, as for example for calculating alpha helical kink angles (`mol.alg.helix_kinks <http://www.openstructure.org/docs/1.4/mol/alg/molalg/>`_), hbonds (mol.alg.hbond) or basic functions for structure and trajectory analysis (`mol.alg.structure_analysis <http://www.openstructure.org/docs/1.4/mol/alg/molalg/>`_ and `mol.alg.trajectory_analysis <http://www.openstructure.org/docs/1.4/mol/alg/molalg/>`_). Some of the modules or functions in the modules presented here might also migrate into *OpenStructure* in the future.

Installation
--------------------------------------

Be aware that these modules depend on `OpenStructure <http://www.openstructure.org>`_ and on several standard scientific python packages.
Help on installation and basic usage can be found :doc:`here <installation>`.

Available modules
--------------------------------------

* Several modules are aimed at analyzing lipidic systems.
   * :doc:`Aligning trajectories <align_traj_on_density>` without mtaching atoms of different frames. This is done using densities instead and is typically used to align membranes or complex lipidic systems.
   * Calculating :doc:`elastic properties <lipid_analysis>` of membranes.
   * Calculating :doc:`normals and curvatures <surface_alg>` on surfaces. 
   * A module to calculate :doc:`lipid order parameters <lipid_order_parameters>`
* Modules aimed at analyzing trajectories of proteins.
   * :doc:`Principal components analysis <principal_components>`
   * Other functions to :doc:`analyze trajectories <trajectory_utilities>` such as covariance, mutual information, correlations and so on.
   * Calculate and plot :doc:`protein secondary structure <secondary_structure>`, both for structures and trajectories. This module depends on the binding of openstructure to dssp.
   * Calculate the :doc:`size of a pore <pore>` in a channel or other protein.
* Graphical output:
   * Helper functions for :doc:`plotting <plot>`
   * Helper functions to make nice :doc:`movies <movies>`
* Some functions to work on :doc:`sequences and alignments <sequence_alg>` (e.g. search for motifs)
* Miscellaneous:
   * Helper functions to :doc:`read/write files <file_utilities>` and parse some information from PDB files.
   * Some functions to perform  :doc:`clustering <clustering>`.
   * Correctly build bonds for :doc:`Martini coarse-grained systems <coarse_grained_conop>`.
   * Some helper functions to work with :doc:`periodic boundary conditions <pbc>`.
   * Some support for reading :doc:`NAMD COLVAR <colvar>` files.


Cite Us
---------------------------------------
If you use any of this code, please cite ref. [1]_ for the use of OpenStructure.
Moreover if you use this code to align trajectories using the densities, please cite ref. [2]_, for calculations of curvature please cite ref. [3]_ and for calulations of lipid elastic properties please cite ref. [4]_. The theory behind the calculation of elastic properties and some applications can be found in refs. [5]_ and [6]_.

Moreover some of these modules also depend on other python libraries, notably
`numpy <http://www.numpy.org>`_, `scipy <http://www.scipy.org>`_ and `matplotlib <http://matplotlib.org>`_. If you use any function from this library in your work, please cite numpy (for example ref. [7]_). If you use any plots generated by these modules please cite matplotlib (ref. [8]_). Finally if you use these modules to align trajectories, calculate membrane elastic properties or for clustering please cite scipy (ref. [9]_). Other appropriate citations for these tools can be found on the `scipy website <http://www.scipy.org/citing.html>`_

.. [1] Marco Biasini, T. Schmidt, S. Bienert, V. Mariani, G. Studer, J. Haas, N. Johner, A.D. Schenk, A. Philippsen and T. Schwede, "OpenStructure: an integrated software framework for computational structural biology", Acta Cryst., 2013
.. [2] George Khelashvili, P. Blecua Carrillo Albornoz, N. Johner, S. Mondal, M. Caffrey, and H. Weinstein. "Why GPCRs Behave Differently in Cubic and Lamellar Lipidic Mesophases.", Journal of the American Chemical Society 134, no. 38 (September 26, 2012): 15858–68.
.. [3] Niklaus Johner, S. Mondal, G. Morra, M. Caffrey, H. Weinstein, and G. Khelashvili. "Protein and Lipid Interactions Driving Molecular Mechanisms of in Meso Crystallization." Journal of the American Chemical Society 136, no. 8 (February 26, 2014): 3271–84.
.. [4] Niklaus Johner, D. Harries and G. Khelashvili, "Implementation of a methodology for determining elastic properties of lipid assemblies from molecular dynamics simulations", submitted to BMC Bioinformatics (2015)
.. [5] George Khelashvili, N. Johner, G. Zhao, D. Harries, and H. L. Scott. "Molecular Origins of Bending Rigidity in Lipids with Isolated and Conjugated Double Bonds: The Effect of Cholesterol." Chemistry and Physics of Lipids 178 (February 2014): 18–26.
.. [6] Niklaus Johner, D. Harries, and G. Khelashvili. "Curvature and Lipid Packing Modulate the Elastic Properties of Lipid Assemblies: Comparing HII and Lamellar Phases." The Journal of Physical Chemistry Letters 5, no. 23 (December 4, 2014): 4201–6.
.. [7] Stéfan van der Walt, S. C. Colbert and G. Varoquaux. "The NumPy Array: A Structure for Efficient Numerical Computation", Computing in Science & Engineering, 13, 22-30 (2011)
.. [8] John D. Hunter. Matplotlib: A 2D Graphics Environment, Computing in Science & Engineering, 9, 90-95 (2007)
.. [9] Jones E, Oliphant E, Peterson P, et al. "SciPy: Open Source Scientific Tools for Python", 2001, http://www.scipy.org/

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

