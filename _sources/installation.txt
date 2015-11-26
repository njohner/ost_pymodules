
Installation
============================================================

Dependencies
--------------

These python modules all depend on the computational framework
`OpenStructure <http://www.openstructure.org>`_ (*ost*), so that you will first
have to install OpenStructure either by `downloading an available binary <http://www.openstructure.org/download/>`_ or `from source <http://www.openstructure.org/docs/1.4/install/>`_. 


Moreover some of these modules also depend on other python libraries, notably
`numpy <http://www.numpy.org>`_, `scipy <http://www.scipy.org>`_ and `matplotlib <http://matplotlib.org>`_. If you are installing `OpenStructure <http://www.openstructure.org>`_ from source, you should therefore make sure to first add these modules to your python installation.



The modules
---------------

These modules are installed and used as any python module. Simply download the modules and either add the *path/to/the/python/files* to your *PYTHONPATH*, for example in your *.bashrc*:
::
  export PYTHONPATH=path/to/the/python/files:$PYTHONPATH

Alternatively you can add the path directly in python:
::
  import sys
  sys.path.append("path/to/the/python/files")


Usage
============================================================

Once the path to the python modules has been set, you can simply import them as any other python module, e.g:
::
  import sys
  sys.path.append("path/to/the/python/files") 
  import align_traj_on_density
  help(align_traj_on_density)

As these modules depend on OpenStructure (*ost*), you will also need to make sure your python knows about *ost*. Several options are available here. First *ost* can be seen as a collection of python modules, which can simply be imported in python, provided that you add *path/to/ost/stage/lib/python2.7/site-packages* either to your *PYTHONPATH* or directly in python to your system path as described above. Second, *ost* can be used as a standalone software, either with (*path/to/ost/stage/bin/dng*) or without (*path/to/ost/stage/bin/ost*) graphical user interface. Running one of these binaries will basically give you a python interpreter with all the *ost* modules already imported (see the `ost documentation <http://www.openstructure.org>`_ for more information).

A typical usage case could therefore look something like:
::
  import sys
  from ost import io
  sys.path.append("path/to/the/python/files") 
  import align_traj_on_density as align
  
  eh=io.LoadPDB("path/to/pdb")
  t=io.LoadCHARMMTraj(eh,"path/to/dcd")
  xl=align.AlignTrajectoryOnFirstFrame(t,eh.Select("rname=TIP3"))

As illustrated in the above example, these python modules mainly take as inputs *ost* objects, so that you will need to have or aquire prior knowledge of *OpenStructure* to use these modules. A good place to start is the `OpenStructure <http://www.openstructure.org>`_ website, which has tutorials and documentation. Moreover, these python modules are not meant as a fool-proof software, rather as a collection of functions, so that you will absolutely need to know python and write more than one line of code to use them.

Nevertheless, we provide examples for the use of many of the functions provided in this set of modules, which can be found in the *examples* subdirectory of this distribution. In an interactive python session, the *help* function will give you access to most of this documentation.



