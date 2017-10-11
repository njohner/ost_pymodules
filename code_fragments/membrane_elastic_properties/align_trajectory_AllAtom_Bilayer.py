"""
Code fragment written by Niklaus Johner (niklaus.johner@a3.epfl.ch)

This example presents the typical steps used to align a trajectory of a lipidic bilayer.
It simply translates the center of mass of the membrane on its position in the first frame.
There is no need for any rotation as the periodic boundary conditions impose
that the membrane plane is in the xy plane.

If you are working with a complex lipidic system (such as cubic or hexagonal phase),
you should consider using the optimization approach to align
your trajectory (see align_trajectory_AllAtom_ComplexLipidicPhase.py).

PBC are treated by replicating the simulation box (only in the x and y directions 
for bilayers). The new chains created for
this purpose will have the same chain name as the ones from the simulation box
with a different number added at the end of the name for each replicate of the box.
For example if you have a chain A, and we double the linear size of the system,
giving 8 replicates of the simulation box, you will end up with chains A,A1,A2,A3,A4,A5,A6,A7.
It is therefore important that all your chains have names shorter than 4 characters,
as otherwise you will get chain names with 5 characters which will not be handled properly
when saving and loading a PDB.

This script will output files to the outdir, replacing files with the same name if present.
The output files will all start with the "fname_basis"
The output files are:
-extended pdb and trajectory, which will starte with "fname_basis", complemented by "_extended".
-Aligned pdb and trajectory, which will starte with "fname_basis", complemented by "_aligned"..
"""
from ost import *
import os,sys
sys.path.append("../") #Set path to the python modules
import align_traj_on_density
import trajectory_utilities


indir="path/to/input directory"
outdir="path/to/output directory"
pdbname="pdbname"
trajname="trajname"
fname_basis="MySystem_"

#Selection used to reduce the system size and keep only useful part
filter_sele="ele!=H"
#Selection whose center of mass is used as center for the wrapping of the extended simulation.
wrap_sele="cname=A"
#Selection whose center of mass is used as center for the wrapping of the extended simulation.
align_sele="cname=A"
#Selection whose center of mass will be aligned.


#First we load the structure and trajectory
p=io.IOProfile(dialect='CHARMM',fault_tolerant=True)
eh=io.LoadPDB(os.path.join(indir,pdbname+".pdb"),profile=p)
t=io.LoadCHARMMTraj(eh,os.path.join(indir,trajname+".dcd"),stride=1)

#It is generally a good idea to keep only what is needed of the trajectory
#to reduce the number of atoms and speed up the calculations.
#For an all-atom simulation one could for example remove the hydrogen atoms.
t=t.Filter(eh.Select(filter_sele))
eh=t.GetEntity()


###################################
#We extend the trajectory
#1. We define the directions in which we want to extend the trajectory
#   Here we want to extend it in the x,y and x+y directions,
#   i.e end up with a system containing 4 unit cells. We do not extend in the z directions
#   as this is a planar bilayer and extending in the z direction is of no use for
#   calculation of membrane properties in that case.
extension_directions=[[1,0,0],[0,1,0],[1,1,0]]
#2. We extend the trajectory to its neighboring unit cells. The last argument is the mutliplicative
#   factors that will be applied to cell dimensions in the information in each frame
extended_t=trajectory_utilities.ExtendTrajectoryToNeighboringUnitCells(t,extension_directions,(2,2,1))
#3. We wrap the trajectory around the central unit cell
cm=mol.alg.AnalyzeCenterOfMassPos(t,eh.Select(wrap_sele))
trajectory_utilities.WrapTrajectoryInPeriodicCell(extended_t,cm)
#4. We extract the entity linked to the trajectory
#   and then set its positions to the first frame of the trajectory.
#   This way, when we save the trajectory and structure, the structure will have the positions of the first frame.
extended_eh=extended_t.GetEntity()
extended_t.CopyFrame(0)
io.SaveCHARMMTraj(extended_t,os.path.join(outdir,fname_basis+'extended.pdb'),os.path.join(outdir,fname_basis+'extended.dcd'),profile=p)


####################################
#Now we make the alignment.
cm=mol.alg.AnalyzeCenterOfMassPos(extended_t,extended_eh.Select(align_sele))
cm=geom.Vec3List([-el for el in cm])
trajectory_utilities.TranslateFrames(extended_t,cm)
io.SaveCHARMMTraj(extended_t,os.path.join(outdir,fname_basis+'aligned.pdb'),os.path.join(outdir,fname_basis+'aligned.dcd'),profile=p)





