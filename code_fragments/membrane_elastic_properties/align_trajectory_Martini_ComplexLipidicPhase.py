"""
Code fragment written by Niklaus Johner (niklaus.johner@a3.epfl.ch)

This example presents the typical steps used to align a trajectory of a complex
lipidic phase. This is a slow process as an optimization is required for each frame.
If you are working with a lipid bilayer, there are much faster ways of aligning 
your trajectory (see align_trajectory_bilayer.py).

PBC are treated by replicating the simulation box. The new chains created for
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
-extended pdb and trajectory, which will have the same names as the loaded trajectory and pdb,
complemented by "_extended".
-Aligned pdb and trajectory, which will have the same names as the loaded trajectory and pdb,
complemented by "_aligned".
-several other files...
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


#Selection whose center of mass is used as center for the wrapping of the extended simulation.
wrap_sele="cname=A"
#selection used to calculate the density map of waters.
#This should contain all the water molecules of the extended simulation box
density_sele="rname=W" 
#selection for which the overlap with the water density will be maximised.
#This should contain only the water molecules from the initial simulation box
water_sele2="rname=W and cname=A"
#selection for which the overlap with the water density will be minimised.
#This should contain everything in the initial simulation box except the water molecules.
memb_sele="rname!=W and cname=A"


#First we load the structure and trajectory
p=io.IOProfile(processor=None,dialect='CHARMM',fault_tolerant=True)
eh=io.LoadPDB(os.path.join(indir,pdbname+".pdb"),profile=p)
t=io.LoadCHARMMTraj(eh,os.path.join(indir,trajname+".dcd"),stride=1)

# As this is a coarse grained trajectory, atoms are not properly recognised so that masses
# were not set when loading the structure. We need to set the masses and element to make sure
# that center of masses as well as the algorithm calculating the density used to determine 
# the membrane interface will work properly.
for a in eh.atoms:
  a.SetElement('Ge')#Has the correct atomic weight
  a.SetRadius(2.35)
  a.SetMass(72)

###################################
#We extend the trajectory
#1. We get the unit cell vectors
vecs_list=[t.GetFrame(i).GetCellVectors() for i in range(t.GetFrameCount())]
#2. vecs_list contains for each frame a Vec3List with 3 Vec3 in it corresponding
#   to the three unit cell vectors. We add the vectors to the other neighboring unit cells
#   for each frame.
for vl in vecs_list:vl.extend([vl[0]+vl[1],vl[0]+vl[2],vl[1]+vl[2],vl[0]+vl[1]+vl[2]])
#3. We extend the trajectory to its neighboring unit cells. The last argument is the mutliplicative
#   factors that will be applied to cell dimensions in the information in each frame
extended_t=trajectory_utilities.ExtendTrajectoryToNeighboringUnitCells(t,vecs_list,(2,2,2))
#4. We wrap the trajectory around the central unit cell
cm=mol.alg.AnalyzeCenterOfMassPos(t,eh.Select(wrap_sele))
trajectory_utilities.WrapTrajectoryInPeriodicCell(extended_t,cm)
#5. We extract the entity linked to the trajectory
#   and then set its positions to the first frame of the trajectory.
#   This way, when we save the trajectory and structure, the structure will have the positions of the first frame.
extended_eh=extended_t.GetEntity()
extended_t.CopyFrame(0)
io.SaveCHARMMTraj(extended_t,os.path.join(outdir,fname_basis+'extended.pdb'),os.path.join(outdir,fname_basis+'extended.dcd'),profile=p)


####################################
#Now we make the alignment. This will take some time as it is an optimization procedure.
cm=extended_eh.Select(wrap_sele).GetCenterOfMass()
xl=align_traj_on_density.AlignTrajectoryOnFirstFrame(extended_t,density_sele,water_sele,memb_sele,\
                                -cm,True,cm,outdir,fname_basis,0,False)




