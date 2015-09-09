"""
Example written by Niklaus Johner (niklaus.johner@a3.epfl.ch)

This is an example of how to align a trajectory for complex lipidic phases.
It should be run from within the "example" directory, otherwise the path
to the python modules will have to be set.

Everything that gets generated will be put in the "tmp" directory
"""
from ost import *
import os,sys

indir="./data"
outdir="./tmp"
sys.path.append("../")
import align_traj_on_density
import trajectory_utilities

#First we load the structure and trajectory
p=io.IOProfile(processor=conop.HeuristicProcessor(),dialect='CHARMM',fault_tolerant=True)
eh=io.LoadPDB(os.path.join(indir,"hex_phase.pdb"),profile=p)
t=io.LoadCHARMMTraj(eh,os.path.join(indir,'hex_phase.dcd'),stride=1)

#It is generally a good idea to keep only what is needed of the trajectory
#to reduce the number of atoms and speed up the calculations.
#For an all-atom simulation one could for example remove the hydrogen atoms:
#t=t.Filter(eh.Select("ele!=H"))
#eh=t.GetEntity()

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
#4. We extract the entity linked to the trajectory
#   and then set its positions to the first frame of the trajectory.
#   This way, when we save the trajectory and structure, the structure will have the positions of the first frame.
extended_eh=extended_t.GetEntity()
extended_t.CopyFrame(0)
io.SaveCHARMMTraj(extended_t,os.path.join(outdir,'hex_phase_extended.pdb'),os.path.join(outdir,'hex_phase_extended.dcd'),profile=p)

#######################################
#This trajectory is already aligned, so we will shift different frames around before realigning it
#Also we will keep fewer frames as alignment takes time.
extended_t=extended_t.Filter(extended_eh.Select(""),stride=5)
shift_list=geom.Vec3List([geom.Vec3(),geom.Vec3(0,2,0),geom.Vec3(0,4,0),geom.Vec3(0,8,0)])
trajectory_utilities.TranslateFrames(extended_t,shift_list)
io.SaveCHARMMTraj(extended_t,os.path.join(outdir,'hex_phase_shifted.pdb'),os.path.join(outdir,'hex_phase_shifted.dcd'),profile=p)

####################################
#Now we make the alignment
#1. As this is a coarse grained trajectory, atoms are not properly recognised so that masses
#   were not set when loading the structure. We need to set the masses and element to make sure
#   that the algorithm calculating the density for the alignment will work properly.
for a in extended_eh.atoms:
  a.SetElement('Ge')#Has the correct atomic weight
  a.SetRadius(2.35)
  a.SetMass(72)
#2. Now we make the alignment. This will take some time as it is an optimization procedure.
cm=extended_eh.Select("cname=A").GetCenterOfMass()
xl=align_traj_on_density.AlignTrajectoryOnFirstFrame(extended_t,'aname=W','cname=A and aname=W','cname=A and aname!=W',\
                                -cm,True,cm,os.path.join(outdir),"hex_phase_",0,False)




