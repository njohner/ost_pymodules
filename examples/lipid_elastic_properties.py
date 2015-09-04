"""
Example written by Niklaus Johner (niklaus.johner@a3.epfl.ch)

This is an example of how to calculate the membrane elastic properties from a simulation.
We start here with an aligned trajectory as making the alignment (see the align_trajectory.py example)
is slow.

This example should be run from within the "example" directory, otherwise the path
to the python modules and other files will have to be set.
Everything that gets generated will be put in the "tmp" directory
"""
from ost import *
import os,sys
import numpy as npy
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

indir="./data"
outdir="./tmp"
sys.path.append("../")
import lipid_analysis,trajectory_utilities

#First we load the structure and trajectory
p=io.IOProfile(dialect='CHARMM',fault_tolerant=True)
eh=io.LoadPDB(os.path.join(indir,"hex_phase.pdb"),profile=p)
t=io.LoadCHARMMTraj(eh,os.path.join(indir,'hex_phase.dcd'),stride=1)


####################################
# As this is a coarse grained trajectory, atoms are not properly recognised so that masses
# were not set when loading the structure. We need to set the masses and element to make sure
# that the algorithm calculating the density used to determine the membrane interface will work properly.
for a in eh.atoms:
  a.SetElement('Ge')#Has the correct atomic weight
  a.SetRadius(2.35)
  a.SetMass(72)

####################################
# For this example the trajectory used was already aligned, but it contains only one unit cell
# To properly treat boundary conditions we replicate the unit cell all around the original one.
# And then wrap it around the central unit cell
# See the aign_trajectory.py example for more information.
vecs_list=[t.GetFrame(i).GetCellVectors() for i in range(t.GetFrameCount())]
for vl in vecs_list:vl.extend([vl[0]+vl[1],vl[0]+vl[2],vl[1]+vl[2],vl[0]+vl[1]+vl[2]])
t=trajectory_utilities.ExtendTrajectoryToNeighboringUnitCells(t,vecs_list,(2,2,2))
eh=t.GetEntity()
centers=mol.alg.AnalyzeCenterOfMassPos(t,eh.Select("cname=A"))
trajectory_utilities.WrapTrajectoryInPeriodicCell(t,centers)

#We save the trajectory for visualization
io.SaveCHARMMTraj(t,os.path.join(outdir,"extended.pdb"),os.path.join(outdir,"extended.dcd"),profile=p)

###################################
# We prepare the dictionaries with information about the selections defining headgroups and tails
# for the different lipids. They will be used in the calculation of tilts and splays.

# 1. The system in the example has two types of lipids: cholesterol (CHO) and DOPE (DOP)
#    We need a list with residue names of all lipids that will be analysed
lipid_names=['CHO','DOP']

# 2. Define the residue name of the water molecules.
#    This is used together with the lipid names to calculate the density maps
#    Of water and lipids used to extract the membrane interface and orient the normals on the interface.
water_name='W'

# 3. For each lipid type, we define the selection for headgroups and tails and the selection
#    that will be used to calculate distances between lipids. This should be atoms lying at
#    the neutral plane.
head_group_dict={'DOP':'aname=PO4,GL1,GL2','CHO':'aname=R1'}
tail_dict={'DOP':'aname=C5A,C5B','CHO':'aname=R5'}
distance_sele_dict={'DOP':'aname=C1*','CHO':'aname=R1'}

# 4. It is generally a good idea to keep only what is needed of the trajectory
#    to reduce the number of atoms and speed up the calculations.
sele="rname={0}".format(water_name)
for l in lipid_names:sele+=" or (rname={0} and {1})".format(l,head_group_dict[l])
for l in lipid_names:sele+=" or (rname={0} and {1})".format(l,tail_dict[l])
for l in lipid_names:sele+=" or (rname={0} and {1})".format(l,distance_sele_dict[l])
t=t.Filter(eh.Select(sele))
eh=t.GetEntity()

# 5. We define the different cutoffs used by the algorithms. More information can be found
#    In the documentation
#max distance between two lipids to be considered in the splay calculation
distance_cutoff=10.0 
#max tilt angle with respect to normal for splay calculations
angle_cutoff=0.175 
# radius of area used to calculate the normals to the membrane interface. If too small, normals will be noisy.
within_size_normals=10.0
#Density cutoff used when calculating the interface.
#This can be set to 0, but calculations will be slower.
density_cutoff=0.3
#Stride used when calculating the water and lipid densities.
#When set to 1, all frames are considered, which can be slow.
density_stride=1

#6. Periodic boundaries are best treated by replicating the simulation box around
#   the original unit cell and then calculating tilts and splays only for lipids from
#   the original unit cell, but using the surrounding unit cells to find all neighbors 
#   of a lipid for the splay calculation and to avoid boundary effects on the interfaces
#   and hence on the normal vectors used both for the tilt and splay calculations.
#   We therefore need to tell the function which lipids belong to the central unit cell.
#   This is done by setting a bool property for the corresponding residues.

tilt_bool_prop='do_tilt'
splay_bool_prop='do_splay'
v=eh.Select('cname=A and rname=DOP,CHO')
for r in v.residues:
  r.SetBoolProp(tilt_bool_prop,True)
  r.SetBoolProp(splay_bool_prop,True)
v=eh.Select('cname!=A and rname=DOP,CHO')
for r in v.residues:
  r.SetBoolProp(tilt_bool_prop,False)
  r.SetBoolProp(splay_bool_prop,False)

#7. Other parameters
prot_sele=None
sele_dict={}
filename_basis='tilt&splay_'
#8. We calculate the tilts and splays:
(lipid_tilt_dict,lipid_normal_dict,splay_dict,b_eh)=lipid_analysis.AnalyzeLipidTiltAndSplay(t,
  lipid_names,head_group_dict,tail_dict,distance_cutoff,within_size_normals,distance_sele_dict,water_name,
  outdir,density_cutoff,prot_sele,density_stride,tilt_bool_prop,splay_bool_prop,filename_basis,sele_dict)


########################################################
# Now we analyze the lipid tilts and splays, fitting the corresponding analytical functions
# to extract the elastic constants. See documentation and papers cited therein for more information.
# We will make one fit for each lipid type, and then calculate the overall tilt modulus
# By taking a weighted average of the individual contributions.

# The number of bins used when making the histograms for the tilts
nbins=100

# Name of the file where the tilt constants will be written
outfile=open(os.path.join(outdir,'tilt_constants.txt'),'w')

# Lists where we will store the tilt moduli, the associated uncertainty and the number of lipids of each type.
k_list=FloatList()
deltak_list=FloatList()
nl_list=FloatList()
for lipid_name in lipid_names:
  fname=lipid_name
  
  #We make a list containing all the tilts for one type of lipids.
  #In the dictionary, there is one entry for each lipid so we have to concatenate all of them into one list.
  tilt_list=FloatList()
  for el in lipid_tilt_dict['all'][lipid_name][0]:tilt_list.extend(el)

  #Now we make the fit. This will generate plots written to outdir.
  k,dk,kl=lipid_analysis.FitTiltDistribution(tilt_list,outdir=outdir,filename_basis=fname,title_complement='for '+lipid_name,nbins=nbins)
  print "Tilt modulus for {0} is k={1:1.2f} +/- {2:1.2f}.".format(lipid_name,k,dk)

  #We write out the tilt modulus into the file
  outfile.write(' '.join([lipid_name,str(round(kl[0],1)),str(round(dk,1))])+'\n')
  
  #Keep the tilt modulus, error and number of lipids of that kind to calculate the overall modulus later.
  k_list.append(k)
  deltak_list.append(dk)
  nl_list.append(len(tilt_list))

#Now we calculate the overall tilt modulus
ntot=npy.sum(nl_list)
k=npy.sum([ni/ki for ni,ki in zip(nl_list,k_list)])/ntot
k=1./k
dk=npy.sum([((k**2.0)/(ki**2.0)*ni/ntot)**2.0*(dki**2.0) for ni,ki,dki in zip(nl_list,k_list,deltak_list)])
dk=npy.sqrt(dk)
#write it out
outfile.write(' '.join(['agregated',str(round(k,1)),str(round(dk,1))])+'\n')
outfile.close()
#Close the generated plots
for i in range(10):plt.close()

# Now we fit the splays to extract the bending rigidities.
# The procedure is the same as for the tilts
# Except that we need the area per lipid to calculate the bending rigidity from the splays
# In this example the area per lipid is 49.5 A^2. For a lipid bilayer the area can be obtained
# Using the AnalyzeAreaPerLipid function.
lipid_area=49.5

outfile=open(os.path.join(outdir,'splay_constants.txt'),'w')
k_list=FloatList()
nl_list=FloatList()
deltak_list=FloatList()
nbins=100
for key in splay_dict['all']:
  fname=key
  splay_list=[el[0] for el in splay_dict['all'][key]]
  w=npy.std(splay_list)
  x_range=[-4*w,4*w]
  k,dk,kl=lipid_analysis.FitSplayDistribution(splay_list,lipid_area,outdir=outdir,filename_basis=fname,title_complement='for '+key,nbins=nbins,x_range=x_range)
  print "Splay modulus for {0} is k={1:1.2f} +/- {2:1.2f}.".format(key,k,dk)
  outfile.write(' '.join([key,str(round(kl[0],1)),str(round(dk,1))])+'\n')
  k_list.append(k)
  deltak_list.append(dk)
  nl_list.append(len(splay_list))
ntot=npy.sum(nl_list)
k=npy.sum([ni/ki for ni,ki in zip(nl_list,k_list)])/ntot
k=1./k
dk=npy.sum([((k**2.0)/(ki**2.0)*ni/ntot)**2.0*(dki**2.0) for ni,ki,dki in zip(nl_list,k_list,deltak_list)])
dk=npy.sqrt(dk)
outfile.write(' '.join(['agregated',str(round(k,1)),str(round(dk,1))])+'\n')
outfile.close()

#Close the generated plots
for i in range(10):plt.close()

