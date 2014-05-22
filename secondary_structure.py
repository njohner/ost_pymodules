"""
Module written by Niklaus Johner (niklaus.johner@a3.epfl.ch) 01.2013
This module contains functions to analyze the secondary structure of proteins
It uses DSSP, which has to be installed and in the path for OpenStructure
DSSP does not allow CHARMM formated PDBs and therefore all chains in the
EntityView of interest need to have a single letter name.
"""
import ost
import ost.gfx as _gfx
import numpy as npy
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys

if not hasattr(ost.bindings,'dssp'):print 'DSSP is not found'

__all__=('AssignSecondaryStructure','AnalyzeSecondaryStructure','PlotSecStructureList',\
        'SimplifySecStructure','SimplifySecStructureList','FindSecStrucElements')


def AssignSecondaryStructure(prot,nmax=50):
  """
  This function calls DSSP to assign the secondary structure to a protein.
  """
  flag=True
  i=0
  while flag and i<nmax:
    try:ost.bindings.dssp.AssignDSSP(prot)
    except IOError:
      i+=1
      continue
    except:
      print sys.exc_info()[1]
      break
    flag=False
  if flag==True:print 'could not assign secondary structure'
  return

def AnalyzeSecondaryStructure(t,prot,first=0,last=-1,stride=1):
  """
  This function calculates the secondary structure of a protein for each frame in
  a trajectory.
  Inputs:
  t     : Trajectory
  prot  : EntityView for which the secondary structure will be computed
  first : First frame to be analyzed
  last  : Last frame to be analyzed
  stride: Number of frames skipped between two analyzed frames
  Outputs: 
  ss_list: A list of lists. Each element of ss_list is a list of letters corresponding
  to the secondary structure of a residue for the different frames of the trajectory.
  """
  if last==-1:last=t.GetFrameCount()
  ss_list=[[] for r in prot.residues]
  for f in range(first,last,stride):
    t.CopyFrame(f)
    AssignSecondaryStructure(prot)
    for i,r in enumerate(prot.residues):
      ss_list[i].append(str(r.GetSecStructure()))
  return ss_list

def SecStrucColorMap(ss_list,color_dict={}):
  """
  This function takes a list of lists of secondary structures and returns
  a matrix with corresponding colors that can be used to plot the secondary
  structure, for example by using matplotlib.imshow()
  Inputs:
  ss_list: a list of lists. Each element of ss_list is a list of letters corresponding
  to the secondary structure of a residue for different frames of a trajectory.
  Such a list can be obtained from the AnalyzeSecondaryStructure function
  color_dict: a dicitonary mapping a color given in rgb to each secondary structure type 
  """
  if not color_dict:
    color_dict={'H':_gfx.MAGENTA.rgb,'E':_gfx.GREEN.rgb,'T':_gfx.CYAN.rgb,'B':_gfx.DARKGREEN.rgb,'G':_gfx.ORANGE.rgb,'I':_gfx.RED.rgb,'C':_gfx.WHITE.rgb,'S':(0.8, 0.8, 0.8)}
  color_map=npy.zeros([len(ss_list),len(ss_list[0]),3])
  for j,ssl in enumerate(ss_list):
    for i,el in enumerate(ssl):
      color_map[j,i]=color_dict[el]
  return color_map,color_dict
  
def PlotSecStructureList(ss_list,color_dict={},title='',labels=[],l=0.08,b=0.08,r=0.95,t=0.95):
  """
  This function takes a list of lists of secondary structures and plots it.
  Inputs:
  ss_list: A list of lists of secondary structures as obtained from AnalyzeSecondaryStructure, see SecStrucColorMap.
  color_dict: a dicitonary mapping a color given in rgb to each secondary structure type 
  title: The title of the plot
  labels:the labels shown in the colorbar
  l,b,r,t:the left, bottom,top and right relative positions for the axes in the plot
  """
  color_map,color_dict=SecStrucColorMap(ss_list,color_dict)
  f=plt.figure()
  w=r-l
  h=t-b
  ax1=f.add_axes([l,b,0.9*w,h])
  ax2=f.add_axes([l+0.93*w,b,0.05*w,h])
  ax1.imshow(color_map,interpolation='nearest',aspect='auto')
  if not labels:labels=['H','I','G','B','E','T','S','C']
  colors=[color_dict[key] for key in labels]
  cmap=mpl.colors.ListedColormap(colors,'SecStruc')
  ncolors=len(colors)
  bounds = range(ncolors+1)
  ticks=[(i+0.5)/float(ncolors) for i in range(ncolors)]
  cb2 = mpl.colorbar.ColorbarBase(ax2,cmap=cmap,ticks=ticks,orientation='vertical')
  cb2.set_ticklabels(labels)
  f.suptitle(title)
  return f,color_map

def SimplifySecStructure(ssl):
  """
  This function takes in a list ssl of standard dssp secondary structure letters (G,H,I,E,B,T,C,S)
  and returns a simplified list where each element is assigned to the broader category of helix, beta or coil (H,E,C).
  So it basically maps G,H and I to H; E and B to E; T,C and S to C.
  """
  ss_dict={}
  for el in ['G','H','I']:ss_dict[el]='H'
  for el in ['E','B']:ss_dict[el]='E'
  for el in ['T','C','S']:ss_dict[el]='C'  
  return [ss_dict[ss] for ss in ssl]

def SimplifySecStructureList(ss_list):
  """
  This function takes a list of lists of secondary structures as obtained from AnalyzeSecondaryStructure
  and reassigns each element to the broader categories of Helix, Extended and Coil (see SimplifySecStructure)
  Inputs:
  ss_list: A list of lists of secondary structures as obtained from AnalyzeSecondaryStructure, see SecStrucColorMap.
  Outpu:
  ss_list: Simplified list of lists of secondary structures
  """
  return [SimplifySecStructure(ssl) for ssl in ss_list]

def FindSecStrucElements(prot):
  """
  This function assigns the secondary structure to the EntityView prot
  and then searches for structural elements and returns a list of index pairs
  corresponding to the begining and end of the secondary structure elements.
  Inputs:
  prot: EntityView for which the secondary structure should be analyzed
  Output: A list of tuples. Each tuple corresponds to a structural element and
  has three elements: the start index, the end index and a letter indicating the
  secondary structure of that structural element. For example: [(1,5,'H'),(8,15,'E')]
  means that there is a helix spanning prot.residues[1:6] and a beta strand prot.residues[8:16].
  """
  AssignSecondaryStructure(prot)
  ssl=[str(r.GetSecStructure()) for r in prot.residues]
  ssl=SimplifySecStructure(ssl)
  ss_ele_list=[]
  current='C'
  for i,el in enumerate(ssl):
    if el!=current:
      if current=='C':
        start=i
        current=el
      else:
        stop=i-1
        ss_ele_list.append((start,stop,current))
        current=el
        if el!='C':start=i
  return ss_ele_list
