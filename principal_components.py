#------------------------------------------------------------------------------------------------
#This file is part of the ost_pymodules project (https://github.com/njohner/ost_pymodules).
#
#Copyright 2015 Niklaus Johner
#
#ost_pymodules is free software: you can redistribute it and/or modify
#it under the terms of the GNU Lesser General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#ost_pymodules is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Lesser General Public License for more details.
#
#You should have received a copy of the GNU Lesser General Public License
#along with ost_pymodules.  If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------------------------
"""
.. codeauthor: Niklaus Johner <niklaus.johner@a3.epfl.ch>

This module contains functions for principal components analysis. 
It can be used to calculate principal components from a trajectory,
graphically display the principal components, reconstruct a trajectory from
a chosen set of components and so on.
"""
from ost import *
import time
import numpy as npy
import os,math

__all__=('RepresentPrincipalComponentOnStruccture', 'ProjectOnPrincipalComponentsAtomWise',\
         'ProjectOnPrincipalComponent','ReconstructTrajFromPrincipalComponents',\
        'CalculatePrincipalComponents')

def _import_numpy():
  try:
    import numpy as npy
  except: 
    raise ImportError("This function uses numpy, which could not be imported")
    return False
  return True


def CalculatePrincipalComponents(t,pc_sele,superposition_sele=None,first=0,last=-1,stride=1):
  """
  This function calculates the principal components for the positions
  of the N atoms in ev, i.e the 3N eigenvalues and eigenvectors.
  Specifically it performs an svd decomposition A=U*S*V
  Return: a tuple (U,S,V,mean_atom_pos,atom_pos_list) containing:
    -The unitary matrix U containing the eigenvectors in its columns
    -The singular values S, so that S*S are the eigenvalues
    -The unitary matrix V
    -A list of the average positions for each atoms in pc_sele
    -A list with all the positions for each atom and each frame.
  """
  if not _import_numpy:return False
  if last==-1:last=t.GetFrameCount()
  nframes=(last-first)/stride
  if superposition_sele!=None:t2=mol.alg.SuperposeFrames(t,t.GetEntity().Select(superposition_sele))
  else:t2=t
  ev=t2.GetEntity().Select(pc_sele)
  natoms=ev.GetAtomCount()
  atom_pos_list=npy.zeros([3*natoms,nframes])
  mean_atom_pos=npy.zeros(3*natoms)
  for i,a in enumerate(ev.atoms):
    vl=mol.alg.AnalyzeAtomPos(t2,a.handle,stride)[first:last]
    for j in range(3):
      vj=[v[j] for v in vl]
      vjm=npy.mean(vj)
      mean_atom_pos[3*i+j]=vjm
      atom_pos_list[3*i+j]=vj-vjm
  (U,S,VT)=npy.linalg.svd(atom_pos_list,full_matrices=False)
  return (U,S,VT.T,mean_atom_pos,atom_pos_list)

def ReconstructTrajFromPrincipalComponents(ev,U,S,V,mean_atom_pos,pc_indices=[0,1]):
  natoms=ev.GetAtomCount()
  nframes=V.shape[0]
  atom_pos_list=npy.dot(U[:,pc_indices],npy.dot(npy.diag(S)[npy.ix_(pc_indices,pc_indices)],V[:,pc_indices].T))
  atom_pos_list=npy.array([atom_pos_list[i]+mean_atom_pos[i] for i in range(3*natoms)])
  eh=mol.CreateEntityFromView(ev,1)
  t=mol.CreateCoordGroup(eh.atoms)
  for j in range(nframes):t.AddFrame(geom.Vec3List([geom.Vec3(atom_pos_list[3*i,j],atom_pos_list[3*i+1,j],atom_pos_list[3*i+2,j]) for i in range(natoms)]))
  return t

def ProjectOnPrincipalComponent(atom_pos_list,U,pc_index=0,first=0,last=-1):
  if last==-1:last=atom_pos_list.shape[1]
  return npy.dot(U[:,pc_index],atom_pos_list[:,first:last])

def ProjectOnPrincipalComponentsAtomWise(t,ev,first=0,last=-1,stride=1):
  """
  Calculates the principal components for each atom and projects the trajectory
  on them.
  It returns:
    - The list of principal axes in which the ith element is a 3x3 matrix
    with in its columns the principal axes for the ith atom in ev.
    - The positions projected on the principal axes
  """
  if last==-1:last=t.GetFrameCount()
  pvl=[]
  pcl=geom.Vec3List()
  for a in ev.atoms:
    vl=mol.alg.AnalyzeAtomPos(t,a.handle,stride)[first:last]
    pc=vl.principal_axes
    c=vl.center
    vl=geom.Vec3List([el-c for el in vl])
    pc1=pc.GetRow(0)
    pc2=pc.GetRow(1)
    pc3=pc.GetRow(2)
    pvl.append(geom.Vec3List([geom.Vec3(geom.Dot(pc1,v),geom.Dot(pc1,v),geom.Dot(pc1,v)) for v in vl]))
    pcl.append(pc)
  return (pcl,pvl)

def RepresentPrincipalComponentOnStruccture(ev,U,pc_index=0,go_name='pc0',color=gfx.RED):
  go=gfx.PrimList(go_name)
  for i in range(ev.GetAtomCount()):
    d=geom.Vec3(U[3*i,pc_index],U[3*i+1,pc_index],U[3*i+2,pc_index])
    p=ev.atoms[i].pos
    go.AddLine(p,p+50*d,color)
  return go
  
