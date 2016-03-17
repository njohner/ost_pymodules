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
.. codeauthor:: Niklaus Johner <niklaus.johner@a3.epfl.ch>

This module is used to generate densities from trajectories and extract
surfaces from these densities in the form of sets of points.
Some algebra can be done on these surfaces, found in the surface_alg module
"""
try:
  import time
  import math
except: print 'module needs time and math module but could not import them'

from ost import *

__all__=('CreateDensityMapFromEntityView','CreateDensityMapFromTrajectory','CreateDensityMapFromTrajectoryWithPBC',\
        'CalculateVolumeFromDensityMap','GetBoundaryBetweenDensities','GetMaxDensitySurface')

def _DetermineMaxBounds(sele_view,t,update_sele=False,update_sele_chain=''):
  if update_sele:bounds=sele_view.Select(update_sele_chain).bounds
  else:bounds=sele_view.bounds
  vmin=bounds.min
  vmax=bounds.max
  for f in range(t.GetFrameCount()):
    t.CopyFrame(f)
    if update_sele:bounds=sele_view.Select(update_sele_chain).bounds
    else:bounds=sele_view.bounds
    for i in range(3):
      if bounds.min[i]<vmin[i]:vmin[i]=bounds.min[i]
      if bounds.max[i]>vmax[i]:vmax[i]=bounds.max[i]
  return (vmin,vmax)

def _VecToNeighborCells(basis_vec):
  v1=geom.Vec3(basis_vec.x,0,0)
  v2=geom.Vec3(0,basis_vec.y,0)
  v3=geom.Vec3(0,0,basis_vec.z)
  vec_list=[v1,v2,v3,v1+v2,v1+v3,v2+v3,v1-v2,v1-v3,v2-v3]
  vec_list.extend([v1+v2+v3,v1+v2-v3,v1-v2+v3,-v1+v2+v3])
  vec_list.extend([-el for el in vec_list])
  return vec_list

def _PositiveVecToNeighborCells(basis_vec):
  v1=geom.Vec3(basis_vec.x,0,0)
  v2=geom.Vec3(0,basis_vec.y,0)
  v3=geom.Vec3(0,0,basis_vec.z)
  vec_list=[v1,v2,v3,v1+v2,v1+v3,v2+v3,v1-v2,v1-v3,v2-v3]
  vec_list.extend([v1+v2+v3,v1+v2-v3,v1-v2+v3,-v1+v2+v3])
  return vec_list

def CreateDensityMapFromEntityView(view,sampling=1,resolution=5,margin=0,ucell_center=None,ucell_vecs=None,cut_back=True):
  basis_vec=view.bounds.size
  if not ucell_center:ucell_center=view.bounds.center
  if margin!=0:
    try: import pbc_utilities
    except:'needs the pbc_utilities module to use PBC'
    if not ucell_vecs:
      print 'unit cell vectors needed when margin>0 to extend the entity'
      return None
    view=pbc_utilities.ExtendEntityWithPBC(view,ucell_center,ucell_vecs,margin)
  basis_vec_ext=view.bounds.size
  box_center=view.bounds.center
  origin=box_center-basis_vec_ext/2.
  real_size=basis_vec_ext
  #we prepare the map
  map_size=img.Size(int(real_size.x)*sampling, int(real_size.y)*sampling, int(real_size.z)*sampling)
  den_map=img.CreateImage(map_size)
  den_map.SetAbsoluteOrigin(origin)
  den_map.SetSpatialSampling(geom.Vec3(1./sampling,1./sampling,1./sampling))
  mol.alg.EntityToDensityRosetta(view.Select(''), den_map, mol.alg.HIGH_RESOLUTION, resolution)
  #Now we cut the map back to the size of the entity
  if not cut_back:return den_map
  ext=img.Extent()
  ext.SetEnd(img.Point(den_map.CoordToIndex(box_center+basis_vec/2.)))
  ext.SetStart(img.Point(den_map.CoordToIndex(box_center-basis_vec/2.)))
  den_map=den_map.Extract(ext)
  den_map.SetAbsoluteOrigin(origin)
  return den_map
    
def CreateDensityMapFromTrajectory(traj, sele_view, sampling=0.5, stride=1,resolution=2.0,update_sele=False,update_sele_chain='',margin=10):
  """
  This function creates a density map for the position of the atoms of an EntityView during a trajectory
  It returns a map
  If update_sele is set to true, then for each frame the density is calculated from sele_view.Select(update_sele_chain)
  Sampling was changed to 1/sampling, careful!
  """
  #For speed we filter the trajectory
  t=traj.Filter(sele_view,stride=stride)
  sele_view=t.GetEntity().CreateFullView()
  #We setup the map boundaries and create the map object
  (vmin,vmax)=_DetermineMaxBounds(sele_view,t,update_sele,update_sele_chain)
  real_size=vmax-vmin+geom.Vec3(margin,margin,margin)
  print 'map size',real_size
  map_size=img.Size(int(real_size.x*sampling), int(real_size.y*sampling), int(real_size.z*sampling))
  den_map=img.CreateImage(map_size)
  den_map.SetSpatialSampling(geom.Vec3(1./sampling,1./sampling,1./sampling))
  origin=vmin-geom.Vec3(0.5*margin,0.5*margin,0.5*margin)
  #origin=geom.Vec3(int(origin.x),int(origin.y),int(origin.z))
  den_map.SetAbsoluteOrigin(origin)

  #Now we go through all frames and for each frame we extract the density and add it up
  N=0
  sv=sele_view
  for i in range(0,t.GetFrameCount()):
    t.CopyFrame(i)
    if update_sele:sv=sele_view.Select(update_sele_chain)
    mol.alg.EntityToDensityRosetta(sv, den_map, mol.alg.HIGH_RESOLUTION, resolution)
    N+=1
  for p in den_map:
    den_map.SetReal(p,den_map.GetReal(p)/float(N))
  return den_map

def CreateDensityMapFromTrajectoryWithPBC(traj, sele_view, sampling=2, stride=1,resolution=2.0, cell_centers=None,cell_vectors=None,margin=10):
  """
  This function creates a density map for the position of the atoms of an EntityView during a trajectory
  It returns a map and uses periodic boundary conditions while generating the map
  """
  try:import pbc_utilities,trajectory_utilities
  except:print 'could not load pbc_utilities and trajectory_utilities'
  if not cell_centers:
    print 'Extracting cell centers from the center of mass of the sele_view over the trajectory'
    cell_centers=mol.alg.AnalyzeCenterOfMassPos(traj,sele_view)
  if not cell_vectors:
    print 'Extracting cell vectors from the unit cell information in the trajectory'
    cell_vectors=trajectory_utilities.GetCellVectorsList(traj)
  #For speed we filter the trajectory
  t=traj.Filter(sele_view,stride=stride)
  sele_view=t.GetEntity().CreateFullView()
  edi_ref=sele_view.handle.EditXCS()
  T=geom.Mat4()
  edi_ref.SetTransform(T)
  #We setup the map boundaries and create the map object
  (vmin,vmax)=_DetermineMaxBounds(sele_view,t)
  m=geom.Vec3(margin+5,margin+5,margin+5)
  (vmin,vmax)=(vmin-m,vmax+m)
  real_size=vmax-vmin#+geom.Vec3(margin+10,margin+10,margin+10)
  map_size=img.Size(int(real_size.x/sampling), int(real_size.y/sampling), int(real_size.z/sampling))
  den_map=img.CreateImage(map_size)
  den_map.SetSpatialSampling(sampling)
  origin=vmin#-geom.Vec3(5+margin/2., 5+margin/2., 5+margin/2.)
  #origin=geom.Vec3(int(origin.x),int(origin.y),int(origin.z))
  den_map.SetAbsoluteOrigin(origin)

  #Now we go through all frames and for each frame we extract the density and add it up
  N=0
  for i in range(0,t.GetFrameCount()):
    T=geom.Mat4()
    edi_ref.SetTransform(T)
    t.CopyFrame(i)
    extended_eh=pbc_utilities.ExtendEntityWithPBC(sele_view,cell_centers[i],cell_vectors[i],margin)
    mol.alg.EntityToDensityRosetta(extended_eh.Select(''), den_map, mol.alg.HIGH_RESOLUTION, resolution)
    N+=1
  for p in den_map:
    den_map.SetReal(p,den_map.GetReal(p)/float(N))
  return den_map



def CalculateVolumeFromDensityMap(den_map,cutoff=0.1):
  """
  This function calculates the volume of a density map with a certain cutoff
  It simply sums up all the volume where the density is higher than the cutoff
  """
  c=0
  for el in den_map:
    if den_map[el]>=cutoff:c+=1
  s=den_map.GetSpatialSampling()
  v=s[0]*s[1]*s[2]
  return c*v


def GetBoundaryBetweenDensities(den_map1,den_map2,cutoff=0.05,PBC=False,cell_center=None,cell_size=None,sampling=None,density_ratio=1.0):
  """
  This function returns a set of points at the the boundary between
  Two density maps. The boundary is determined as the set of points where
  the ratio between the densities from the two maps changes from >1 to <1. 
  """
  if PBC:
    if not cell_center:
      print 'Entity bounds center used as cell_center'
      cell_center=eh.bounds.center
    if not cell_size:
      print 'Entity bounds size used as cell_size'
      cell_size=eh.bounds.size
    half_cell=cell_size/2.
  boundary=geom.Vec3List()
  origin=den_map1.GetAbsoluteOrigin()
  #if not sampling:step=den_map1.GetSpatialSampling()
  #else:step=sampling
  if not sampling:step=1
  else:step=sampling
  size=[int(math.ceil(den_map1.size[0]/step)),int(math.ceil(den_map1.size[1]/step)),int(math.ceil(den_map1.size[2]/step))]
  #size=den_map1.size
  t1=time.time()
  count=0
  count2=0
  n_tot=(size[0])*(size[1])*(size[2])
  #We add a condition to avoid having a gap larger than sampling of the map at the periodic boundaries
  if PBC:
    v=(cell_center-half_cell)-origin
    v1=geom.Vec3(int(v.x)+step,int(v.y)+step,int(v.z)+1)-v
    v=(cell_center+half_cell)-origin
    v2=v-geom.Vec3(int(v.x),int(v.y),int(v.z))
    v=v1+v2
    v=geom.Vec3(int(v.x),int(v.y),int(v.z))
    print v
  for (i1,j1,k1) in [(i,j,k) for i in range(size[0]) for j in range(size[1]) for k in range(size[2])]:
    i=step*i1
    j=step*j1
    k=step*k1
    count+=1
    count2+=1
    if count==1000000:
      count=0
      print count2,'out of',n_tot,'points assessed in',time.time()-t1,'seconds'
    c=den_map1.IndexToCoord(img.Point(i,j,k))
    if PBC:
      dc=c-cell_center
      if any([abs(di)-hci-vi>0.0 for di,hci,vi in zip(dc.data,half_cell.data,v.data)]):continue
      if any([abs(di)-hci>0.0 and di<0.0 for di,hci in zip(dc.data,half_cell.data)]):continue
    d1=den_map1.GetReal(img.Point(den_map1.CoordToIndex(c)))
    if d1<cutoff:continue
    d2=den_map2.GetReal(img.Point(den_map2.CoordToIndex(c)))
    if (not density_ratio*d1>d2) or (d2<cutoff):continue
    for (l,m,n) in [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]:
      if PBC:p=geom.WrapVec3(den_map1.IndexToCoord(img.Point(i+l,j+m,k+n)),cell_center,cell_size)
      else:p=den_map1.IndexToCoord(img.Point(i+l,j+m,k+n))
      d1=den_map1.GetReal(img.Point(den_map1.CoordToIndex(p)))
      d2=den_map2.GetReal(img.Point(den_map2.CoordToIndex(p)))
      if d2>density_ratio*d1:
        boundary.append(c)
        break
  return boundary



def GetMaxDensitySurface(den_map,cutoff=0.1,PBC=False,cell_center=None,cell_size=None,sampling=None):
  """
  This function calculates a set of points that have local maximal density 
  """
  if PBC:
    if not cell_center:
      print 'Entity bounds center used as cell_center'
      cell_center=eh.bounds.center
    if not cell_size:
      print 'Entity bounds size used as cell_size'
      cell_size=eh.bounds.size
    half_cell=cell_size/2.
  boundary=geom.Vec3List()
  origin=den_map.GetAbsoluteOrigin()
  #if not sampling:step=den_map.GetSpatialSampling()
  #else:step=sampling
  if not sampling:step=1
  else:step=sampling
  size=[int(math.ceil(den_map.size[0]/step)),int(math.ceil(den_map.size[1]/step)),int(math.ceil(den_map.size[2]/step))]
  t1=time.time()
  count=0
  count2=0
  n_tot=(size[0])*(size[1])*(size[2])
  vec_list=_PositiveVecToNeighborCells(geom.Vec3(1,1,1))
  #We add a condition to avoid having a gap larger than sampling of the map at the periodic boundaries
  if PBC:
    v=(cell_center-half_cell)-origin
    v1=geom.Vec3(int(v.x)+step,int(v.y)+step,int(v.z)+step)-v
    v=(cell_center+half_cell)-origin
    v2=v-geom.Vec3(int(v.x),int(v.y),int(v.z))
    v=v1+v2
    v=geom.Vec3(int(v.x),int(v.y),int(v.z))
    print v
  for (i1,j1,k1) in [(i,j,k) for i in range(size[0]) for j in range(size[1]) for k in range(size[2])]:
    i=step*i1
    j=step*j1
    k=step*k1
    count+=1
    count2+=1
    if count==1000000:
      count=0
      print count2,'out of',n_tot,'points assessed in',time.time()-t1,'seconds'
    p1=den_map.IndexToCoord(img.Point(i,j,k))
    if PBC:
      dc=p1-cell_center
      if any([abs(di)-hci-vi>0.0 for di,hci,vi in zip(dc.data,half_cell.data,v.data)]):continue
      if any([abs(di)-hci>0.0 and di<0.0 for di,hci in zip(dc.data,half_cell.data)]):continue
    d1=den_map.GetReal(img.Point(den_map.CoordToIndex(p1)))
    if d1<cutoff:continue
    #for [(l1,m1,n1),(l2,m2,n2)] in [[(1,0,0),(-1,0,0)],[(0,1,0),(0,-1,0)],[(0,0,1),(0,0,-1)]]:
    for v1 in vec_list:
      if PBC:
        p2=geom.WrapVec3(den_map.IndexToCoord(img.Point(geom.Vec3(i,j,k)+v1)),cell_center,cell_size)
        p3=geom.WrapVec3(den_map.IndexToCoord(img.Point(geom.Vec3(i,j,k)-v1)),cell_center,cell_size)
      else:
        p2=den_map.IndexToCoord(img.Point(geom.Vec3(i,j,k)+v1))
        p3=den_map.IndexToCoord(img.Point(geom.Vec3(i,j,k)-v1))
      d2=den_map.GetReal(img.Point(den_map.CoordToIndex(p2)))
      if d2>d1:continue
      
      d3=den_map.GetReal(img.Point(den_map.CoordToIndex(p3)))
      if d3>d1:continue
      boundary.append(p1)
      break
  return boundary

            
            