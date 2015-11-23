"""
.. codeauthor:: Niklaus Johner <niklaus.johner@a3.epfl.ch>

This module contains functions to determine lipid tilt and splay angles and
Calculate the elastic properties of the membrane from there.

If you use this code please cite ref. [1]_, where the method and its implementation 
is described in details or ref. [2]_, where the method is applied to complex lipidic phases.

To ensure proper treatment of the periodic boundary conditions, the trajectory
should first be extended to neighboring unit cells, and then aligned.
The tilts and splays are then calculated only for the central cell, the others 
being only used to ascertain correct treatment of the PBCs.

References
-------------
.. [1] Niklaus Johner, D. Harries and G. Khelashvili, 
       "Release note of the lipid tilt and splay code"
.. [2] Niklaus Johner, D. Harries and G. Khelashvili,
       "Curvature and lipid packing modulate the elastic properties of lipid assemblies: comparing the HII and lamellar phases."
       The Journal of Physical Chemistry Letters 5, no. 23 (2014), 4201-6.
"""
try:
  from ost import *
  import time
  import numpy as npy
  import os,math
  from scipy.optimize import curve_fit
  import matplotlib as mpl
  try:gui.dng
  except:mpl.use('Agg')
  import matplotlib.pyplot as plt
  import entity_alg,trajectory_utilities,surface_alg,file_utilities
except:
  print 'could not import at least one of the modules nedded: ost, time, numpy, os, math, entity_alg,trajectory_utilities,surface_alg,file_utilities'


__all__=('GetBoundaryBetweenViews','AssignNormalsToLipids','AnalyzeLipidTilts',\
        'AnalyzeLipidSplays','AnalyzeLipidTiltAndSplay',\
        'WriteTiltDict','WriteSplayDict',\
        'FitTiltDistribution','FitSplayDistribution','AnalyzeAreaPerLipid')


def _CalculateSplayAngle(v11,v12,v21,v22,v1p,v2p,n1,n2,distance_cutoff):
  """
  Calculates the splay angle for a pair of lipids.

  :param v11: headgroup of the first lipid
  :param v12: terminal tail atoms of the first lipid
  :param v21: headgroup of the second lipid
  :param v22: terminal tail atoms of the second lipid
  :param v1p: Atoms of the first lipid situated at the neutral/pivotal plane
  :param v1p: Atoms of the second lipid situated at the neutral/pivotal plane
  :param distance_cutoff: The maximal distance between lipids considered for splay calculation.

  :type v11: :class:`~ost.mol.EntityView`
  :type v12: :class:`~ost.mol.EntityView`
  :type v21: :class:`~ost.mol.EntityView`
  :type v22: :class:`~ost.mol.EntityView`
  :type v1p: :class:`~ost.mol.EntityView`
  :type v2p: :class:`~ost.mol.EntityView`
  :type distance_cutoff: :class: `float`

  :return: returns a tuple of floats composed of the splay angle and the distance between the lipids
  :rtype: (:class:`float`,:class:`float`)
  """
  if geom.Dot(n1,n2)<0:return None
  x=v2p.GetCenterOfMass()-v1p.GetCenterOfMass()
  d=geom.Length(x)
  x=x-geom.Dot(x,n1)*n1
  #d=geom.Length(x)
  if d>distance_cutoff:return None
  x=geom.Normalize(x)
  v1=geom.Normalize(v11.GetCenterOfMass()-v12.GetCenterOfMass())
  v2=geom.Normalize(v21.GetCenterOfMass()-v22.GetCenterOfMass())
  #x=geom.Normalize(v2p.GetCenterOfMass()-v1p.GetCenterOfMass())
  #d=geom.Length(v2p.GetCenterOfMass()-v1p.GetCenterOfMass())
  return ((geom.Dot(v2,x)-geom.Dot(v1,x)-geom.Dot(n2,x)+geom.Dot(n1,x))/d,d)
  

CalculateInterfaceFromTraj=trajectory_utilities.CalculateInterfaceFromTraj

def _AssignNormalsFromSurfaceToResidues(t,sele,surface,within_size=10):
  """
  This function assigns a normal vector to each residue for each frame of a trajectory
  from the closest point on the surface.
  
  :param t: the trajectory
  :param sele: the selection to which normals will be assigned.
  :param surface: the surface. Each atom of the surface should have an associated
                  normal vector as Vec3 property "n".
  :param within_size: Size of surrounding used to find the closest atom on the surface.
                      This parameter only optimizes the speed of the calculation

  :type t: :class:`~ost.mol.CoordGroupHandle`
  :type sele: :class:`~ost.mol.EntityView`
  :type surface: :class:`~ost.mol.EntityView`
  :type within_size: :class:`float`

  :return: A list containing one **Vec3List** for each residue in **sele**.
  """
  surface=surface.Select('')#Make sure we don't have to make that selection each time
  ln=[geom.Vec3List() for r in sele.residues]
  res_view_list=[r.Select('') for r in sele.residues]
  for i in range(t.GetFrameCount()):
    t.CopyFrame(i)
    for j,r in enumerate(res_view_list):
      within=surface.FindWithin(r.center_of_atoms,within_size)
      if not len(within)==0:a=entity_alg.FindClosestAtom(r,mol.CreateViewFromAtoms(within))
      else:a=entity_alg.FindClosestAtom(r,surface)
      ln[j].append(a.GetVec3Prop('n'))
  return ln

def _CalculateTilts(t,lipids,normals,head_sele,tail_sele,prot_cm=None,bool_prop=''):
  """
  This function calculates the lipid tilts.
  :param t: the trajectory
  :param lipids: the selection for which tilts will be calculated.
  :param head_sele: the selection used for the headgroups.
  :param tail_sele: the selection used for the tails. 
  :param prot_cm: A point from which the distance will be calculated for each tilt.
                  This is typically the position of the center of mass of a protein.

  :type t: :class:`~ost.mol.CoordGroupHandle`
  :type lipids: :class:`~ost.mol.EntityView`
  :type head_sele: :class:`str`
  :type tail_sele: :class:`str`
  :type bool_prop: :class:`str`
  :type prot_cm: :class:`~ost.mol.EntityView`
  
  :return: A tuple of arrays **(tilts, prot_dist)**. Each array has the shape N\ :subscript:`Lipids`\ x N\ :subscript:`Frames`
   If **prot_cm** is **None**, the second list is empty.
  :rtype: (:class:`npy.array`,:class:`npy.array`)
  """
  tilts=[]
  prot_dist=[]
  for i,r in enumerate(lipids.residues):
    if bool_prop and not r.GetBoolProp(bool_prop):continue
    p1=mol.alg.AnalyzeCenterOfMassPos(t,r.Select(head_sele))
    p2=mol.alg.AnalyzeCenterOfMassPos(t,r.Select(tail_sele))
    n=normals[i]
    tilts.append([geom.Angle(ni,p1i-p2i) for ni,p1i,p2i in zip(n,p1,p2)])
    if prot_cm:
      prot_dist.append([min(geom.Length(geom.Vec2(cm-el1)),geom.Length(geom.Vec2(cm-el2))) for cm,el1,el2 in zip(prot_cm,p1,p2)])
  return (npy.array(tilts),npy.array(prot_dist))

def GetBoundaryBetweenViews(t,waters,lipids,outdir='',density_cutoff=None,stride=1,within_size_normals=5.0,filename_basis=''):
  """
  This function determines the interface between two views, typically the lipid-water interface and assigns normals 
  to every point on the surface.

  :param t: the trajectory
  :param waters: First view
  :param lipids: Second view
  :param outdir: Path to output directory
  :param density_cutoff: Interface will not be calculated for regions where the density is lower than this cutoff.
  :param stride: stride used to calculate the average density from the trajectory.
  :param within_size_normals: radius of the patch used to determine the normals on the water-lipid interface.
  :param filename_basis: used as first part in the name of all the files generated.

  :type t: :class:`~ost.mol.CoordGroupHandle`
  :type waters: :class:`~ost.mol.EntityView`
  :type lipids: :class:`~ost.mol.EntityView`
  :type outdir: :class:`str`
  :type density_cutoff: :class:`float`
  :type stride: :class:`int`
  :type within_size: :class:`float`
  :type filename_basis: :class:`str`
  
  :return: A tuple **(water_filtered,lipid_filtered,b_eh)** containing the density for 
    the first and second views and an entity for the boundary between the two views.
    Every atom in the boundary has an associated normal vector set as a Vec3 property 'n'.

  WARNING: Interface was changed. Taking 2 views instead of a list of lipid names and water name.
  The order of the two views was also inversed. I also removed the PBC, cell_center and cell_size parameters

  """

  t0=time.time()
  eh=t.GetEntity()
  t2=t.Filter(eh.Select(''),stride=stride)
  print 'Calculating densities on',t2.GetFrameCount(),'frames'
  (water_filtered,lipid_filtered,boundary_filtered)=CalculateInterfaceFromTraj(t2,waters,lipids,False,None,None,10,density_cutoff,outdir,filename_basis)
  b_eh=entity_alg.CreateEntityFromVec3List(boundary_filtered)
  print 'calculating the normals',time.time()-t0
  surface_alg.CalculateNormals(b_eh,within_size_normals,False,None,None)
  surface_alg.OrientNormalsAlongDensity(b_eh,lipid_filtered)
  if outdir:
    io.SavePDB(b_eh,os.path.join(outdir,filename_basis+'boundary.pdb'))
    entity_alg.WriteVec3Prop(b_eh,"n",os.path.join(outdir,filename_basis+'boundary_normals.txt'),index=True)
  return (water_filtered,lipid_filtered,b_eh)
      
def AssignNormalsToLipids(t,eh,b_eh,lipid_names,head_group_dict):
  """
  :param t: The trajectory
  :param eh: The associated entity
  :param b_eh: The surface from which the normals will get assigned to the lipids.
  :param lipid_names: a list of the lipid residue names to which normals will be assigned.
  :param head_group_dict: a dictionary containing the selections defining the headgroup 
   for each lipid type. There should be one entry for each residue name in **lipid_names**

  :type t: :class:`~ost.mol.CoordGroupHandle`
  :type eh: :class:`~ost.mol.Entity`
  :type b_eh: :class:`~ost.mol.Entity`
  :type lipid_names: :class:`list`
  :type head_group_dict: :class:`dict`
  
  :return: A dictionary with one entry for each residue name in **lipid_names**.
   Each element in the dictionary is a :class:`list`\ (:class:`~ost.geom.Vec3List`\ ). Each element in the list
   corresponds to one residue and each :class:`~ost.geom.Vec3` in the :class:`~ost.geom.Vec3List` is the normal for one frame
   of the trajectory.
  """
  t0=time.time()
  #Now we assign a normal for each lipid in each frame
  print 'assigning normals to lipids',lipid_names
  lipid_normal_dict={}
  for lipid_name in lipid_names:
    print 'rname='+lipid_name+' and '+head_group_dict[lipid_name]
    t0=time.time()
    lipids=eh.Select('rname='+lipid_name+' and '+head_group_dict[lipid_name])
    print lipids.GetResidueCount()
    lipid_normal_dict[lipid_name]=_AssignNormalsFromSurfaceToResidues(t,lipids,b_eh)
    print 'done for',lipid_name,'in',time.time()-t0,'seconds'
  return lipid_normal_dict

def AnalyzeLipidTilts(t,eh,lipid_names,lipid_normal_dict,head_group_dict,tail_dict,prot_cm=None,bool_prop=''):
  """
  This function calculates the lipid tilts from a trajectory.

  :param t: The trajectory
  :param eh: The associated entity
  :param lipid_names: List of the residue names of the different lipids in the system
  :param lipid_normal_dict: Dictionary of normal vectors. One entry for every lipid type (element in lipid_names)
   Every entry is a :class:`list`\ (:class:`~ost.geom.Vec3List`\ ) of normals for every frame for every lipid of that type 
   (size of list:N\ :subscript:`Lipids`\ x N\ :subscript:`Frames`).
  :param head_group_dict: Dictionary containing a selection string for each lipid type
   that is used to determine the position of the lipid headgroups (center of mass of the selection).
  :param tail_dict: Dictionary containing a selection string for each lipid type
   that is used to determine the position of the lipid tails (center of mass of the selection).
  :param prot_cm: A list of position (one for each frame). If specified the each tilt the distance between this position
    and the lipid in question will also be returned. This is typically used to calculate local properties of the membrane around an insertion.
  :param bool_prop: Boolean property assigned to lipids to determine whether they should be considered in the tilt calculations.
   This is typically used to treat the periodic boundary conditions, to differentiate lipids from the central unit cell, for which tilt and 
   splay are calculated, from the lipids from neighboring unit cells, used only to ensure correct treatment of PBC.

  :type t: :class:`~ost.mol.CoordGroupHandle`
  :type eh: :class:`~ost.mol.EntityHandle`
  :type lipid_normal_dict: :class:`dict`
  :type lipid_names: :class:`str`
  :type head_group_dict: :class:`dict`
  :type tail_dict: :class:`dict`
  :type prot_cm: :class:`~ost.geom.Vec3List`
  :type bool_prop: :class:`bool`

  :return: Dictionary of lipid tilts. One entry for every lipid type (element in **lipid_names**).
           Every entry is a :class:`list` with two elements. The first one is a :class:`list`\ (:class:`~ost.geom.FloatList`\ ) of tilts 
           for every frame for every lipid of that type (size of list:N\ :subscript:`Lipids`\ x N\ :subscript:`Frames`).
           The second element contains the list of distances to **prot_cm** if it was defined and is empty otherwise.

  WARNING: Removed parameters PBC, cell_center, cell_size
  """
  lipid_tilt_dict={}
  print "analyzing lipid tilts for",lipid_names
  for ln in lipid_names:
    print "starting for",ln
    lipids=eh.Select('rname='+ln)
    lipid_tilt_dict[ln]=_CalculateTilts(t,lipids,lipid_normal_dict[ln],head_group_dict[ln],tail_dict[ln],prot_cm,bool_prop)
  return lipid_tilt_dict
  

def AnalyzeLipidSplays(t,eh,lipid_names,head_group_dict,tail_dict,lipid_normal_dict,lipid_tilt_dict,distance_sele_dict,distance_cutoff=10,bool_prop=''):
  """
  This function calculates the lipid splays from a trajectory.

  :param t: The trajectory
  :param eh: The associated entity
  :param lipid_names: A list of the residue names of the different lipids in the system
  :param head_group_dict: Dictionary containing a selection string for each lipid type
   that is used to determine the position of the lipid headgroups (center of mass of the selection).
  :param tail_dict: Dictionary containing a selection string for each lipid type
   that is used to determine the position of the lipid tails (center of mass of the selection).
  :param lipid_normal_dict: Dictionary of normal vectors. One entry for every lipid type (element in lipid_names)
                            Every entry is a :class:`list`\ (:class:`~ost.geom.Vec3List`\ ) of normals for every frame 
                            for every lipid of that type (size of list:N\ :subscript:`Lipids`\ x N\ :subscript:`Frames`).
  :param lipid_tilt_dict: Dictionary of lipid tilts. One entry for every lipid type (element in lipid_names)
                          Every entry is a :class:`list`\ (:class:`~ost.geom.FloatList`\ ) of tilts for every frame 
                          for every lipid of that type (size of list:N\ :subscript:`Lipids`\ x N\ :subscript:`Frames`).
  :param distance_sele_dict: Dictionary containing a selection string for each lipid type that is used
   to calculate the distance between lipids (center of mass distance). The center of mass of these selections should lie
   on the neutral plane.
  :param distance_cutoff: Lipid pairs further apart than this distance will not be considered for splay calculation.
  :param bool_prop: Boolean property assigned to lipids to determine whether they should be considered in the splay calculations.
   This is typically used to treat the periodic boundary conditions, to differentiate lipids from the central unit cell, for which tilt and 
   splay are calculated, from the lipids from neighboring unit cells, used only to ensure correct treatment of PBC.

  :return: Dictionary of splays, containing one entry for each possible lipid pairs
  
  WARNING: Changed the name of that function from AnalyzeLipidSplay to AnalyzeLipidSplays.
  """
  nframes=t.GetFrameCount()
  t0=time.time()
  splay_dict={}
  mean_dist_dict={}
  area_dict={}
  n_dict={}
  
  lipid_sele_string='rname='+','.join([rname for rname in distance_sele_dict.keys()])
  distance_sele_string=' or '.join(['(rname={0} and {1})'.format(rname,distance_sele_dict[rname]) for rname in distance_sele_dict.keys()])
  print 'selecting all the lipids using',lipid_sele_string
  print 'selection for distance calculation',distance_sele_string
  lipids=eh.Select(lipid_sele_string)
  distance_sele=lipids.Select(distance_sele_string)
  print bool_prop,lipid_names
  head_group_view_list_dict={}
  tail_view_list_dict={}
  distance_view_list_dict={}
  for ln in lipid_names:
    v=lipids.Select('rname='+ln)
    head_group_view_list_dict[ln]=[r.Select(head_group_dict[ln]) for r in v.residues]
    tail_view_list_dict[ln]=[r.Select(tail_dict[ln]) for r in v.residues]
    distance_view_list_dict[ln]=[r.Select(distance_sele_dict[ln]) for r in v.residues]
  splay_dict={}#We prepare the dictionary that will contain the splays for each type of lipid pair
  for ln1 in lipid_names:
    for ln2 in lipid_names:splay_dict[ln1+'-'+ln2]=[]
  for f in range(t.GetFrameCount()):
    if f%50==0:print "frame {0} in {1} seconds".format(f,time.time()-t0)
    t.CopyFrame(f)
    for r1 in lipids.residues:
      if not r1.GetBoolProp(bool_prop):continue
      ln1=r1.GetName()
      i=r1.GetIntProp('index')
      within=distance_sele.Select("{0}<>[cname={1} and rnum={2}]".format(distance_cutoff,r1.chain.name,r1.number.num))
      within=within.Select("not (cname={0} and rnum={1})".format(r1.chain.name,r1.number.num))
      #if within.GetResidueCount()<2:continue
      v11=head_group_view_list_dict[ln1][i]
      v12=tail_view_list_dict[ln1][i]
      v13=distance_view_list_dict[ln1][i]
      n1=lipid_normal_dict[ln1][i][f]
      a1=lipid_tilt_dict[ln1][0][i][f]
      for r2 in within.residues:
        ln2=r2.name
        j=r2.GetIntProp('index')
        v21=head_group_view_list_dict[ln2][j]
        v22=tail_view_list_dict[ln2][j]
        v23=distance_view_list_dict[ln2][j]
        n2=lipid_normal_dict[ln2][j][f]
        a2=lipid_tilt_dict[ln2][0][j][f]
        if a1>0.8 or a2>0.8:continue
        s=_CalculateSplayAngle(v11,v12,v21,v22,v13,v23,n1,n2,distance_cutoff)
        if s:splay_dict[ln1+'-'+ln2].append(s)
  #We finalize the dictionary by merging equivalent entries
  nsplays=0
  splay_dict_def={}
  for ln1 in lipid_names:
    for ln2 in lipid_names:
      key=ln1+'-'+ln2
      key2=ln2+'-'+ln1
      nsplays+=len(splay_dict[key])
      if key2 in splay_dict_def.keys():
        splay_dict_def[key2].extend(splay_dict[key])
      else:splay_dict_def[key]=splay_dict[key]
  print 'done in {0} seconds. Computed {1} splays'.format(time.time()-t0,nsplays)
  return splay_dict_def

def AnalyzeLipidTiltAndSplay(t,lipid_names,head_group_dict,tail_dict,distance_cutoff=10.0,within_size_normals=10.0
                      ,distance_sele_dict={},water_name='TIP3',outdir='',density_cutoff=None
                      ,prot_sele=None,density_stride=10,tilt_bool_prop='',splay_bool_prop='',filename_basis='',sele_dict={}):
  """
  This function is a wrapper to determine the membrane elastic moduli from the lipid titls and splays.
  Periodic boundary conditions are not treated explicitely here and should be treated as suggested in the
  description of this module.

  :param t: The trajectory
  :param lipid_names: List of the residue names of the different lipids in the system
  :param head_group_dict: Dictionary containing a selection string for each lipid type
   that is used to determine the position of the lipid headgroups (center of mass of the selection).
  :param tail_dict: Dictionary containing a selection string for each lipid type
   that is used to determine the position of the lipid tails (center of mass of the selection).
  :param distance_cutoff: Lipid pairs further apart than this distance will not be considered for splay calculation.
  :param within_size_normals: radius of the patch used to determine the normals on the water-lipid interface.
  :param distance_sele_dict: Dictionary containing a selection string for each lipid type that is used
   to calculate the distance between lipids (center of mass distance). The center of mass of these selections should lie
   on the neutral plane.
  :param water_name: Residue name of the waters (used to calculate the water-lipid interface).
  :param outdir: Path to output directory. If none, no files will be written.
  :param density_cutoff: Interface will not be calculated for regions where the density is lower than this cutoff.
  :param prot_sele: Selection string used to determine the position of the protein. This is used to calculate the distance
   between lipids and the protein which will be returned together with splays and tilts and can be used to determine the vatiation
   in membrane properties around a protein
  :param density_stride: Stride to be used for the calculation of the average lipid and water densities used to determine the interface.
   Using every single frame can slow down the calculation.
  :param tilt_bool_prop: Boolean property assigned to lipids to determine whether they should be considered in the titl calculations.
   This is typically used to treat the periodic boundary conditions, to differentiate lipids from the central unit cell, for which tilt and 
   splay are calculated, from the lipids from neighboring unit cells, used only to ensure correct treatment of PBC.
  :param splay_bool_prop: Same as **tilt_bool_prop** but for the calculation of splays
  :param filename_basis: used as first part in the name of all the files generated.
  :param sele_dict: Dictionary containing selection strings used to separate the system in several parts for the calculation
   This can be used for example to make the tilt and splay calculations separately for each leaflet of a bilayer.
   In such a case **sele_dict** would be something like **sele_dict={"upper":"z>0","lower":"z<0"}**

  :type t: :class:`~ost.mol.CoordGroupHandle`
  :type lipid_names: :class:`str`
  :type head_group_dict: :class:`dict`
  :type tail_dict: :class:`dict`
  :type distance_cutoff: :class:`float`
  :type within_size_normals: :class:`float`
  :type distance_sele_dict: :class:`dict`
  :type outdir: :class:`str`
  :type density_cutoff: :class:`float`
  :type prot_sele: :class:`str`
  :type density_stride: :class:`int`
  :type tilt_bool_prop: :class:`bool`
  :type splay_bool_prop: :class:`bool`
  :type filename_basis: :class:`str`
  :type sele_dict: :class:`dict`

  :return: A tuple **(lipid_tilt_dict,lipid_normal_dict,splay_dict,b_eh)**, where **lipid_tilt_dict**,
          **lipid_normal_dict** and **lipid_splay_dict** are dictionaries with keys corresponding to the elements
          in **sele_dict**. For more information about **lipid_tilt_dict** and **lipid_splay_dict**, refer to the
          the documentation for the **AnalyzeLipidTilts** and **AnalyzeLipidSplays** functions.

  WARNING: Removed parameters PBC, cell_center, cell_size
  """
  import time
  t0=time.time()
  eh=t.GetEntity()
  t.CopyFrame(0)
  #first we build the interafce between water and membrane
  #and assign a normal for each lipid in each frame
  print 'Generating the water and lipid densities and boundary surface'
  lipid_sele='rname='+','.join(lipid_names)
  water_sele='rname='+water_name
  lipids=eh.Select(lipid_sele)
  waters=eh.Select(water_sele)
  (water_filtered,lipid_filtered,b_eh)=GetBoundaryBetweenViews(t,waters,lipids,outdir,density_cutoff,
                                                               density_stride,within_size_normals,filename_basis)
  #Protein center of mass
  if prot_sele:prot_cm=mol.alg.AnalyzeCenterOfMassPos(t,eh.Select(prot_sele))
  else:prot_cm=None
  #From here we separate into the sub categories defined in sele_dict, which could be upper and lower leaflets
  s=['(rname='+ln+' and '+head_group_dict[ln]+')' for ln in lipid_names]
  lipid_sele_dict={}
  if not sele_dict:
    sele_dict['all']=''
  for sele_name in sele_dict:
    if sele_dict[sele_name]:
      lipid_sele_dict[sele_name]=' or '.join(s)+' and '+sele_dict[sele_name]
    else:
      lipid_sele_dict[sele_name]=' or '.join(s)
  print sele_dict,lipid_sele_dict
  lipid_normal_dict={}
  lipid_tilt_dict={}
  if tilt_bool_prop:lipid_tilt_dict_sele={}
  splay_dict={}
  for sele_name in sele_dict:
    lipid_sele=lipid_sele_dict[sele_name]
    sele=sele_dict[sele_name]
    sele_ev=eh.Select(lipid_sele,mol.MATCH_RESIDUES)
    for ln in lipid_names:
      v=sele_ev.Select('rname={0}'.format(ln))
      for i,r in enumerate(v.residues):r.SetIntProp('index',i)
    print 'Assigning normals for',sele_name,sele,lipid_sele
    lipid_normal_dict[sele_name]=AssignNormalsToLipids(t,sele_ev,b_eh.Select(sele),lipid_names,head_group_dict)
    print 'lipid normal dict',lipid_normal_dict[sele_name].keys()
    print 'Done in',time.time()-t0,'seconds'
    t0=time.time()
    print 'calculating tilts for',sele
    lipid_tilt_dict[sele_name]=AnalyzeLipidTilts(t,eh.Select(lipid_sele,mol.MATCH_RESIDUES),lipid_names,lipid_normal_dict[sele_name],
                                                  head_group_dict,tail_dict,prot_cm)
    if tilt_bool_prop:
      print "Caalculating tilts again only for the lipids with bool prop",tilt_bool_prop
      lipid_tilt_dict_sele[sele_name]=AnalyzeLipidTilts(t,eh.Select(lipid_sele,mol.MATCH_RESIDUES),lipid_names,lipid_normal_dict[sele_name],
                                                  head_group_dict,tail_dict,prot_cm,tilt_bool_prop)
    print 'Done in',time.time()-t0,'seconds'
    t0=time.time()
    print 'calculating splay for',sele
    splay_dict[sele_name]=AnalyzeLipidSplays(t,eh.Select(lipid_sele,mol.MATCH_RESIDUES),lipid_names,head_group_dict,tail_dict,
                                            lipid_normal_dict[sele_name],lipid_tilt_dict[sele_name],distance_sele_dict,distance_cutoff,splay_bool_prop)
    print 'Done in',time.time()-t0,'seconds'
    if outdir:
      try:WriteSplayDict(splay_dict[sele_name],outdir,filename_basis+sele_name+'_')
      except:print 'could not write dict'
      if tilt_bool_prop:WriteTiltDict(lipid_tilt_dict_sele[sele_name],outdir,filename_basis+sele_name+'_')
      else:WriteTiltDict(lipid_tilt_dict[sele_name],outdir,filename_basis+sele_name+'_')
  if tilt_bool_prop: return (lipid_tilt_dict_sele,lipid_normal_dict,splay_dict,b_eh)
  else: return (lipid_tilt_dict,lipid_normal_dict,splay_dict,b_eh)

def WriteTiltDict(lipid_tilt_dict,outdir,filename_basis=''):
  """
  This function writes out the *lipid_tilt_dict* to the firectory *outdir*.
  It writes out one file for every key in *lipid_tilt_dict*, i.e. for every lipid type.
  Filenames are preceded by the *filename_basis*.

  :param lipid_tilt_dict: the lipid tilt dictionary as produced by :func:`AnalyzeLipidTiltAndSplay`.
  :param outdir: the directory to which the files will be written
  :param filename_basis: Will be prepended to all file names.

  :type lipid_tilt_dict: :class:`dict`
  :type outdir: :class:`str`
  :type filename_basis: :class:`str`
  """
  for ln in lipid_tilt_dict:
    tl=[180.*el/math.pi  for lt in lipid_tilt_dict[ln][0] for el in lt]
    dl=[el  for lt in lipid_tilt_dict[ln][1] for el in lt]
    if len(dl)!=0:
      if len(dl)!=len(tl):print 'not the same number of tilts and distances',len(tl),len(dl)
      ll=[tl,dl]
      titles=['tilt','dist from prot']
    else:
      ll=[tl]
      titles=['tilt']
    file_utilities.WriteListOfListsInColumns(titles,ll,os.path.join(outdir,filename_basis+'tilt_'+ln+'.txt'))

def WriteSplayDict(splay_dict,outdir,filename_basis):
  """
  This function writes out the *splay_dict* to the firectory *outdir*.
  It writes out one file for every key in *splay_dict*, i.e. for every lipid type pair.
  Filenames are preceded by the *filename_basis*.

  :param splay_dict: the lipid splay dictionary as produced by :func:`AnalyzeLipidTiltAndSplay`.
  :param outdir: the directory to which the files will be written
  :param filename_basis: Will be prepended to all file names.

  :type splay_dict: :class:`dict`
  :type outdir: :class:`str`
  :type filename_basis: :class:`str`
  """
  if len(splay_dict.keys())==0:return
  for ln in splay_dict:
    if len(splay_dict[ln])==0:continue
    if len(splay_dict[ln][0])==3:prot_dist_flag=True
    else:prot_dist_flag=False
    if prot_dist_flag:
      sl=[[el[0],el[1],el[2]] for el in splay_dict[ln]]
      file_utilities.WriteListOfListsInLines(['splay','dist','ProtDist'],sl,os.path.join(outdir,filename_basis+'splay_'+ln+'.txt'))
    else:
      sl=[[el[0],el[1]] for el in splay_dict[ln]]
      file_utilities.WriteListOfListsInLines(['splay','dist'],sl,os.path.join(outdir,filename_basis+'splay_'+ln+'.txt'))
  return


def _gauss(x, *p):
  A, mu, sigma = p
  return A*npy.exp(-(x-mu)**2/(2.*sigma**2))

def _parabole(x,*p):
  a,b,x0=p
  return a+b*(x-x0)**2.0
  
def _centered_parabole(x,*p):
  a,b=p
  return a+b*x**2.0

def _FindIndexOfClosestValue(l,v):
  return min(enumerate(l), key=lambda x: abs(x[1]-v))[0]


def _FitGaussian(bincenters,pa):
  mu0=npy.sum(bincenters*pa)/npy.sum(pa)
  A0=max(pa)
  sigma0=npy.sqrt(npy.sum(((bincenters-mu0)**2.0)*pa)/npy.sum(pa))
  (A,mu,sigma),v=curve_fit(_gauss,bincenters,pa,[A0,mu0,sigma0])
  #print A0,mu0,sigma0,A,mu,abs(sigma)
  return A,mu,abs(sigma)

def _PlotGaussian(bincenters,pa,A,mu,sigma,outfile,title='',xlabel=''):
  plt.figure()
  plt.plot(bincenters,pa,'o')
  plt.plot(bincenters,_gauss(bincenters,A,mu,sigma),'-',color='g',label='$\mu={0}\ ;\ \sigma={1}$'.format(round(mu,int(1/mu)+1),round(sigma,int(1/sigma)+1)))
  plt.vlines(mu,0,plt.ylim()[1],linestyle='--',color='r')
  plt.vlines([mu+sigma,mu-sigma],0,plt.ylim()[1],linestyle='--',color='c')
  plt.title(title)
  plt.xlabel(xlabel)
  plt.ylabel('Probability')
  #plt.show()
  plt.savefig(outfile)
  plt.close()
  
def _FitParabole(bincenters,fa,fitting_range,centered=False):
  first=_FindIndexOfClosestValue(bincenters,fitting_range[0])
  last=_FindIndexOfClosestValue(bincenters,fitting_range[1])
  mask=fa!=npy.inf
  a=min(fa)
  x0=bincenters[npy.argmin(fa)]
  xm=bincenters[mask][npy.argmax(fa[mask])]
  fm=max(fa[mask])
  b=(fm-a)/(xm-x0)**2.0
  if centered:r,v=curve_fit(_centered_parabole,bincenters[first:last],fa[first:last],[a,b])
  else:r,v=curve_fit(_parabole,bincenters[first:last],fa[first:last],[a,b,x0])
  return r

def _PlotParabola(bincenters,fa,a,b,x0,fitting_range,outfile,title='',xlabel='',ylabel='',fit_label=''):
  plt.figure()
  plt.plot(bincenters,fa,'o')
  ymax=plt.ylim()[1]
  ymin=plt.ylim()[0]
  plt.plot(bincenters,[_parabole(xi,a,b,x0) for xi in bincenters],'--',color='r',label=fit_label)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.title(title)
  if fit_label:plt.legend(loc='best')
  plt.vlines(fitting_range,ymin,ymax,linestyle='--',color='c')
  plt.ylim(ymin,ymax)
  #plt.show()
  plt.savefig(outfile)
  plt.close()
  return

def ExtractTiltAndSplayModuli(tilt_dict,splay_dict,lipid_area,outdir,nbins=100):
  for sele_name in tilt_dict.keys():
    outfile=open(os.path.join(outdir,sele_name+'_tilt_constants.txt'),'w')
    k_list=FloatList()
    deltak_list=FloatList()
    nl_list=FloatList()
    for lipid_name in tilt_dict[sele_name].keys():
      fname="_".join([sele_name,lipid_name])
      tilt_list=FloatList()
      for el in tilt_dict[sele_name][lipid_name][0]:tilt_list.extend(el)
      k,dk,kl=FitTiltDistribution(tilt_list,nbins,outdir=outdir,filename_basis=fname,title_complement='for '+sele_name+" "+lipid_name)
      print "Tilt modulus for {0} {1} is k={2:1.2f} +/- {3:1.2f}.".format(sele_name,lipid_name,k,dk)
      outfile.write(' '.join([lipid_name,str(round(kl[0],1)),str(round(dk,1))])+'\n')
      k_list.append(k)
      deltak_list.append(dk)
      nl_list.append(len(tilt_list))
    ntot=npy.sum(nl_list)
    k=npy.sum([ni/ki for ni,ki in zip(nl_list,k_list)])/ntot
    k=1./k
    dk=npy.sum([((k**2.0)/(ki**2.0)*ni/ntot)**2.0*(dki**2.0) for ni,ki,dki in zip(nl_list,k_list,deltak_list)])
    dk=npy.sqrt(dk)
    outfile.write(' '.join(['agregated',str(round(k,1)),str(round(dk,1))])+'\n')
    outfile.close()
    for i in range(2*len(k_list)):plt.close()

  for sele_name in splay_dict.keys():
    outfile=open(os.path.join(outdir,sele_name+'_splay_constants.txt'),'w')
    k_list=FloatList()
    nl_list=FloatList()
    deltak_list=FloatList()
    for key in splay_dict[sele_name]:
      fname="_".join([sele_name,key])
      splay_list=[el[0] for el in splay_dict[sele_name][key]]
      k,dk,kl=FitSplayDistribution(splay_list,lipid_area,nbins,outdir=outdir,filename_basis=fname,title_complement='for '+sele_name+" "+key)
      print "Splay modulus for {0} {1} is k={2:1.2f} +/- {3:1.2f}.".format(sele_name,key,k,dk)
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
    for i in range(2*len(k_list)):plt.close()

def FitSplayDistribution(splay_list,lipid_area,nbins=100,x_range=None,outdir='',filename_basis='',title_complement=''):  
  """
  This function extracts the bending modulus from a list of splays.
  It will first fit a gaussian y=A exp[(x-mu)/sigma^2] to the distribution of splays to determine
  the fitting range used to fit the potential of mean force (PMF).
  Different fitting ranges are used to estimate the error on the extracted bending rigidity

  :param splay_list: A list of lipid splays.
  :param lipid_area: The area per lipid.
  :param nbins:   The number of bins used when determining the distribution of splays
  :param x_range:  The range in which the distribution will be calculated. Defaults to [mean-3*std,mean+3*std].
  :param outdir:  the directory to which output files will be written, i.e. plots of splay distribution and PMF.
                  if it is not defined, plots will not be generated
  :param filename_basis: Will be prepended to all file names.
  :param title_complement: will be added in the title of the plots

  :type splay_list: :class:`list`
  :type lipid_area: :class:`float`
  :type nbins: :class:`int`
  :type x_range: :class:`tuple` of 2 floats
  :type outdir: :class:`str`
  :type filename_basis: :class:`str`
  :type title_complement: :class:`str`

  :return: A tuple **(K, DeltaK,K_list)**, containing the bending rigidity *K*,
           the estimated uncertainty on *K*, and the list of *K* values obtained from the different fitting ranges.
  :rtype: (:class:`float`,:class:`float`,:class:`list`)
  """
  if not x_range:
    w=npy.std(splay_list)
    m=npy.average(splay_list)
    x_range=[m-3*w,m+3*w]
  pa=npy.histogram(splay_list,nbins,range=x_range,density=True)
  bincenters=0.5*(pa[1][1:]+pa[1][:-1])
  fa=-npy.log(pa[0])
  A,mu,sigma=_FitGaussian(bincenters,pa[0])
  ranges=[(mu-i*sigma,mu+i*sigma) for i in [1,1.25,1.5,1.75,2]]
  print "Gaussian y=A exp[(x-mu)/sigma^2] fitted to tilt distribution with A={0},mu={1} and sigma={2}".format(round(A,2),round(mu,2),round(sigma,2))
  print "Using the following ranges to fit the PMF:",ranges
  res_list=[]
  for fitting_range in ranges:
    try:r=_FitParabole(bincenters,fa,fitting_range)
    except:r=[0,0]
    res_list.append(r)
  K_list=[2.*r[1]/lipid_area for r in res_list]
  DeltaK=npy.std([el for el in K_list if not el==0.0])
  K=K_list[0]
  if outdir:
    file_utilities.WriteListOfListsInColumns(['bin','distribution','pmf'],[bincenters,pa[0],fa],os.path.join(outdir,'_'.join([filename_basis,'splay','distribution'])+'.txt'),separator=' ')
    title='Splay distribution {0}'.format(title_complement)
    outfile=os.path.join(outdir,'_'.join([filename_basis,'splay','distribution'])+'.png')
    _PlotGaussian(bincenters,pa[0],A,mu,sigma,outfile,title=title,xlabel='Splay')
    outfile=os.path.join(outdir,'_'.join([filename_basis,'splay','fit'])+'.png')
    a,b,x0=res_list[0]
    title='Splay PMF {0}'.format(title_complement)
    _PlotParabola(bincenters,fa,a,b,x0,ranges[0],outfile,title,'Splay',r'$-\log\left[P(\alpha)\right]$','$\chi^{{12}}={0}\pm {1}$'.format(round(K,int(-math.log10(K))+1),round(DeltaK,int(-math.log10(DeltaK))+1)))
  return K,DeltaK,K_list
  
def FitTiltDistribution(tilt_list,nbins=90,x_range=None,outdir='',filename_basis='',title_complement='',degrees=False):
  """
  This function extracts the tilt modulus from a list of lipid tilts.
  It will first fit a gaussian y=A exp[(x-mu)/sigma^2] to the distribution of tilts to determine
  the fitting range used to fit the potential of mean force (PMF).
  Different fitting ranges are used to estimate the error on the extracted tilt modulus

  :param tilt_list: A list of lipid splays.
  :param nbins:   The number of bins used when determining the distribution of splays
  :param x_range:  The range in which the distribution will be calculated. Defaults to [mean-3*std,mean+3*std].
  :param outdir:  the directory to which output files will be written, i.e. plots of splay distribution and PMF.
                  if it is not defined, plots will not be generated
  :param filename_basis: Will be prepended to all file names.
  :param title_complement: will be added in the title of the plots
  :param degrees: Whether plots should be in degrees or radians.

  :type tilt_list: :class:`list`
  :type nbins: :class:`int`
  :type x_range: :class:`tuple` of 2 floats
  :type outdir: :class:`str`
  :type filename_basis: :class:`str`
  :type title_complement: :class:`str`
  :type degrees: :class:`bool`

  :return: A tuple **(Chi, DeltaChi,Chi_list)**, containing the tilt modulus *Chi*,
           the estimated uncertainty on *Chi*, and the list of *Chi* values obtained from the different fitting ranges.
  :rtype: (:class:`float`,:class:`float`,:class:`list`)
  """
  if degrees:ac=math.pi/180.
  else:ac=1
  if None:range=[0,math.pi/(2.*ac)]
  pa=npy.histogram(tilt_list,nbins,range=x_range,density=True)
  bincenters=0.5*(pa[1][1:]+pa[1][:-1])
  pa=pa[0]
  pa2=pa/npy.sin(bincenters*ac)
  fa=-npy.log(pa2)
  A,mu,sigma=_FitGaussian(bincenters,pa)
  ranges=[(max(mu-i*sigma,0),mu+i*sigma) for i in [1,1.25,1.5,1.75,2.0]]
  print "Gaussian y=A exp[(x-mu)/sigma^2] fitted to tilt distribution with A={0},mu={1} and sigma={2}".format(round(A,2),round(mu,2),round(sigma,2))
  print "Using the following ranges to fit the PMF:",ranges
  res_list=[_FitParabole(bincenters,fa,fitting_range,centered=True) for fitting_range in ranges]
  K_list=[2.*r[1]/(ac*ac) for r in res_list]
  DeltaK=npy.std(K_list)
  K=K_list[0]
  if outdir:
    file_utilities.WriteListOfListsInColumns(['bin','distribution','pmf'],[bincenters,pa,fa],os.path.join(outdir,'_'.join([filename_basis,'tilt','distribution'])+'.txt'),separator=' ')
    title='Tilt distribution {0}'.format(title_complement)
    outfile=os.path.join(outdir,'_'.join([filename_basis,'tilt','distribution'])+'.png')
    _PlotGaussian(bincenters,pa,A,mu,sigma,outfile,title=title,xlabel='Tilt')
    outfile=os.path.join(outdir,'_'.join([filename_basis,'tilt','fit'])+'.png')
    a,b=res_list[0]
    r=[el for el in ranges[0]]
    title='Tilt PMF {0}'.format(title_complement)
    _PlotParabola(bincenters,fa,a,b,0.0,r,outfile,title,'Tilt',r'$-\log\left[\frac{P(\alpha)}{\sin(\alpha)}\right]$','$\chi={0}\pm {1}$'.format(round(K,int(-math.log10(K))+1),round(DeltaK,int(-math.log10(DeltaK))+1)))
  return K,DeltaK,K_list


def AnalyzeAreaPerLipid(t,lipids):
  """
  This function calculates the area per lipid simply from the number of lipids
  and the size of the simulation box.
  The area per lipid is calculated for each frame in the simulation and the average
  is returned. **This is only suitable for bilayers**

  :param t: The trajectory
  :param lipids: Selection of the lipids in the system

  :type t: :class:`~ost.mol.CoordGroupHandle`
  :type lipid: :class:`~ost.mol.EntityView`

  :return: The area per lipid
  :rtype: :class:`float`
  """
  n=lipids.GetResidueCount()
  Al=FloatList()
  for i in range(t.GetFrameCount()):
    f=t.GetFrame(i)
    cell_size=f.GetCellSize()
    Al.append(cell_size[0]*cell_size[1])
  return 2.*npy.average(Al)/float(n)

