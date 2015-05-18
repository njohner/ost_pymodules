"""
Module written by Niklaus Johner (niklaus.johner@a3.epfl.ch) 01.2013

This module contains basic functions to work with structures (Entities and EntityViews).
It notably includes several functions to work with periodic conditions in orthogonal
and non-orthogonal unit cells
"""
__all__=('FindClosestAtom','CreateEntityFromVec3List','ProjectVectorOnUnitCellVectors',\
         'GetUCellVectorSizeAndAngles','FindWithinWithPBC','FindWithinWithNonOrthogonalPBC',\
         'ExtendEntityToNeighboringUnitCells','ExtendEntityWithPBC','WriteFloatPropList',\
         'AssignFloatPropToEntityResidues','AssignFloatPropToEntity','BuildUnitCellVectors',\
         'VectorsToNeighboringUnitCells','GenerateSearchVectors','GenerateCrystalPackingFromPDB')
from ost import *
import math
import numpy as npy
import scipy
from scipy import optimize
import time
import sys
sys.path.append('/work/python_modules')
import file_utilities

def FindClosestAtom(v1,v2):
  """
  v1 and v2 need to be an Entity, EntityView, ChainView, ResidueView, ChainHandle or ResidueHandle
  Returns the atom of v2 that is closest to v1 (minimal distance to any atom in v1).
  """
  if not type(v2)==mol.EntityView:v2=mol.CreateViewFromAtoms(v2.atoms)
  if not type(v1)==mol.EntityView:v1=mol.CreateViewFromAtoms(v1.atoms)
  (i,j)=geom.MinDistanceIndices(mol.alg.GetPosListFromView(v1),mol.alg.GetPosListFromView(v2))
  return v2.atoms[j]


def CreateEntityFromVec3List(pos_list,chain_name='A',rname='b',aname_base='C',atom_element='C'):
  """
  Creates an entity with atoms positioned according to vectors in the pos_list (Vec3List).
  All atoms are generated in a single chain, with 100 atoms per residue.
  """
  eh=mol.CreateEntity()
  edi=eh.EditXCS(mol.BUFFERED_EDIT)
  c=edi.InsertChain(chain_name)
  r=edi.AppendResidue(c,rname)
  ac=0
  for v in pos_list:
    edi.InsertAtom(r,aname_base+str(ac),v,atom_element)
    ac+=1
    if ac==100:
      ac=0
      r=edi.AppendResidue(c,'b')
  edi.ForceUpdate()
  return eh


def ProjectVectorOnUnitCellVectors(v,ucell_vecs):
  """
  This function project the vector v on the unit cell vectors
  ucell_vecs, of which the first vector is assumed to be along x
  and the second in the xy plane.
  """
  v_projected=geom.Vec3()
  v_projected[2]=v[2]/ucell_vecs[2][2]
  v_projected[1]=(v[1]-v_projected[2]*ucell_vecs[2][1])/ucell_vecs[1][1]
  v_projected[0]=(v[0]-v_projected[2]*ucell_vecs[2][0]-v_projected[1]*ucell_vecs[1][0])/ucell_vecs[0][0]
  return v_projected
  
def GetUCellVectorSizeAndAngles(ucell_vecs):
  """
  Returns the length of the unit cell vectors and the angles alpha, beta and gamma
  The first vector has to be along x and the second in the xy plane.
  First angle is between v1 and v2, the second one between v1 and v3 and the third one between v2 and v3.
  """
  if any([ucell_vecs[0][1]!=0,ucell_vecs[0][2]!=0,ucell_vecs[1][2]!=0]):
    print 'Problem with unit cell vectors, first vector has to be along x and second in the xy plane'
    return None
  ucell_size=geom.Vec3(geom.Length(ucell_vecs[0]),geom.Length(ucell_vecs[1]),geom.Length(ucell_vecs[2]))
  ucell_angles=geom.Vec3(geom.Angle(ucell_vecs[0],ucell_vecs[1]),geom.Angle(ucell_vecs[0],ucell_vecs[2]),geom.Angle(ucell_vecs[1],ucell_vecs[2]))
  return(ucell_size,ucell_angles)
  
def FindWithinWithPBC(eh,pos,radius,cell_center,cell_size):
  """
  Returns a view containing all the atoms within *radius* of a point
  *pos* (Vec3), using periodic boundary conditions (orthogonal unit cell).
  """
  alist=eh.FindWithin(pos,radius)
  within=mol.CreateViewFromAtomList(alist)
  vl=VectorsToNeighboringUnitCells(cell_size)
  vl2=geom.Vec3List()
  for v in vl:
    d=cell_center+v-pos
    if all (di<csi/2.0+radius for di,csi in zip(d.data,cell_size.data)):
      vl2.append(v)
  for v in vl2:
    p=pos+v
    alist.extend(eh.FindWithin(p,radius))
  if len(alist)>0:
    within=mol.CreateViewFromAtomList(alist)
    return within
  else:return eh.CreateEmptyView()


def FindWithinWithNonOrthogonalPBC(eh,pos,radius,ucell_vecs,vecs_to_neighbor_ucells=None):
  """
  Returns a view containing all the atoms within *radius* of a point
  *pos* (Vec3), using periodic boundary conditions (orthogonal and non-orthogonal unit cell).
  By default it searches in all neighboring unit cells (first shell of neighboring unit cells),
  except if a list of vectors to the neighboring unit cells in which to search is passed.
  """
  alist=eh.FindWithin(pos,radius)
  if not vecs_to_neighbor_ucells:vecs_to_neighbor_ucells=VectorsToNeighboringUnitCells(ucell_vecs)
  for v in vecs_to_neighbor_ucells:alist.extend(eh.FindWithin(pos+v,radius))
  if len(alist)>0:
    within=mol.CreateViewFromAtomList(alist)
    return within.Select('')
  else:return eh.CreateEmptyView()

def ExtendEntityToNeighboringUnitCells(eh,vecs_to_neighbor_ucells):
  """
  Extends an Entity to surrounding unit cells by replicating it and translating
  it for each vector in vecs_to_neighbor_ucells. 
  For each replica the chains will have the same name as in the original Entity followed by a number
  corresponding to the index of the vector used for the translation (starting at 1). 
  The translations can be obtained for example from the function VectorsToNeighboringUnitCells().
  
  :param eh: The *Entity* or *EntityView* to be extended
  :param vecs_to_neighbor_ucells: List of translations used to extend the entity.

  :type eh: :class:`~ost.mol.EntityView`
  :type vecs_to_neighbor_ucells: :class:`~ost.geom.Vec3List`
  """
  if type(eh)==mol.EntityView:eh=mol.CreateEntityFromView(eh,1)
  eh_extended=eh.Copy()
  eh_t=eh.Copy()
  edi=eh_t.EditXCS()
  edi2=eh_extended.EditXCS()
  T=geom.Transform()
  for i,v in enumerate(vecs_to_neighbor_ucells):
    T.SetTrans(v)
    edi.SetTransform(T)
    for c in eh_t.chains:
      edi2.InsertChain(c.name+str(i+1),c,1)
  return eh_extended
  

def ExtendEntityWithPBC(eh,cell_center,ucell_vecs,extension_size=10):
  """
  This function extends an entity to its neighboring unit cell vectors,
  by a size given by the *extension_size* in each direction.
  It's slower than ExtendEntityToNeighboringUnitCells but allows
  to extend to small parts of the surrounding unit cells.
  """
  vl=VectorsToNeighboringUnitCells(ucell_vecs)
  eh_extended=ExtendEntityToNeighboringUnitCells(eh,vl)
  pos_list=eh_extended.GetAtomPosList()
  (ucell_size,ucell_angles)=GetUCellVectorSizeAndAngles(ucell_vecs)
  ucell_size+=geom.Vec3(extension_size,extension_size,extension_size)
  mol.alg.WrapEntityInPeriodicCell(eh_extended,cell_center,ucell_size,ucell_angles)
  edi=eh_extended.EditXCS(mol.BUFFERED_EDIT)
  for p,a in zip(pos_list,eh_extended.atoms):
    if not p==a.pos:edi.DeleteAtom(a)
  for r in eh_extended.residues:
    if r.GetAtomCount()==0:edi.DeleteResidue(r)
  return eh_extended
    
  
def WriteFloatPropList(eh,prop_list,filename,index=False,delimiter=' ',precision='4',column_width=10):
  with open(filename,'w') as f: 
    if index:f.write(('{:^'+str(column_width)+'} ').format('aindex'))
    format=' '.join(['{:^'+str(column_width)+'}' for key in prop_list])
    format+='\n'
    f.write(format.format(*prop_list))
    if index:format='{:^'+str(column_width)+'} '+' '.join(['{:'+str(column_width)+'.'+str(precision)+'g}' for key in prop_list])
    else:format=' '.join(['{:'+str(column_width)+'.'+str(precision)+'g}' for key in prop_list])
    format+='\n'
    for a in eh.atoms:
      l=[a.GetFloatProp(key) for key in prop_list]
      if index:l.insert(0,a.index)
      f.write(format.format(*l))
  return

def AssignFloatPropToEntityResidues(eh,fl,key):
  """
  That function will not work if the atom indices are not the same as the order of atoms in the entity
  """
  for r,f in zip(eh.residues,fl):
    r.SetFloatProp(key,f)
  return

   
def AssignFloatPropToEntity(eh,fl,key,aindex_list=None):
  """
  That function will not work if the atom indices are not the same as the order of atoms in the entity
  """
  if not aindex_list:
    for a,f in zip(eh.atoms,fl):
      a.SetFloatProp(key,f)
  else:
    for ai,f in zip(aindex_list,fl):
      eh.atoms[ai].SetFloatProp(key,f)
  return


def BuildUnitCellVectors(a,b,c,alpha,beta,gamma):
  """
  Returns the unit cell vectors from their length and angles
  """
  v1=geom.Vec3(a,0.0,0.0)
  v2=geom.Vec3(b*math.cos(gamma),b*math.sin(gamma),0.0)
  v31=c*math.cos(beta)
  v32=c*((math.cos(alpha)-math.cos(beta)*math.cos(gamma))/math.sin(gamma))
  v33=c*(math.sqrt(1-math.pow(math.cos(alpha),2.0)-math.pow(math.cos(beta),2.0)-math.pow(math.cos(gamma),2.0)+2*math.cos(alpha)*math.cos(beta)*math.cos(gamma))/math.sin(gamma))
  v3=geom.Vec3(v31,v32,v33)
  return (v1,v2,v3)

def VectorsToNeighboringUnitCells(ucell_vecs,level=1):
  """
  Returns the vectors to all the neighboring unit cells.
  level=1 will the vectors to the 26 cells surrounding the central cell,
  level=2 all the vectors to the first and second shell (124 vectors), etc.
  """
  vl=[i*ucell_vecs[0]+j*ucell_vecs[1]+k*ucell_vecs[2] for i in range(-level,level+1) for j in range(-level,level+1) for k in range(-level,level+1)]
  vl.remove(geom.Vec3())
  return vl
  
def GenerateSearchVectors(v1,v2,v3,level=4,dist_cutoff=None,delta=0.0001):
  """
  Returns all the vectors to unit cells within a certain distance of the central cell.
  The distance is the minimal distance between the neighboring cell and the central cell.
  """
  if not dist_cutoff:dist_cutoff=0.5*max(geom.Length(v1),geom.Length(v2),geom.Length(v3))
  corners=[i*v1+j*v2+k*v3 for i in [-0.5,0.5] for j in [-0.5,0.5] for k in [-0.5,0.5]]
  vl_b=[v1,v2,v3,-v1,-v2,-v3]
  vl=[el for el in vl_b]
  for i in range(level):
    for v1 in vl_b:
      for v2 in vl_b:
        v=v1+v2
        if all([geom.Length(v+c1+c2)>dist_cutoff for c1 in corners for c2 in corners]):continue
        if any(geom.Length(v3-v)<delta for v3 in vl):continue
        vl.append(v1+v2)
    vl_b=[el for el in vl]
  v=geom.Vec3(0,0,0)
  if all(geom.Length(v3-v)>delta for v3 in vl_b):vl_b.append(geom.Vec3(0,0,0))
  return vl_b

def GenerateBioUnit(pdb_filename,biounit_id=1):
  """
  This function generates the biological unit from a pdb file
  by applying the transformations found in the REMARK 350 record of the PDB
  """
  pdb=open(pdb_filename,'r')
  Tl,cnames=file_utilities.FindBioUnitTransformations(pdb,biounit_id)
  eh=io.LoadPDB(pdb_filename)
  eh2=mol.CreateEntity()
  edi=eh2.EditXCS()
  for i,T in enumerate(Tl):
    for cname in cnames:
      if i>0:cname_new=cname+str(i)
      else:cname_new=cname
      c=eh.FindChain(cname)
      eh.SetTransform(T)
      edi.InsertChain(cname_new,c,deep=True)
  return eh2

def GenerateCrystalPackingFromPDB(pdb_filename,distance_cutoff=20,vec_search_level=4,vec_search_cutoff=None,superpose_sele='aname=CA'):
  """
  This function loads a pdb and generates the crystal packing around the initial structure, including
  all structures within *distance_cutoff* of the initial structure.
  ref_sele is used to test if a certain structure coming form symmetry operations and translations
  is already present in the crystal packing being generated.
  """
  reload(file_utilities)
  pdb=open(pdb_filename,'r')
  trans=file_utilities.ReadSymmetryFromPDB(pdb)
  uc=file_utilities.ReadUnitCellFromPDB(pdb)
  (v1,v2,v3)=BuildUnitCellVectors(uc[0],uc[1],uc[2],uc[3],uc[4],uc[5])
  vl=GenerateSearchVectors(v1,v2,v3,dist_cutoff=vec_search_cutoff,level=vec_search_level)
  ref_eh=io.LoadPDB(pdb_filename)
  ref_sele=ref_eh.Select(superpose_sele)
  ref_CM=ref_eh.GetCenterOfMass()
  bound_size=max(ref_eh.GetBoundarySize())
  pdb_list=[]
  pdb_list.append(ref_eh)
  eh=ref_eh.Copy()
  sele=eh.Select(superpose_sele)
  edi=eh.EditXCS()
  n_eh=1
  for j,v in enumerate(vl):  
    for i,t2 in enumerate(trans):
      t=geom.Mat4()
      t.PasteRotation(t2.ExtractRotation())
      t.PasteTranslation(t2.ExtractTranslation()+v)
      edi.SetTransform(t)
      if geom.Length(ref_CM-eh.GetCenterOfMass())>1.5*bound_size:continue
      flag=0
      for eh2 in pdb_list:
        rmsd=mol.alg.CalculateRMSD(eh2.Select('aname=CA'),sele)
        if rmsd<0.1:
          flag=1
          break
      if flag==1:continue
      if geom.MinDistance(mol.alg.GetPosListFromView(ref_sele),mol.alg.GetPosListFromView(sele))>20.0:continue
      n_eh+=1
      for c in eh.chains:
        edi.RenameChain(c,c.name+str(n_eh))
      pdb_list.append(eh)
      eh=ref_eh.Copy()
      sele=eh.Select('aname=CA')
      edi=eh.EditXCS()
  #Now we build the final view
  view=pdb_list[0].CreateEmptyView()
  for eh in pdb_list:
    for c in eh.chains:
      view.AddChain(c,7)
  eh=mol.CreateEntityFromView(view,1)
  return eh
  