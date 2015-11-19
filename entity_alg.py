"""
.. codeauthor:: Niklaus Johner <niklaus.johner@a3.epfl.ch>

This module contains basic functions to work with structures (Entities and EntityViews).
It notably functions to generate the biounit or crystal packing from a pdb file
"""
__all__=('FindClosestAtom','CreateEntityFromVec3List',"GenerateBioUnit"\
         'WriteFloatPropList','AssignFloatPropToEntityResidues','AssignFloatPropToEntity','GenerateCrystalPackingFromPDB')
from ost import *
import math
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

def WriteVec3Prop(eh,prop_name,filename,index=False,delimiter=' ',precision='4',column_width=10):
  with open(filename,'w') as f: 
    if index:f.write(('{:^'+str(column_width)+'} ').format('aindex'))
    col_names=[prop_name+el for el in [".x",".y",".z"]]
    format=' '.join(['{:^'+str(column_width)+'}' for key in col_names])
    format+='\n'
    f.write(format.format(*col_names))
    if index:format='{:^'+str(column_width)+'} '+' '.join(['{:'+str(column_width)+'.'+str(precision)+'g}' for key in col_names])
    else:format=' '.join(['{:'+str(column_width)+'.'+str(precision)+'g}' for key in col_names])
    format+='\n'
    for a in eh.atoms:
      v=a.GetVec3Prop(prop_name)
      l=[v.x,v.y,v.z]
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
  try: import pbc_utilities
  except:
    print "could not load the pbc_utilities module"
    print "This module is needed in the function GenerateCrystalPackingFromPDB"
    return None
  reload(file_utilities)
  pdb=open(pdb_filename,'r')
  trans=file_utilities.ReadSymmetryFromPDB(pdb)
  uc=file_utilities.ReadUnitCellFromPDB(pdb)
  (v1,v2,v3)=pbc_utilities.BuildUnitCellVectors(uc[0],uc[1],uc[2],uc[3],uc[4],uc[5])
  vl=pbc_utilities.GenerateSearchVectors(v1,v2,v3,dist_cutoff=vec_search_cutoff,level=vec_search_level)
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
  