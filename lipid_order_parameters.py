"""
Module written by Niklaus Johner (niklaus.johner@a3.epfl.ch) 01.2013
This module contains functions to determine lipid tilt and splay angles and
Calculate the elastic properties of the membrane from there
"""
try:
  from ost import *
  import time
  import numpy as npy
  import os,math
  import entity_alg,trajectory_utilities,surface_alg,file_utilities
except:
  print 'could not import at least one of the modules nedded: ost, time, numpy, os, math, entity_alg,trajectory_utilities,surface_alg,file_utilities'


__all__=('AnalyzeMolecularOrderParameters')

def _MolecularOrderParameter(v1,v2):
  a=geom.Angle(v1,v2)
  return 0.5*(3.*math.cos(a)**2.0-1)

def CalculateMolecularOrderParameter(res,aname1,aname2,aname3):
  a1=res.FindAtom(aname1)
  a2=res.FindAtom(aname2)
  a3=res.FindAtom(aname3)
  if not (a1.IsValid() and a2.IsValid() and a3.IsValid()):return
  v1=a1.pos-a2.pos
  v2=a2.pos-a3.pos
  return _MolecularOrderParameter(v1,v2)

def CalculateMolecularOrderParameters(res,aname_list):
  natoms=len(aname_list)
  if natoms<3:return
  return [CalculateMolecularOrderParameter(res,*aname_list[i:i+3]) for i in range(natoms-2)]

def AnalyzeMolecularOrderParameters(t,lipids,aname_list,return_average=True):
  """
  This function calculates the order parameter for each successive triplet of atoms
  for each lipid over a trajectory. The order parameter is calculated from the angle
  between the director vectors between two successive pairs of atoms.
  Inputs:
  t : Trajectory
  lipids : EntityView containing the lipids to be analyzed
  aname_list: An ordered list of atom names. An order parameter is calculated for each successive triple of atoms
  """
  atom_triplet_list=[]
  order_parameter_list=[]
  natoms=len(aname_list)
  if natoms<3:return
  for i in range(natoms-2):
    [aname1,aname2,aname3]=aname_list[i:i+3]
    atom_triplet_list.append([aname1,aname2,aname3])
    op=FloatList()
    for res in lipids.residues:
      a1=res.FindAtom(aname1)
      a2=res.FindAtom(aname2)
      a3=res.FindAtom(aname3)
      if not (a1.IsValid() and a2.IsValid() and a3.IsValid()):
        print res,'is missing one of',aname1,aname2,aname3
        continue
      pl1=mol.alg.AnalyzeAtomPos(t,a1.handle)
      pl2=mol.alg.AnalyzeAtomPos(t,a2.handle)
      pl3=mol.alg.AnalyzeAtomPos(t,a3.handle)
      vl1=pl1-pl2
      vl2=pl2-pl3
      op.extend([_MolecularOrderParameter(v1,v2) for v1,v2 in zip(vl1,vl2)])
    order_parameter_list.append(op)
  if return_average:return atom_triplet_list,FloatList([npy.average(el) for el in order_parameter_list])
  else: return atom_triplet_list,order_parameter_list
  
  