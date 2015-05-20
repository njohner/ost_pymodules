"""
.. codeauthor: Niklaus Johner <niklaus.johner@a3.epfl.ch>

This module contains basic functions to work with structures (Entities and EntityViews).
It notably includes several functions to work with periodic conditions in orthogonal
and non-orthogonal unit cells
"""
__all__=('')

import ost as _ost
import math as _math
import ost.mol as _mol
import ost.geom as _geom
import numpy as _npy


def DistanceBetweenTwoAngles(angle1,angle2,period=2.0*_math.pi):
  return abs(angle1-angle2)%(0.5*period)

def _WrapAngle2(angle,center,period=2.0*_math.pi,half_period=_math.pi):
  if angle>center+half_period:angle=angle-period
  elif angle<center-half_period:angle=angle+period
  return angle

def WrapAngle(angle,center,period=2.0*_math.pi):
  half_period=period/2.
  return _WrapAngle2(angle,center,period,half_period)

def WrapAngles(angles,center,period=2.0*_math.pi):
  half_period=period/2.
  for i,a in enumerate(angles):
    angles[i]=_WrapAngle2(a,center,period,half_period)
  return

def ResidueSymmetryDihedrals():
  dihedral_dict={}
  dihedral_dict.update({'PHE':[[('CA','CB','CG','CD1'),_math.pi,[('CD1','CD2'),('CE1','CE2')]]]})
  dihedral_dict.update({'TYR':[[('CA','CB','CG','CD1'),_math.pi,[('CD1','CD2'),('CE1','CE2')]]]})
  dihedral_dict.update({'ASP':[[('CA','CB','CG','OD1'),_math.pi,[('OD1','OD2')]]]})
  dihedral_dict.update({'GLU':[[('CB','CG','CD','OE1'),_math.pi,[('OE1','OE2')]]]})
  dihedral_dict.update({'LEU':[[('CA','CB','CG','CD1'),_math.pi,[('CD1','CD2')]]]})
  dihedral_dict.update({'ARG':[[('CD','NE','CZ','NH1'),_math.pi,[('NH1','NH2')]]]})
  dihedral_dict.update({'VAL':[[('N','CA','CB','CG1'),_math.pi,[('CG1','CG2')]]]})
  return dihedral_dict


def MatchResidueSymmetries(ref_ev,ev,dihedral_dict=None):
  if not dihedral_dict:dihedral_dict=ResidueSymmetryDihedrals()
  edi=ev.handle.EditXCS(_mol.BUFFERED_EDIT)
  for r1,r2 in zip(ref_ev.residues,ev.residues):
    r1=r1.handle
    r2=r2.handle
    if not r1.name in dihedral_dict.keys():continue
    for (an1,an2,an3,an4),period,atom_pairs in dihedral_dict[r1.name]:
      p1=r1.FindAtom(an1).pos
      p2=r1.FindAtom(an2).pos
      p3=r1.FindAtom(an3).pos
      p4=r1.FindAtom(an4).pos
      a1=r2.FindAtom(an1)
      a2=r2.FindAtom(an2)
      a3=r2.FindAtom(an3)
      a4=r2.FindAtom(an4)
      center=_geom.DihedralAngle(p1,p2,p3,p4)
      a=_geom.DihedralAngle(a1.pos,a2.pos,a3.pos,a4.pos)
      wa=WrapAngle(a,center,period)
      if a!=wa:
        for anp1,anp2 in atom_pairs:
          a1=r2.FindAtom(anp1)
          a2=r2.FindAtom(anp2)
          p1=a1.pos
          p2=a2.pos
          edi.SetAtomPos(a1,p2)
          edi.SetAtomPos(a2,p1)
  edi.ForceUpdate()
  return

def _CreateViewWithMostProbableDihedrals(t,ev,dihedral_dict):
  ref_ev=_mol.CreateEntityFromView(ev,True)
  edi=ref_ev.EditICS()
  for r1,r2 in zip(ev.resiudes,ref_ev.residues):
    if not r1.name in dihedral_dict.keys():continue
    for (an1,an2,an3,an4),period in dihedral_dict[r1.name]:
      a11=r1.FindAtom(an1)
      a21=r1.FindAtom(an2)
      a31=r1.FindAtom(an3)
      a41=r1.FindAtom(an4)
      angles=_mol.alg.AnalyzeDihedralAngle(t,a1,a2,a3,a4)
      h=_npy.histogram(angles,range=(0,2*_math.pi),bins=20)
      imax=_npy.argmax(h[0])
      center=0.5*(h[1][imax]+h[1][imax+1])
      a12=r2.FindAtom(an1)
      a22=r2.FindAtom(an2)
      a32=r2.FindAtom(an3)
      a42=r2.FindAtom(an4)
      edi.SetTorsionAngle(a12,a22,a32,a42,center)
  edi.UpdateXCS()
  return ref_ev

def CorrectResidueSymmetries(t,ev,ref_ev=None,dihedral_dict=None):
  if not dihedral_dict:dihedral_dict=ResidueSymmetryDihedrals()
  if ref_ev:ref_ev=_mol.CreateEntityFromView(ref_ev,True)
  else:ref_ev=_CreateViewWithMostProbableDihedrals(t,ev,dihedral_dict)
  
  for i in range(t.GetFrameCount()):
    t.CopyFrame(i)
    MatchResidueSymmetries(ref_ev,ev,dihedral_dict)
    t.Capture(i)
  return
    

"""
def ResidueSymmetryDihedrals():
  dihedral_dict={}
  dihedral_dict.update({'PHE':[[('CA','CB','CG','CD1'),_math.pi]]})
  dihedral_dict.update({'TYR':[[('CA','CB','CG','CD1'),_math.pi]]})
  dihedral_dict.update({'ASP':[[('CA','CB','CG','OD1'),_math.pi]]})
  dihedral_dict.update({'GLU':[[('CB','CG','CD','OE1'),_math.pi]]})
  dihedral_dict.update({'LEU':[[('CA','CB','CG','CD1'),_math.pi]]})
  dihedral_dict.update({'ARG':[[('CD','NE','CZ','NH1'),_math.pi]]})
  dihedral_dict.update({'VAL':[[('N','CA','CB','CG1'),_math.pi]]})
  return dihedral_dict
"""

"""
def CorrectResidueSymmetry(t,a1,a2,a3,a4,period=2.0*_math.pi,center=None):
  angles=_mol.alg.AnalyzeDihedralAngle(t,a1,a2,a3,a4)
  if not center:
    h=_npy.histogram(angles,range=(0,2*_math.pi),bins=20)
    imax=_npy.argmax(h[0])
    center=0.5*(h[1][imax]+h[1][imax+1])
  WrapAngles(angles,center,period)
  eh=t.GetEntity()
  edi=eh.EditICS()
  for i in range(t.GetFrameCount()):
    t.CopyFrame(i)
    edi.SetTorsionAngle(a1,a2,a3,a4,angles[i],True)
    edi.UpdateXCS()
    t.Capture(i)
  return
"""
"""  
def CorrectResidueSymmetries(t,ev,ref_ev=None,dihedral_dict=None):
  if not dihedral_dict:dihedral_dict=ResidueSymmetryDihedrals()
  if ref_ev:
    centers=[]
    for r in ref_ev.residues:
      r=r.handle
      if not r.name in dihedral_dict.keys():
        centers.append(_ost.FloatList())
        continue
      fl=_ost.FloatList()
      for (an1,an2,an3,an4),period in dihedral_dict[r.name]:
        a1=r.FindAtom(an1)
        a2=r.FindAtom(an2)
        a3=r.FindAtom(an3)
        a4=r.FindAtom(an4)
        fl.append(_geom.DihedralAngle(a1.pos,a2.pos,a3.pos,a4.pos))
      centers.append(fl)
  for i,r in enumerate(ev.residues):
    r=r.handle
    if not r.name in dihedral_dict.keys():continue
    for j,((an1,an2,an3,an4),period) in enumerate(dihedral_dict[r.name]):
      a1=r.FindAtom(an1)
      a2=r.FindAtom(an2)
      a3=r.FindAtom(an3)
      a4=r.FindAtom(an4)
      if ref_ev:CorrectResidueSymmetry(t,a1,a2,a3,a4,period,centers[i][j])
      else:CorrectResidueSymmetry(t,a1,a2,a3,a4,period)
  return
"""
"""
def CorrectResidueSymmetries2(ev1,ev2,dihedral_dict=None):
  if not dihedral_dict:dihedral_dict=ResidueSymmetryDihedrals()
  edi=ev2.handle.EditICS(_mol.BUFFERED_EDIT)
  for r1,r2 in zip(ev1.residues,ev2.residues):
    r1=r1.handle
    r2=r2.handle
    if not r1.name in dihedral_dict.keys():continue
    for (an1,an2,an3,an4),period in dihedral_dict[r1.name]:
      p1=r1.FindAtom(an1).pos
      p2=r1.FindAtom(an2).pos
      p3=r1.FindAtom(an3).pos
      p4=r1.FindAtom(an4).pos
      a1=r2.FindAtom(an1)
      a2=r2.FindAtom(an2)
      a3=r2.FindAtom(an3)
      a4=r2.FindAtom(an4)
      center=_geom.DihedralAngle(p1,p2,p3,p4)
      a=WrapAngle(_geom.DihedralAngle(a1.pos,a2.pos,a3.pos,a4.pos),center,period)
      edi.SetTorsionAngle(a1,a2,a3,a4,a,True)
  edi.UpdateXCS()
  return
"""