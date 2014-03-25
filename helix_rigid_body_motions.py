"""
Module written by Niklaus Johner (niklaus.johner@a3.epfl.ch) 01.2014
This module contains functions to analyze the rigid body motions of helices in proteins
"""
from ost import geom as _geom
from ost import mol as _mol
import ost as _ost
import numpy as _npy

def FindLineSuperpositionTransformation(line1,line2):
  trans=line1.origin-line2.origin
  rot_axis=_geom.Cross(line2.direction,line1.direction)
  rot_angle=_geom.Angle(line1.direction,line2.direction)
  rot_q=_geom.Quat(rot_angle,rot_axis)
  return (trans,rot_q)
  
def FindVectorFromPointToLine(v,l):
  d=v-l.origin
  d=d-_geom.Dot(d,l.direction)*l.direction
  return d

def TransformLine(line,trans,rot_q):
  origin=line.origin+trans
  direction=rot_q.ToRotationMatrix()*line.direction
  return geom.Line3(origin,origin+direction)

def RotateVec3List(vl,rot_q):
  R=rot_q.ToRotationMatrix()
  vl2=_geom.Vec3List([R*v for v in vl])
  return vl2
  
def FindHelixSuperpositionTransformation(helix1,helix2,return_transform=False):
  l1=_mol.alg.structure_analysis.CalculateHelixAxis(helix1)
  l2=_mol.alg.structure_analysis.CalculateHelixAxis(helix2)
  (helix_trans,helix_rot)=FindLineSuperpositionTransformation(l1,l2)
  vl1=_geom.Vec3List([FindVectorFromPointToLine(a.pos,l1) for a in helix1.atoms])
  vl2=_geom.Vec3List([FindVectorFromPointToLine(a.pos,l2) for a in helix2.atoms])
  vl2=RotateVec3List(vl2,helix_rot)
  (transz,rotz)=FindLineSuperpositionTransformation(_geom.Line3(_geom.Vec3(),_geom.Vec3(0,0,1)),l1)
  vl1=RotateVec3List(vl1,rotz)
  vl2=RotateVec3List(vl2,rotz)
  angles=_ost.FloatList([_geom.SignedAngle(_geom.Vec2(v1),_geom.Vec2(v2)) for v1,v2 in zip(vl1,vl2)])
  axis_rot_angle=-_npy.average(angles)  
  axis_rot=_geom.Quat(axis_rot_angle,l1.direction)
  if return_transform:
    T1=_geom.Transform()
    T1.SetTrans(-l2.origin)
    T2=_geom.Transform()
    T2.SetTrans(l2.origin+helix_trans)
    HR=_geom.Transform()
    HR.SetRot(helix_rot.ToRotationMatrix())
    AR=_geom.Transform()
    AR.SetRot(axis_rot.ToRotationMatrix())
    T=_geom.Transform()
    T.SetMatrix(T2.matrix*AR.matrix*HR.matrix*T1.matrix)
    return (helix_trans,helix_rot,axis_rot,T)
  else:return (helix_trans,helix_rot,axis_rot)

def AnalyzeHelixRigidBodyMotion(t,helix,ref_helix,stride=1):
  helix_rot_angles=_ost.FloatList()
  helix_rot_axes=_geom.Vec3List()
  axis_rot_angles=_ost.FloatList()
  axes=_geom.Vec3List()
  transform_list=[]
  for i in range(0,t.GetFrameCount(),stride):
    t.CopyFrame(i)
    (helix_trans,helix_rot,axis_rot,T)=FindHelixSuperpositionTransformation(ref_helix,helix,True)
    helix_rot_angles.append(helix_rot.GetAngle())
    helix_rot_axes.append(helix_rot.GetAxis())
    axis_rot_angles.append(axis_rot.GetAngle())
    axes.append(axis_rot.GetAxis())
    transform_list.append(T)
  axis_rot_axis=_mol.alg.structure_analysis.CalculateHelixAxis(ref_helix).direction
  axis_rot_angles=[_geom.Dot(axis_rot_axis,v)*a for a,v in zip(axis_rot_angles,axes)]
  return helix_rot_angles,axis_rot_angles,helix_rot_axes,axis_rot_axis,transform_list

def CalculateHelixRigidBodyMotion(helix,ref_helix,stride=1):
  (helix_trans,helix_rot,axis_rot,T)=FindHelixSuperpositionTransformation(ref_helix,helix,True)
  axis_rot_axis=_mol.alg.structure_analysis.CalculateHelixAxis(ref_helix).direction
  axis_rot_angle=_geom.Dot(axis_rot_axis,axis_rot.GetAxis())*axis_rot.GetAngle()
  return helix_rot.GetAngle(),axis_rot_angle,helix_rot.GetAxis(),axis_rot_axis,T


def GetPossibleSeparationList(n_ele,n_sep):
  """
  Returns all the possible separations of a list
  with n_ele elements into n_sep+1 pieces
  """
  if n_sep<1:return []
  indices_list=[[i] for i in range(1,n_ele+1-n_sep)]
  for j in range(1,n_sep):
    indices_new=[]
    for el in indices_list:
      for i in range(j,n_ele+1+j-n_sep):
        if i<=el[-1]:continue
        indices_new.append([k for k in el])
        indices_new[-1].append(i)
    indices_list=[el for el in indices_new]
  return indices_list

def CreateViewFromResidueViewList(rl):
  al=[]
  for r in rl:
    al.extend(r.GetAtomList())
  view=_mol.CreateViewFromAtoms(al)
  return view
 
def FindBestStructuralDecomposition(t,ev,n_pieces):
  """
  Finds the decomposition into n_pieces (continuous in sequence)
  that minimize the total rmsf.
  t : trajectory
  ev : EntityView on which the decomposition is performed
  n_pieces: The number of pieces into which to decompose ev.
  """
  n_ele=ev.GetResidueCount()
  separation_list=GetPossibleSeparationList(n_ele,n_pieces-1)
  for el in separation_list:
    el.insert(0,0)
    el.append(n_ele)
  if len(separation_list)==0:separation_list=[[0,n_ele]]
  rmsf_list=_ost.FloatList()
  for indices in separation_list:
    sele_list=[CreateViewFromResidueViewList(ev.residues[indices[i]:indices[i+1]]) for i in range(len(indices)-1)]
    rmsf=0
    for sele in sele_list:
      t2=_mol.alg.SuperposeFrames(t,sele)
      rmsfi=_mol.alg.AnalyzeRMSF(t2,sele)
      rmsf+=rmsfi*rmsfi*sele.GetResidueCount()
    rmsf_list.append(rmsf)
  min_index=_npy.argmin(rmsf_list)
  indices=separation_list[min_index]
  sele_list=[CreateViewFromResidueViewList(ev.residues[indices[i]:indices[i+1]]) for i in range(len(indices)-1)]
  rmsf=(rmsf_list[int(min_index)]/float(n_ele))**0.5
  return sele_list,rmsf
  
  
  
 