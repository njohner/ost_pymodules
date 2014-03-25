"""
Module written by Niklaus Johner (niklaus.johner@a3.epfl.ch) 01.2014
This module contains functions to analyze the rigid body motions of helices in proteins
"""
from ost import geom as _geom
from ost import mol as _mol
import ost as _ost
import numpy as _npy

def SignedAngle(v1,v2,plane_normal):
  a=_geom.Angle(v1,v2)
  if a==0.0:return a
  n=_geom.Cross(v1,v2)
  s=_geom.Dot(n,plane_normal)
  return s/abs(s)*a

def GetAxisRotAngleBetweenVec3(v1,v2,axis):
  n1=_geom.Cross(v1,axis)
  n2=_geom.Cross(v2,axis)
  return SignedAngle(n2,n1,axis)
"""
def FindVec3SuperpositionRotations(v1,v2,axis1,axis2=None):
  if not axis2:axis2=_geom.Cross(v1,axis1)
  a1=GetAxisRotAngleBetweenVec3(v1,v2,axis1)
  q1=_geom.Quat(a1,axis1)
  vt=q1.ToRotationMatrix()*v2
  a2=GetAxisRotAngleBetweenVec3(v1,vt,axis2)
  q2=_geom.Quat(a2,axis2)
  return (a1,a2,q1,q2)
"""
def FindVec3SuperpositionRotations(v1,v2,axis1,axis2=None):
  if not axis2:axis2=_geom.Cross(v1,axis1)
  a2=GetAxisRotAngleBetweenVec3(v1,v2,axis2)
  q2=_geom.Quat(a2,axis2)
  vt=q2.ToRotationMatrix()*v2
  a1=GetAxisRotAngleBetweenVec3(v1,vt,axis1)
  q1=_geom.Quat(a1,axis1)
  return (a1,a2,q1,q2)


def FindVectorFromPointToLine(v,l):
  d=v-l.origin
  d=d-_geom.Dot(d,l.direction)*l.direction
  return d

def RotateVec3List(vl,rot_q):
  R=rot_q.ToRotationMatrix()
  vl2=_geom.Vec3List([R*v for v in vl])
  return vl2

def FindHelixSuperpositionTransformation(helix1,helix2,axis1,axis2=None,return_transform=False):
  l1=_mol.alg.structure_analysis.CalculateHelixAxis(helix1)
  if not axis2:axis2=_geom.Cross(l1.direction,axis1)
  l2=_mol.alg.structure_analysis.CalculateHelixAxis(helix2)
  helix_trans=l1.origin-l2.origin
  (axis1_rot_angle,axis2_rot_angle,axis1_rot,axis2_rot)=FindVec3SuperpositionRotations(l1.direction,l2.direction,axis1,axis2)
  #print axis1,axis1_rot.GetAxis(),axis1_rot.GetAngle(),axis2,axis2_rot.GetAxis(),axis2_rot.GetAngle()
  vl1=_geom.Vec3List([FindVectorFromPointToLine(a.pos,l1) for a in helix1.atoms])
  vl2=_geom.Vec3List([FindVectorFromPointToLine(a.pos,l2) for a in helix2.atoms])
  #vl2=RotateVec3List(vl2,axis1_rot)
  #vl2=RotateVec3List(vl2,axis2_rot)
  vl2=RotateVec3List(vl2,axis2_rot)
  vl2=RotateVec3List(vl2,axis1_rot)
  #angles=_ost.FloatList([_geom.SignedAngle(_geom.Vec2(v1),_geom.Vec2(v2)) for v1,v2 in zip(vl1,vl2)])
  angles=_ost.FloatList([SignedAngle(v1,v2,l1.direction) for v1,v2 in zip(vl1,vl2)])
  helix_axis_rot_angle=-_npy.average(angles)  
  helix_axis_rot=_geom.Quat(helix_axis_rot_angle,l1.direction)
  if return_transform:
    T1=_geom.Transform()
    T1.SetTrans(-l2.origin)
    T2=_geom.Transform()
    T2.SetTrans(l2.origin+helix_trans)
    HR=_geom.Transform()
    HR.SetRot(helix_axis_rot.ToRotationMatrix())
    A1R=_geom.Transform()
    A1R.SetRot(axis1_rot.ToRotationMatrix())
    A2R=_geom.Transform()
    A2R.SetRot(axis2_rot.ToRotationMatrix())
    T=_geom.Transform()
    T.SetMatrix(T2.matrix*HR.matrix*A2R.matrix*A1R.matrix*T1.matrix)
    return (helix_trans,axis1_rot_angle,axis2_rot_angle,helix_axis_rot_angle,T)
  else:return (helix_trans,axis1_rot_angle,axis2_rot_angle,helix_axis_rot_angle)

def AnalyzeHelixRigidBodyMotion(t,helix,ref_helix,axis1,axis2=None,stride=1):
  l1=_mol.alg.structure_analysis.CalculateHelixAxis(ref_helix)
  if not axis2:axis2=_geom.Cross(l1.direction,axis1)
  tilt=(180*_geom.Angle(l1.direction,axis1)/3.1415)
  tilt=min(tilt,180-tilt)
  print 'tilt',tilt
  axis1_rot_angles=_ost.FloatList()
  axis2_rot_angles=_ost.FloatList()
  helix_axis_rot_angles=_ost.FloatList()
  translation_list=_geom.Vec3List()
  transform_list=[]
  for i in range(0,t.GetFrameCount(),stride):
    t.CopyFrame(i)
    (helix_trans,axis1_rot_angle,axis2_rot_angle,helix_axis_rot_angle,T)=FindHelixSuperpositionTransformation(ref_helix,helix,axis1,axis2,True)
    axis1_rot_angles.append(axis1_rot_angle)
    axis2_rot_angles.append(axis2_rot_angle)
    helix_axis_rot_angles.append(helix_axis_rot_angle)
    translation_list.append(helix_trans)
    transform_list.append(T)
  return translation_list,axis1_rot_angles,axis2_rot_angles,helix_axis_rot_angles,transform_list

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
  
  
  
 