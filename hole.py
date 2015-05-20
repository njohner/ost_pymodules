"""
.. codeauthor:: Niklaus Johner <niklaus.johner@a3.epfl.ch>

This module contains functions to calculate the radius of a pore 
from a structure or a trajectory
"""
from ost import *
import time
import numpy as npy
import os,math
import random

__all__=("CalculatePoreRadii","AnalyzePoreRadii")

def _MinDistanceWithRadii(ev,p,maxr):
  al=ev.FindWithin(p,maxr)
  return min([geom.Distance(a.pos,p)-a.radius for a in al])


def _CalculatePoreRadius(ev,p,v1,v2,step,nsteps=1000,kT=0.1,cooling_factor=0.9):
  maxVDW=max([a.radius for a in ev.atoms])
  maxr=max(geom.Distance(a.pos,p) for a in ev.atoms)
  r=_MinDistanceWithRadii(ev,p,maxr+1.0)
  plane_center=geom.Vec3(p)
  n=geom.Cross(v1,v2)
  #print "starting with",r,p,n
  for i in range(nsteps):
    keep=False
    new_p=p+2.*step*((random.random()-0.5)*v1+(random.random()-0.5)*v2)
    new_p-=n*geom.Dot(new_p-plane_center,n)
    new_r=_MinDistanceWithRadii(ev,new_p,r+maxVDW+step)
    if new_r>r:keep=True
    elif kT>0.00001:
      if random.random()<math.exp((new_r-r)/kT):keep=True
    kT*=cooling_factor
    if keep:
      #print i,p,new_p,r,new_r,new_r-r,kT
      p=geom.Vec3(new_p)
      r=new_r
  return p,r

def CalculatePoreRadii(ev,pore_center,pore_direction,pore_lim,pore_step=0.25,mini_nsteps=1000,mini_step=0.05,kT=0.1,cooling_factor=0.9):
  #atom_pos_list=geom.Vec3List([a.pos for a in ev.atoms])
  #atom_radii=[a.radius for a in ev.atoms]
  v1=geom.OrthogonalVector(pore_direction)
  v2=geom.Cross(pore_direction,v1)
  step2=mini_step/math.sqrt(2.)
  pore_centers=geom.Vec3List()
  pore_radii=FloatList()
  p=geom.Vec3(pore_center)
  for i in npy.arange(0,pore_lim[0],-pore_step):
    new_p,r=_CalculatePoreRadius(ev,p,v1,v2,step2,mini_nsteps,kT,cooling_factor)
    p=new_p-pore_step*pore_direction
    #print i,new_p,r
    pore_centers.append(new_p)
    pore_radii.append(r)
  p=geom.Vec3(pore_center)
  for i in npy.arange(pore_step,pore_lim[1],pore_step):
    new_p,r=_CalculatePoreRadius(ev,p,v1,v2,step2)
    #print i,new_p,r
    p=new_p+pore_step*pore_direction
    pore_centers.append(new_p)
    pore_radii.append(r)
  return pore_centers,pore_radii


"""

def _MinDistanceWithRadii(atom_pos_list,atom_radii,p):
  return min([geom.Distance(ap,p)-r for ap,r in zip(atom_pos_list,atom_radii)])

def _CalculatePoreRadius(atom_pos_list,atom_radii,p,v1,v2,step,nsteps=1000,kT=0.1,cooling_factor=0.9):
  r=_MinDistanceWithRadii(atom_pos_list,atom_radii,p)
  plane_center=geom.Vec3(p)
  n=geom.Cross(v1,v2)
  #print "starting with",r,p,n
  for i in range(nsteps):
    keep=False
    new_p=p+2.*step*((random.random()-0.5)*v1+(random.random()-0.5)*v2)
    new_p-=n*geom.Dot(new_p-plane_center,n)
    new_r=_MinDistanceWithRadii(atom_pos_list,atom_radii,new_p)
    if new_r>r:keep=True
    elif kT>0.00001:
      if random.random()<math.exp((new_r-r)/kT):keep=True
    kT*=cooling_factor
    if keep:
      #print i,p,new_p,r,new_r,new_r-r,kT
      p=geom.Vec3(new_p)
      r=new_r
  return p,r

def CalculatePoreRadii(ev,pore_center,pore_direction,pore_lim,pore_step=0.25,mini_nsteps=1000,mini_step=0.05,kT=0.1,cooling_factor=0.9):
  atom_pos_list=geom.Vec3List([a.pos for a in ev.atoms])
  atom_radii=[a.radius for a in ev.atoms]
  v1=geom.OrthogonalVector(pore_direction)
  v2=geom.Cross(pore_direction,v1)
  step2=mini_step/math.sqrt(2.)
  pore_centers=geom.Vec3List()
  pore_radii=FloatList()
  p=geom.Vec3(pore_center)
  for i in npy.arange(0,pore_lim[0],-pore_step):
    new_p,r=_CalculatePoreRadius(atom_pos_list,atom_radii,p,v1,v2,step2,mini_nsteps,kT,cooling_factor)
    p=new_p-pore_step*pore_direction
    #print i,new_p,r
    pore_centers.append(new_p)
    pore_radii.append(r)
  p=geom.Vec3(pore_center)
  for i in npy.arange(pore_step,pore_lim[1],pore_step):
    new_p,r=_CalculatePoreRadius(atom_pos_list,atom_radii,p,v1,v2,step2)
    #print i,new_p,r
    p=new_p+pore_step*pore_direction
    pore_centers.append(new_p)
    pore_radii.append(r)
  return pore_centers,pore_radii

"""
def AnalyzePoreRadii(t,ev,pore_center,pore_direction,pore_lim,first=0,last=-1,stride=1,pore_step=0.25,mini_nsteps=1000,mini_step=0.05,kT=0.1,cooling_factor=0.9):
  pore_centers_list=[]
  pore_radii_list=[]
  if last==-1:last=t.GetFrameCount()
  for i in range(first,last,stride):
    print "calculating pore for frame {0}".format(i)
    t.CopyFrame(i)
    pore_centers,pore_radii=CalculatePoreRadii(ev,pore_center,pore_direction,pore_lim,pore_step,mini_nsteps,mini_step,kT,cooling_factor)
    pore_centers_list.append(pore_centers)
    pore_radii_list.append(pore_radii)
    pore_center=pore_centers[0]
  nframes=float(len(pore_radii_list))
  nsteps=len(pore_centers)
  average_centers=sum(pore_centers_list)/nframes
  average_radii=[sum([el[i] for el in pore_radii_list])/nframes for i in range(nsteps)]
  return pore_centers_list,pore_radii_list,average_centers,average_radii

def CreatePoreEntity(pore_centers,pore_radii):
  pore=mol.CreateEntity()
  edi=pore.EditXCS(mol.BUFFERED_EDIT)
  c=edi.InsertChain("P")
  for i,(pos,radius) in enumerate(zip(pore_centers,pore_radii)):
    if i%10==0:r=edi.AppendResidue(c,"POR")
    a=edi.InsertAtom(r,"C",pos)
    a.SetRadius(radius)
  return pore




