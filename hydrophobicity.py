"""
.. codeauthor:: Niklaus Johner <niklaus.johner@a3.epfl.ch>

This module contains function to calculate hydrophobicity and hydrophobic
moments from different hydrophobic scales.
"""
try:
  from ost import *
  import math,sys
  import entity_alg
except:
  print 'could not import at least on of the following modules: ost,math,sys,entity_alg'

__all__=('CalculateCharge','CalculateChargeWithSlidingWindow','GetResidueHydrophobicity',\
        'CalculateHydrophobicMoment','CalculateHydrophobicMoment2','CalculateHydrophobicMoment3'\
        ,'CalculateHydrophobicMomentWithSlidingWindow','CalculateHydrophobicity'\
        ,'CalculateHydrophobicityWithSlidingWindow','DetermineAmphipathicity','AnalyzeHydrophobicMoment')

hydrophobicity_dict_kyte={'A':1.8,'C':2.5,'D':-3.5,'E':-3.5,'F':2.8,'G':-0.4,'H':-3.2,'I':4.5,'K':-3.9,'L':3.8,'M':1.9,'N':-3.5,'P':-1.6,'Q':-3.5,'R':-4.5,'S':-0.8,'T':-0.7,'V':4.2,'W':-0.9,'Y':-1.3}
hydrophobicity_dict_eisenberg={'A':0.62,'C':0.29,'D':-0.9,'E':-0.74,'F':1.19,'G':0.48,'H':-0.4,'I':1.38,'K':-1.5,'L':1.06,'M':0.64,'N':-0.78,'P':0.12,'Q':-0.85,'R':-2.53,'S':-0.18,'T':-0.05,'V':1.08,'W':0.81,'Y':0.26}

residue_charge_dict={'A':0,'C':0,'D':-1,'E':-1,'F':0,'G':0,'H':0,'I':0,'K':1,'L':0,'M':0,'N':0,'P':0,'Q':0,'R':1,'S':0,'T':0,'V':0,'W':0,'Y':0}

def CalculateCharge(sequence):
  if type(sequence)==type(seq.SequenceHandle()):sequence=sequence.GetGaplessString()
  elif type(sequence)!=type(str()):sequence=entity_alg.CreateSequenceFromView(sequence,'temp').GetGaplessString()
  c=0
  for s in sequence:
    c+=residue_charge_dict[s]
  return c
    
def CalculateChargeWithSlidingWindow(sequence,w=11):
  if type(sequence)==type(seq.SequenceHandle()):sequence=sequence.GetGaplessString()
  elif type(sequence)!=type(str()):sequence=entity_alg.CreateSequenceFromView(sequence,'temp').GetGaplessString()
  charge_list=[]
  before=int(w/2.)
  after=w-before
  n=len(sequence)
  for i in range(n):
    b=max(0,i-before)
    e=min(n,i+after)
    l=float(e-b)
    charge_list.append(CalculateCharge(sequence[b:e])/l)
  return charge_list

def GetResidueHydrophobicity(res,hydrophobicity_dict=hydrophobicity_dict_eisenberg):
  return hydrophobicity_dict(res.GetOneLetterCode())
  
def CalculateHydrophobicMoment(sequence,angle=100./180.*math.pi,hydrophobicity_dict=hydrophobicity_dict_eisenberg):
  if type(sequence)==type(seq.SequenceHandle()):sequence=sequence.GetGaplessString()
  elif type(sequence)!=type(str()):sequence=entity_alg.CreateSequenceFromView(sequence,'temp').GetGaplessString()
  v=geom.Rotate(geom.Vec2(1.0,0.0),-angle)
  moment=geom.Vec2()
  for r in sequence:
    v=geom.Rotate(v,angle)
    moment+=hydrophobicity_dict[r]*v
  return geom.Length(moment)

def CalculateHydrophobicMoment2(sequence,angle=100./180.*math.pi,hydrophobicity_dict=hydrophobicity_dict_eisenberg):
  if type(sequence)==type(seq.SequenceHandle()):sequence=sequence.GetGaplessString()
  elif type(sequence)!=type(str()):sequence=entity_alg.CreateSequenceFromView(sequence,'temp').GetGaplessString()
  v=geom.Rotate(geom.Vec2(1.0,0.0),-angle)
  moment=geom.Vec2()
  hm=CalculateHydrophobicity(sequence,hydrophobicity_dict)
  for r in sequence:
    v=geom.Rotate(v,angle)
    moment+=(hydrophobicity_dict[r]-hm)*v
  return geom.Length(moment)

def CalculateHydrophobicMoment3(sequence,angle=100./180.*math.pi,hydrophobicity_dict=hydrophobicity_dict_eisenberg):
  if type(sequence)==type(seq.SequenceHandle()):sequence=sequence.GetGaplessString()
  elif type(sequence)!=type(str()):sequence=entity_alg.CreateSequenceFromView(sequence,'temp').GetGaplessString()
  v=geom.Rotate(geom.Vec2(1.0,0.0),-angle)
  moment=geom.Vec2()
  for r in sequence:
    v=geom.Rotate(v,angle)
    moment+=hydrophobicity_dict[r]*v
  vm=moment/geom.Length(moment)
  m_max=moment
  for i in range(-3,4):
    v=geom.Rotate(geom.Vec2(1.0,0.0),-angle)
    moment=geom.Vec2()
    for r in sequence:
      v=geom.Rotate(v,angle)
      v2=v-0.2*i*vm
      v2=geom.Dot(v2,vm)/geom.Length(v2)*vm
      moment+=hydrophobicity_dict[r]*v2
    if geom.Length(moment)>geom.Length(m_max):m_max=moment
  return geom.Length(m_max)


def CalculateHydrophobicMomentWithSlidingWindow(sequence,w=11,angle=100./180.*math.pi,hydrophobicity_dict=hydrophobicity_dict_eisenberg):
  if type(sequence)==type(seq.SequenceHandle()):sequence=sequence.GetGaplessString()
  elif type(sequence)!=type(str()):sequence=entity_alg.CreateSequenceFromView(sequence,'temp').GetGaplessString()
  moment_list=[]
  before=int(w/2.)
  after=w-before
  n=len(sequence)
  for i in range(n):
    b=max(0,i-before)
    e=min(n,i+after)
    l=float(e-b)
    moment_list.append(CalculateHydrophobicMoment3(sequence[b:e],angle,hydrophobicity_dict)/l)
  return moment_list

def CalculateHydrophobicity(sequence,hydrophobicity_dict=hydrophobicity_dict_eisenberg):
  if type(sequence)==type(seq.SequenceHandle()):sequence=sequence.GetGaplessString()
  elif type(sequence)!=type(str()):sequence=entity_alg.CreateSequenceFromView(sequence,'temp').GetGaplessString()
  hydrophobicity=0.0
  for r in sequence:
    hydrophobicity+=hydrophobicity_dict[r]
  return hydrophobicity

def CalculateHydrophobicityWithSlidingWindow(sequence,w=11,hydrophobicity_dict=hydrophobicity_dict_eisenberg):
  if type(sequence)==type(seq.SequenceHandle()):sequence=sequence.GetGaplessString()
  elif type(sequence)!=type(str()):sequence=entity_alg.CreateSequenceFromView(sequence,'temp').GetGaplessString()
  hydrophobicity_list=[]
  before=int(w/2.)
  after=w-before
  n=len(sequence)
  for i in range(n):
    b=max(0,i-before)
    e=min(n,i+after)
    l=float(e-b)
    hydrophobicity_list.append(CalculateHydrophobicity(sequence[b:e],hydrophobicity_dict_eisenberg)/l)
  return hydrophobicity_list

def DetermineAmphipathicity(sequence,w=11,angle=100./180.*math.pi):
  hydrophobicity_list=CalculateHydrophobicityWithSlidingWindow(sequence,w,hydrophobicity_dict_eisenberg)
  moment_list=CalculateHydrophobicMomentWithSlidingWindow(sequence,w,angle,hydrophobicity_dict_eisenberg)
  #return [m/w for m,h in zip(moment_list,hydrophobicity_list)]
  return [m-0.6-0.382*h for m,h in zip(moment_list,hydrophobicity_list)]

def AnalyzeHydrophobicMoment(sequence,w=11,angle_min=1.39,angle_max=3.14,step=0.1,hydrophobicity_dict=hydrophobicity_dict_eisenberg):
  moment_list=[]
  for a in npy.arange(angle_min,angle_max,step):
    moment_list.append(CalculateHydrophobicMomentWithSlidingWindow(sequence,w,a,hydrophobicity_dict_eisenberg))
  return moment_list

