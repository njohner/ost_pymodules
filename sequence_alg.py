"""
Module written by Niklaus Johner (niklaus.johner@a3.epfl.ch) 01.2013

This module contains basic functions to work with sequences
"""

__all__=('FindMotif')

import ost as _ost
import random

def FindMotif(motif,sequence):
  s=sequence.GetGaplessString()
  start_list=[]
  for i in range(0,len(s)-len(motif)):
    start=False
    for el in motif[0]:
      if s[i]==el:
        start=True
        break
    if not start:continue
    for j in range(1,len(motif)):
      if 'X' in motif[j]:
        flag=True
        continue
      flag=False
      for el in motif[j]:
        if s[i+j]==el:
          flag=True
          break
      if not flag:break
    if flag:start_list.append(i)
  return start_list
  

def FindOneOfSeveralMotifs(motif_list,sequence):
  for i,motif in enumerate(motif_list):
    start_list=FindMotif(motif,sequence)
    if len(start_list)>0:return start_list,i
  return None

def ShuffleSequence(sequence,fix_gaps=True):
  if not fix_gaps:
    s=[el for el in sequence]
    random.shuffle(s)
    s="".join(s)
    return seq.CreateSequence(sequence.name,s)
  non_gap_indices=[]
  s=[el for el in sequence if not el in ["-","?"]]
  random.shuffle(s)
  s2=""
  c=0
  for i,el in enumerate(sequence):
    if not el in ["-","?"]:
      s2+=s[c]
      c+=1
    else:s2+="-"
  return seq.CreateSequence(sequence.name,s2)


def RandomizeAlignment(ali,fix_gaps=True):
  ali2=seq.CreateAlignment()
  for s in ali.sequences:
    s2=ShuffleSequence(s,fix_gaps)
    ali2.AddSequence(s2)
  return ali2


