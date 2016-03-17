#------------------------------------------------------------------------------------------------
#This file is part of the ost_pymodules project (https://github.com/njohner/ost_pymodules).
#
#Copyright 2015 Niklaus Johner
#
#ost_pymodules is free software: you can redistribute it and/or modify
#it under the terms of the GNU Lesser General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#ost_pymodules is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Lesser General Public License for more details.
#
#You should have received a copy of the GNU Lesser General Public License
#along with ost_pymodules.  If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------------------------
"""

.. codeauthor:: Niklaus Johner <niklaus.johner@a3.epfl.ch>

This module contains basic functions to work with sequences, notably to find motifs,
and build a position specific scoring matrix from an alignment.
"""

__all__=('FindMotif','CreateSequenceFromView')

import ost as _ost
import random

def CreateSequenceFromView(eh,seq_name):
  """
  Returns the sequence of the view *eh*, named *seq_name*.
  """
  s=''
  for r in eh.residues:
    s+=conop.ResidueNameToOneLetterCode(r.name)
  return seq.CreateSequence(seq_name,s)

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

class PSSM:
  """
  Class to create a position specific scoring matrix from a multiple sequence alignment.
  The object can then be used to score a sequence.
  """
  def __init__(self,ali):#,aa_list=None,pseudocounts=None,aa_frequencies=None):
    self.aa_list=['D','E','R','K','H','S','T','N','Q','G','P','Y','W','V','I','L','M','F','A','C']
    self.aa_index_dict={}
    for i,aa in enumerate(self.aa_list):
      self.aa_index_dict[aa]=i
    self.aa_frequencies={"A": 0.09398206173005932, "C": 0.014640847892538094, "E": 0.06167898538341536, "D": 0.057537026726662464, "G": 0.07661345622483157, "F": 0.04125062078013584, "I": 0.05988964354148949, "H": 0.023180301946789893, "K": 0.048775277841393995, "M": 0.022045733851907126, "L": 0.09741587583495241, "N": 0.03822084621874947, "Q": 0.03472639486143553, "P": 0.043952778781438454, "S": 0.060691166121181965, "R": 0.05494787451396235, "T": 0.05250697453686419, "W": 0.012220035049334168, "V": 0.07323789898435108, "Y": 0.032486199178507216}
    self.ncol=ali.GetLength()
    self.ali_depth=ali.GetCount()
    nres=len(self.aa_list)
    self.matrix=npy.zeros([self.ncol,nres])
    print "generating pssm from alignment with {0} col and {1} sequences".format(self.ncol,self.ali_depth)
    for i,aa in enumerate(self.aa_list):
      self.matrix[:,i]=nres*self.aa_frequencies[aa]
    for i,col in enumerate(ali):
      for j,aa in enumerate(self.aa_list):
        self.matrix[i,j]+=str(col).count(aa)
      self.matrix[i,:]/=float(npy.sum(self.matrix[i,:]))
    for i,col in enumerate(ali):
      for j,aa in enumerate(self.aa_list):
        #print self.matrix[i,j]/self.aa_frequencies[aa]
        self.matrix[i,j]=math.log(self.matrix[i,j]/self.aa_frequencies[aa])

  def ScoreSequence(self,s,start_index):
    score=0
    for i,c in enumerate(s):
      if not c in self.aa_index_dict:continue
      score+=self.matrix[i+start_index,self.aa_index_dict[c]]
    return score




