from ost import *

"""
This module is used to compute hydrogen bonds from entities and trajectories
Criteria used are similar to the HBPlus program
"""

class HBondableAtoms:
	def __init__(self,donors=[],acceptors=[]):
		self.donors=donors
		self.acceptors=acceptors
    
def BuildCHARMMHBondDonorAcceptorDict():
  hb_da_dict={}
  bb_donors=[['N','HN']]
  bb_acceptors=[['O',['C']]]
  
  #hydrophobic
  hb_da_dict['ALA']=HBondableAtoms(bb_donors,bb_acceptors)
  hb_da_dict['GLY']=HBondableAtoms(bb_donors,bb_acceptors)
  hb_da_dict['ILE']=HBondableAtoms(bb_donors,bb_acceptors)
  hb_da_dict['LEU']=HBondableAtoms(bb_donors,bb_acceptors)
  hb_da_dict['PHE']=HBondableAtoms(bb_donors,bb_acceptors)
  hb_da_dict['MET']=HBondableAtoms(bb_donors,bb_acceptors)
  hb_da_dict['VAL']=HBondableAtoms(bb_donors,bb_acceptors)
  #special case
  hb_da_dict['PRO']=HBondableAtoms([],bb_acceptors)
  
  #sidechain donors
  ne=[['NE','HE']]
  cz=[['NH1','HH11'],['NH1','HH12'],['NH2','HH21'],['NH2','HH22']]
  hb_da_dict['ARG']=HBondableAtoms(bb_donors+ne+cz,bb_acceptors)
  nz=[['NZ','HZ1'],['NZ','HZ2'],['NZ','HZ3']]
  hb_da_dict['LYS']=HBondableAtoms(bb_donors+nz,bb_acceptors)
  ne1=[['NE1','HE1']]
  hb_da_dict['TRP']=HBondableAtoms(bb_donors+ne1,bb_acceptors)
  sg=[['SG','HG1']]
  hb_da_dict['CYS']=HBondableAtoms(bb_donors+sg,bb_acceptors)
  
  #sidechain acceptors
  od12=[['OD1',['CG']],['OD2',['CG']]]
  hb_da_dict['ASP']=HBondableAtoms(bb_donors,bb_acceptors+od12)
  oe12=[['OE1',['CD']],['OE2',['CD']]]
  hb_da_dict['GLU']=HBondableAtoms(bb_donors,bb_acceptors+oe12)
  
  #sidechain donor and acceptor
  od1=[['OD1',['CG']]]
  nd2=[['ND2','HD21'],['ND2','HD22']]
  hb_da_dict['ASN']=HBondableAtoms(bb_donors+nd2,bb_acceptors+od1)
  ne2=[['NE2','HE21'],['NE2','HE22']]
  oe1=[['OE1',['CD']]]
  hb_da_dict['GLN']=HBondableAtoms(bb_donors+ne2,bb_acceptors+oe1)
  og_d=[['OG','HG1']]
  og_a=[['OG',['CB']]]
  hb_da_dict['SER']=HBondableAtoms(bb_donors+og_d,bb_acceptors+og_a)
  og1_d=[['OG1','HG1']]
  og1_a=[['OG1',['CB']]]
  hb_da_dict['THR']=HBondableAtoms(bb_donors+og1_d,bb_acceptors+og1_a)
  oh_d=[['OH','HH']]
  oh_a=[['OH',['CZ']]]
  hb_da_dict['TYR']=HBondableAtoms(bb_donors+oh_d,bb_acceptors+oh_a)
  #histidine
  nd1_d=[['ND1','HD1']]
  ne2_a=[['NE2',['CD2','CE1']]]
  hb_da_dict['HSD']=HBondableAtoms(bb_donors+nd1_d,bb_acceptors+ne2_a)
  ne2_d=[['NE2','HE2']]
  nd1_a=[['ND1',['CG','CE1']]]
  hb_da_dict['HSE']=HBondableAtoms(bb_donors+ne2_d,bb_acceptors+nd1_a)
  hb_da_dict['HSP']=HBondableAtoms(bb_donors+nd1_d+ne2_d,bb_acceptors)
  #non-standard protonation state:
  oe12=[['OE1',['CD']],['OE2',['CD']]]
  oe2=[['OE2','HE2']]
  hb_da_dict['GLUP']=HBondableAtoms(bb_donors+oe2,bb_acceptors+oe12)
  od12=[['OD1',['CG']],['OD2',['CG']]]
  od2=[['OD2','HD2']]
  hb_da_dict['ASPP']=HBondableAtoms(bb_donors+od2,bb_acceptors+od12)
  return hb_da_dict

def BuildCHARMMHBondDonorEquivalenceDict():
  donor_swap_dict={}
  donor_swap_dict['ARG']=[[['NH1','HH11'],['NH1','HH12'],['NH2','HH21'],['NH2','HH22']]]
  donor_swap_dict['ASN']=[[['ND2','HD21'],['ND2','HD22']]]
  donor_swap_dict['GLN']=[[['NE2','HE21'],['NE2','HE22']]]
  donor_swap_dict['LYS']=[[['NZ','HZ1'],['NZ','HZ2'],['NZ','HZ3']]]
  return donor_swap_dict

def BuildCHARMMHBondAcceptorEquivalenceDict():
  acceptor_swap_dict={}
  acceptor_swap_dict['ASP']=[[['OD1',['CG']],['OD2',['CG']]]]
  acceptor_swap_dict['GLU']=[[['OE1',['CD']],['OE2',['CD']]]]
  return acceptor_swap_dict

def ListEquivalentDonors(donor,donor_swap_dict):
  if not donor.heavy_atom.residue.name in donor_swap_dict:return [donor]
  donor_list=[donor]
  donor_atom_names=[donor.heavy_atom.name,donor.hydrogen.name]
  res=donor.heavy_atom.residue
  for equivalence_list in donor_swap_dict[donor.heavy_atom.residue.name]:
    if not donor_atom_names in equivalence_list:continue
    for atom_names in equivalence_list:
      if atom_names==donor_atom_names:continue
      else:donor_list.append(HBondDonor.FromResidue(res,atom_names[0],atom_names[1]))
  return donor_list

def ListEquivalentAcceptors(acceptor,acceptor_swap_dict):
  if not acceptor.atom.residue.name in acceptor_swap_dict:return [acceptor]
  acceptor_list=[acceptor]
  acceptor_atom_names=[acceptor.atom.name,[a.name for a in acceptor.antecedent_list]]
  res=acceptor.atom.residue
  for equivalence_list in acceptor_swap_dict[acceptor.atom.residue.name]:
    if not acceptor_atom_names in equivalence_list:continue
    for atom_names in equivalence_list:
      if atom_names==acceptor_atom_names:continue
      else:acceptor_list.append(HBondAcceptor.FromResidue(res,atom_names[0],atom_names[1]))
  return acceptor_list

  
class HBondDonor:
  def __init__(self,donor,hydrogen):
    self.heavy_atom=donor
    self.hydrogen=hydrogen
  def __eq__(self,hbd):
    if not isinstance(hbd,HBondDonor):return False
    else:return self.heavy_atom==hbd.heavy_atom and self.hydrogen==hbd.hydrogen
  def __hash__(self):
    return hash((self.heavy_atom,self.hydrogen))
  @classmethod
  def FromResidue(cls,res,donor_name,hydrogen_name,verbose=True):
    _donor=res.FindAtom(donor_name).handle
    _hydrogen=res.FindAtom(hydrogen_name).handle
    if not _donor.IsValid():
      if verbose:print 'Could not find '+donor_name+' in residue '+str(res)
      return
    if not _hydrogen.IsValid():
      if verbose:print 'Could not find '+hydrogen_name+' in residue '+str(res)
      return
    return cls(_donor,_hydrogen)
  
  def IsHBondedTo(self,acceptor):
    return AreHBonded(self,acceptor)


class HBondAcceptor:
  def __init__(self,acceptor,antecedent_list):
    self.atom=acceptor
    self.antecedent_list=antecedent_list
  def __eq__(self,hba):
    if not isinstance(hba,HBondAcceptor):return False
    else:return self.atom==hba.atom
  def __hash__(self):
    return hash(self.atom)
  @classmethod
  def FromResidue(cls,res,acceptor_name,antecedent_name_list,verbose=True):
    _acceptor=res.FindAtom(acceptor_name).handle
    _antecedent_list=mol.AtomHandleList([res.FindAtom(name).handle for name in antecedent_name_list])
    if not _acceptor.IsValid():
      if verbose:print 'Could not find '+acceptor_name+' in residue '+str(res)
      return
    for i,a in enumerate(_antecedent_list):
      if not a.IsValid():
        if verbose:print 'Could not find '+antecedent_name_list[i]+' in residue '+str(res)
        return
    return cls(_acceptor,_antecedent_list)
  
  def IsHBondedTo(self,donor):
    return AreHBonded(donor,self)

class HBond:
  def __init__(self,donor,acceptor):
    self.donor=donor
    self.acceptor=acceptor
  def __eq__(self,hb):
    if not isinstance(hb,HBond):return False
    else:return self.acceptor==hb.acceptor and self.donor==hb.donor
  def __hash__(self):
    return hash((self.donor,self.acceptor))
  def IsFormed(self):return AreHBonded(self.donor,self.acceptor)
  def DonorAcceptorDistance(self):return geom.Distance(self.donor.heavy_atom.pos,self.acceptor.atom.pos)
  def HydrogenAcceptorDistance(self):return geom.Distance(self.donor.hydrogen.pos,self.acceptor.atom.pos)
  def DonorHydrogenAcceptorAngle(self):
    dp=self.donor.heavy_atom.pos;ap=self.acceptor.atom.pos;hp=self.donor.hydrogen.pos
    return geom.Angle(dp-hp,ap-hp)
  def DonorAcceptorAcceptorAntecedentAngle(self):
    dp=self.donor.heavy_atom.pos;ap=self.acceptor.atom.pos
    a=min([geom.Angle(aa.pos-ap,dp-ap) for aa in self.acceptor.antecedent_list])
    return a
  def HydrogenAcceptorAcceptorAntecedentAngle(self):
    hp=self.donor.hydrogen.pos;ap=self.acceptor.atom.pos
    a=min([geom.Angle(aa.pos-ap,hp-ap) for aa in acceptor.antecedent_list])
    return a

    
def AreHBonded(donor,acceptor,da_dist=3.9,ha_dist=2.5,dha_angle=1.57,daaa_angle=1.57,haaa_angle=1.57):
  """
  determines if a donor/acceptor pair is hydrogen bonded or not
  """
  dp=donor.heavy_atom.pos
  ap=acceptor.atom.pos
  if geom.Distance(dp,ap)>da_dist:return False
  hp=donor.hydrogen.pos
  if geom.Distance(hp,ap)>ha_dist:return False
  if geom.Angle(dp-hp,ap-hp)<dha_angle:return False
  a=min([geom.Angle(aa.pos-ap,dp-ap) for aa in acceptor.antecedent_list])
  if a<daaa_angle:return False
  a=min([geom.Angle(aa.pos-ap,hp-ap) for aa in acceptor.antecedent_list])
  if a<haaa_angle:return False
  return True

def GetHbondDonorAcceptorList(eh,hbond_donor_acceptor_dict={},verbose=True):
  """
  returns a list of hydrogen-bond donors and acceptors from an Entity or EntityView.
  It relies on atom names to determine the list of H-bond donors and acceptors.
  These names are given in a dictionary, which defaults to CHARMM.
  """
  if not hbond_donor_acceptor_dict:
    print 'Using default CHARMM atom namings'
    hbond_donor_acceptor_dict=BuildCHARMMHBondDonorAcceptorDict()
  donor_list=[]
  acceptor_list=[]
  for r in eh.residues:
    if not r.name in hbond_donor_acceptor_dict:
      print 'donors and acceptors for',r,'are not defined in the dictionary and will not be included'
      continue
    res_da_dict=hbond_donor_acceptor_dict[r.name]
    for acceptor in res_da_dict.acceptors:
      a=HBondAcceptor.FromResidue(r,acceptor[0],acceptor[1],verbose)
      if not a==None:acceptor_list.append(a)
    for donor in res_da_dict.donors:
      d=HBondDonor.FromResidue(r,donor[0],donor[1],verbose)
      if not d==None:donor_list.append(d)
  return [donor_list,acceptor_list]  
  
  
def GetHbondListFromDonorAcceptorLists(donor_list,acceptor_list):
  """
  return a list of hydrogen bonds between donors and acceptors from
  a list of donors and a list of acceptors.
  """
  hbond_list=[]
  for donor in donor_list:
    for acceptor in acceptor_list:
      if AreHBonded(donor,acceptor):hbond_list.append(HBond(donor,acceptor))
  return hbond_list

def GetHbondListFromView(eh,hbond_donor_acceptor_dict={},verbose=True):
  """
  return a list of hydrogen bonds from an Entity or EntityView.
  if no dictionary for the hbond donors and acceptors is specified
  it will use the standard CHARMM names to determine them
  """
  [donor_list,acceptor_list]=GetHbondDonorAcceptorList(eh,hbond_donor_acceptor_dict,verbose)
  hbond_list=[]
  for donor in donor_list:
    for acceptor in acceptor_list:
      if AreHBonded(donor,acceptor):hbond_list.append(HBond(donor,acceptor))
  return hbond_list

def GetHbondListFromTraj(t,eh,cutoff=0.7,stride=1,swap=False,donor_swap_dict={},acceptor_swap_dict={},hbond_donor_acceptor_dict={},verbose=True):
  """
  return a list of hydrogen bonds from an Entity or EntityView.
  if no dictionary for the hbond donors and acceptors is specified
  it will use the standard CHARMM names to determine them
  """
  if swap:
    if not donor_swap_dict:
      print 'use of standard CHARMM HBond donor swap dictionary'
      donor_swap_dict=BuildCHARMMHBondDonorEquivalenceDict()
    if not acceptor_swap_dict:
      print 'use of standard CHARMM HBond acceptor swap dictionary'  
      acceptor_swap_dict=BuildCHARMMHBondAcceptorEquivalenceDict()
  [donor_list,acceptor_list]=GetHbondDonorAcceptorList(eh,hbond_donor_acceptor_dict,verbose)
  hb_list=[]
  hb_score=[]
  nframes=0
  for i in range(0,t.GetFrameCount(),stride):
    t.CopyFrame(i)
    nframes+=1
    hbond_list=GetHbondListFromDonorAcceptorLists(donor_list,acceptor_list)
    if swap:
      hbond_list=GetEquivalentHBonds(hbond_list,eh,swap,donor_swap_dict,acceptor_swap_dict,verbose)
      for hb in hbond_list:
        if hb in hbond_list[hbond_list.index(hb)+1:]:hbond_list.pop(hbond_list.index(hb))
    for hb in hbond_list:
      if hb in hb_list:hb_score[hb_list.index(hb)]+=1
      else:
        hb_score.append(1)
        hb_list.append(hb)
  hbond_list=[]
  hbond_score=[]
  for hb,s in zip(hb_list,hb_score):
    if s/float(nframes)>=cutoff:
      if swap:hbond_list.append(list(hb)[0])
      else:hbond_list.append(hb)
      hbond_score.append(s/float(nframes))
  return (hbond_list,hbond_score)

def GetHbondListBetweenViews(eh1,eh2,hbond_donor_acceptor_dict={},verbose=True):
  """
  return the list of hydrogen bonds formed between two Entity or EntityView.
  if no dictionary for the hbond donors and acceptors is specified
  it will use the standard CHARMM names to determine them
  """
  [donor_list1,acceptor_list1]=GetHbondDonorAcceptorList(eh1,hbond_donor_acceptor_dict,verbose)
  [donor_list2,acceptor_list2]=GetHbondDonorAcceptorList(eh2,hbond_donor_acceptor_dict,verbose)
  hbond_list=[]
  for donor in donor_list1:
    for acceptor in acceptor_list2:
      if AreHBonded(donor,acceptor):
        hb=HBond(donor,acceptor)
        if hb in hbond_list:continue
        hbond_list.append(hb)
  for donor in donor_list2:
    for acceptor in acceptor_list1:
      if AreHBonded(donor,acceptor):
        hb=HBond(donor,acceptor)
        if hb in hbond_list:continue
        hbond_list.append(hb)
  return hbond_list

def GetEquivalentHBonds(ref_hbond_list,eh,swap=False,donor_swap_dict={},acceptor_swap_dict={},verbose=True):
  hbond_list=[]
  if not swap:
    for hbond in ref_hbond_list:
      res1=eh.FindResidue(hbond.donor.heavy_atom.chain.name,hbond.donor.heavy_atom.residue.number.num)
      res2=eh.FindResidue(hbond.acceptor.atom.chain.name,hbond.acceptor.atom.residue.number.num)
      donor=HBondDonor.FromResidue(res1,hbond.donor.heavy_atom.name,hbond.donor.hydrogen.name,verbose)
      acceptor=HBondAcceptor.FromResidue(res2,hbond.acceptor.atom.name,[a.name for a in hbond.acceptor.antecedent_list],verbose)
      hbond_list.append(set([HBond(donor,acceptor)]))
  else:
    if not donor_swap_dict:
      print 'use of standard CHARMM HBond donor swap dictionary'
      donor_swap_dict=BuildCHARMMHBondDonorEquivalenceDict()
    if not acceptor_swap_dict:
      print 'use of standard CHARMM HBond acceptor swap dictionary'  
      acceptor_swap_dict=BuildCHARMMHBondAcceptorEquivalenceDict()
    for hbond in ref_hbond_list:
      res1=eh.FindResidue(hbond.donor.heavy_atom.chain.name,hbond.donor.heavy_atom.residue.number.num)
      res2=eh.FindResidue(hbond.acceptor.atom.chain.name,hbond.acceptor.atom.residue.number.num)
      donor=HBondDonor.FromResidue(res1,hbond.donor.heavy_atom.name,hbond.donor.hydrogen.name,verbose)
      acceptor=HBondAcceptor.FromResidue(res2,hbond.acceptor.atom.name,[a.name for a in hbond.acceptor.antecedent_list],verbose)
      donor_list=ListEquivalentDonors(donor,donor_swap_dict)
      acceptor_list=ListEquivalentAcceptors(acceptor,acceptor_swap_dict)
      hbond_list.append(set([HBond(d,a) for d in donor_list for a in acceptor_list]))
  return hbond_list
    
    
def CalculateHBondScore(ref_eh,eh2,ref_eh2=None,hbond_donor_acceptor_dict={},swap=False,donor_swap_dict={},acceptor_swap_dict={},verbose=True):
  """
  Returns the fraction of H-bonds from ref_eh that are also present in eh2.
  If ref_eh2 is specified, it uses as reference the Hbonds between ref_eh and ref_eh2. 
  This allows to look at H-bonds between specific parts of proteins or so.
  Alternatively ref_eh can be a list of H-bonds.
  This function relies on atom names to determine the list of H-bond donors and acceptors.
  These names are given in a dictionary, which defaults to CHARMM.
  If swap is set to True, a dictionary for equivalent donors and one for equivalent acceptors
  (defaults to CHARMM) is used to check for equivalent HBonds in eh2.
  If swap is set to True, if two equivalent hydrogen bonds are present in the reference entity 
  (for example both oxygens of ASP H-bonding the same atom), it suffices that on of these bonds is
  present in eh2 for both of them to be counted as present in eh2.
  """
  if ref_eh2:hbond_list1=GetHbondListBetweenViews(ref_eh,ref_eh2,hbond_donor_acceptor_dict,verbose)
  elif type(ref_eh)==list:hbond_list1=ref_eh
  else:hbond_list1=GetHbondListFromView(ref_eh,hbond_donor_acceptor_dict,verbose)
  nbonds=float(len(hbond_list1))
  if nbonds==0:
    print 'No HBonds in reference view'
    return None
  hbond_list2=GetEquivalentHBonds(hbond_list1,eh2,swap,donor_swap_dict,acceptor_swap_dict,verbose)
  c=0
  for hl in hbond_list2:
    for hbond in hl:
      if HBond(hbond.donor,hbond.acceptor).IsFormed():
        c+=1
        break
  return c/float(len(hbond_list1))


def AnalyzeHBondScore(ref_eh,t,eh2,ref_eh2=None,hbond_donor_acceptor_dict={},swap=False,donor_swap_dict={},acceptor_swap_dict={},first=0,last=-1,stride=1,verbose=True):
  """
  Returns the same score as CalculateHBondScore, but for a trajectory.
  """
  if swap:
    if not donor_swap_dict:
      print 'use of standard CHARMM HBond donor swap dictionary'
      donor_swap_dict=BuildCHARMMHBondDonorEquivalenceDict()
    if not acceptor_swap_dict:
      print 'use of standard CHARMM HBond acceptor swap dictionary'  
      acceptor_swap_dict=BuildCHARMMHBondAcceptorEquivalenceDict()
  if ref_eh2:hbond_list1=GetHbondListBetweenViews(ref_eh,ref_eh2,hbond_donor_acceptor_dict,verbose)
  elif type(ref_eh)==list:hbond_list1=ref_eh
  else:hbond_list1=GetHbondListFromView(ref_eh,hbond_donor_acceptor_dict,verbose)
  nbonds=float(len(hbond_list1))
  if nbonds==0:
    print 'No HBonds in reference view'
    return None
  print 'number of hbonds in ref_eh:',nbonds
  hbond_list2=GetEquivalentHBonds(hbond_list1,eh2,swap,donor_swap_dict,acceptor_swap_dict,verbose)
  if last==-1:last=t.GetFrameCount()
  score=FloatList()
  for f in range(first,last,stride):
    t.CopyFrame(f)
    c=0
    for hl in hbond_list2:
      for hbond in hl:
        if hbond.IsFormed():
          c+=1
          break
    score.append(c/nbonds)
  return score


def GetHBondListIntersection(ref_hbond_list,ref_eh,hbond_list,swap=False,donor_swap_dict={},acceptor_swap_dict={}):
  if swap:
    if not donor_swap_dict:
      print 'use of standard CHARMM HBond donor swap dictionary'
      donor_swap_dict=BuildCHARMMHBondDonorEquivalenceDict()
    if not acceptor_swap_dict:
      print 'use of standard CHARMM HBond acceptor swap dictionary'  
      acceptor_swap_dict=BuildCHARMMHBondAcceptorEquivalenceDict()
  hbond_list=GetEquivalentHBonds(hbond_list,ref_eh,swap,donor_swap_dict,acceptor_swap_dict)
  ref_hbond_list=GetEquivalentHBonds(ref_hbond_list,ref_eh,swap,donor_swap_dict,acceptor_swap_dict)
  out_hbond_list=[]
  for hb in ref_hbond_list:
    if hb in hbond_list:
      out_hbond_list.append(hb)
  return out_hbond_list

  
  
  
  
  
