
import os
import sys
import json
import getopt

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolTransforms
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from numpy  import *
import numpy as np
from numpy import (array, dot, arccos, clip,zeros,degrees)

from itertools import *


def iterAngles(mol):
    """ Return an iterator over all angles in molecule """
    natoms = mol.GetNumAtoms() # Get the coordinates of all atoms
    conf = mol.GetConformer()
    coords = [conf.GetAtomPosition(i) for i in range(natoms)]
    for atom in mol.GetAtoms():
        center = coords[atom.GetIdx()]
        neighbours = [a.GetIdx() for a in atom.GetNeighbors()]
        for a1, a2 in combinations(neighbours, 2):
            yield coords[a1], center, coords[a2]


def cangle4(v1, v2):
    v12 = np.dot(v1, v2)
    v11v22 = np.dot(v1, v1) * np.dot(v2, v2)
    return np.sign(v12) * (v12 * v12) / v11v22

def angleContribution4(mol):
    """ Return the angle contributions """
    a60 = a90 = a102 = aUnmatched = 0
    cos102 = np.cos(102.0 / 180 * np.pi)
    cos102 = -cos102 * cos102
    for a1, c, a2 in iterAngles(mol):
        v1 = np.array([a1.x - c.x, a1.y - c.y, a1.z - c.z])
        v2 = np.array([a2.x - c.x, a2.y - c.y, a2.z - c.z])
        cosAngle = cangle4(v1, v2)
        if 0.5 * 0.5 <= cosAngle:
            a60 += 1
        elif 0 <= cosAngle < 0.5 * 0.5:
            a90 += 1
        elif cos102 <= cosAngle < 0:
            a102 += 1
        else:
            aUnmatched += 1
    return (a60, a90, a102)


def findHBonds(m,confId=-1,possiblePartners='[#8,#7,#9]',possibleHs='[#1][#8,#7,#16]',distThresh=2.5):
    conf = m.GetConformer(confId)
    partners =[x[0] for x in  m.GetSubstructMatches(Chem.MolFromSmarts(possiblePartners))]
    hs=  [x[0] for x in m.GetSubstructMatches(Chem.MolFromSmarts(possibleHs))]
    res = []
    for h in hs:
        ph = conf.GetAtomPosition(h)
        for partner in partners:
            if m.GetBondBetweenAtoms(h,partner) is not None:
                continue
            d = conf.GetAtomPosition(partner).Distance(ph)
            if d<=distThresh:
                res.append((h,partner,d))
    return res

def findintramolecularHH(m,confId=-1,possiblePartners='[#1][*]',possibleHs='[#1]',distThresh=2.5):
    conf = m.GetConformer(confId)
    partners =[x[0] for x in  m.GetSubstructMatches(Chem.MolFromSmarts(possiblePartners))]
    linkers =[x[1] for x in  m.GetSubstructMatches(Chem.MolFromSmarts(possiblePartners))]
    hs=  [x[0] for x in m.GetSubstructMatches(Chem.MolFromSmarts(possibleHs))]
    res = []
    for h in hs:
        ph = conf.GetAtomPosition(h)
        for partner,linker in zip(partners,linkers):
            if m.GetBondBetweenAtoms(h,linker) is not None:
                continue
            d = conf.GetAtomPosition(partner).Distance(ph)
            if d<=distThresh:
                res.append((h,partner,d))
    return res


def angle_HH(m):

    result = {};

    mol=Chem.AddHs(m)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)
    h1=len(findHBonds(mol,distThresh=1.74))  #Donor_Acceptor(mol) # need to check that parameter!

    #Dist=HHdistances(mol)
    h2=len(findintramolecularHH(mol,distThresh=1.999)) #Dist[0]
    h3=len(findintramolecularHH(mol,distThresh=2.3))-h2 #Dist[1]

    # angle issue: need to use the Bond angle not the Torsions!
    BA=angleContribution4(mol)
    #TA=AnglesTorsion(mol)
    result["AngleContrib_1"] = BA[0];
    result["AngleContrib_2"] = BA[1];
    result["AngleContrib_3"] = BA[2];
    result["Angle_H1"] = h1;
    result["Angle_H2"] = h2;
    result["Angle_H3"] = h3;

    return result;

##### RDKit descriptors listing 
_rdkitdescDict = dict(Descriptors.descList)
rdkitdescs = _rdkitdescDict.keys()

def rdkitcalcdesc(rdkitmol, descnames=[]):
        """Calculate descriptor values.

        Optional parameter:
           descnames -- a list of names of descriptors

        If descnames is not specified, all available descriptors are
        calculated. See the descs variable for a list of available
        descriptors.
        """
        if not descnames:
            descnames = rdkitdescs
        ans = {}
        for descname in descnames:
            try:
                desc = _rdkitdescDict[descname]
            except KeyError:
                raise ValueError ("%s is not a recognised RDKit descriptor type" % descname)
            ans[descname] = desc(rdkitmol)
        return ans

####### polarizability equation 
def add_pol_vector(mol):
  NB_fonction =13
  fonction_name = []
  fonction_smarts = []
  fonction_name.append('CH3 sp3 1')
  fonction_smarts.append('[CX4H3]')
  fonction_name.append('CH2 sp3 2')
  fonction_smarts.append('[CX4H2]')
  fonction_name.append('CH sp3 3 ')
  fonction_smarts.append('[CX4H1]')
  fonction_name.append('F desc 5')
  fonction_smarts.append('F')
  fonction_name.append('Cl desc 6')
  fonction_smarts.append('Cl')
  fonction_name.append('Br desc 7')
  fonction_smarts.append('Br')
  fonction_name.append('I desc 8')
  fonction_smarts.append('I')
  fonction_name.append('NITRO 9')
  fonction_smarts.append('[$([NX3](=O)=O),$([NX3+](=O)[O-])]')
  fonction_name.append('N desc 10')
  fonction_smarts.append('[N,n]')
  fonction_name.append('O desc 11')
  fonction_smarts.append('[O,o]')
  fonction_name.append('S sulfone 12')
  fonction_smarts.append('[$([SX4](=[OX1])(=[OX1])([#6])[#6]),$([SX4+2]([OX1-])([OX1-])([#6])[#6])]')
  fonction_name.append('soufre any 13')
  fonction_smarts.append('[S,s]')
  fonction_name.append('phosphore any 14')
  fonction_smarts.append('[P,p]')

  
  desc = zeros(13)
  i=0
  while i < NB_fonction:
    #print str(i+1)
    func = Chem.MolFromSmarts(fonction_smarts[i])
    maps = mol.GetSubstructMatches(func)
    k=len(maps)
    desc[i]=k
    i +=1
  return desc

def countAtom(mol,index):
  k=0
  for atom in mol.GetAtoms():
     if atom.GetAtomicNum()==index:
        k += 1
  return k


def add_V_E(mol):
  #M = obabelreader(smi)
  #d = obcalcdesc(M,["MR"])
  MolMR = _rdkitdescDict["MolMR"](mol)
  mol.SetProp("MolMR", str(MolMR))
  R = _rdkitdescDict["RingCount"](mol)
  mol.SetProp("RingCount", str(R))

  # for key, value in dict.items(d):
  #   v = str(value).replace('nan', '0')
  #   mol.SetProp(key, v)
  #   if key == "RingCount":
  #     R = value
  #   if key =="MolMR":
  #     MolMR = value
  #  enumarate the atomic Elements caution not working with all the elements!!!!
  nC = countAtom(mol,6)
  nB = countAtom(mol,5)
  nO = countAtom(mol,8)
  nN = countAtom(mol,7)
  nS = countAtom(mol,16)
  nI= countAtom(mol,53)
  nCl = countAtom(mol,17)
  nBr = countAtom(mol,35)

  k=0
  k2=0
  m3 = Chem.RemoveHs(mol)
  k= m3.GetNumAtoms()
  m2 = Chem.AddHs(m3)
  k2= m2.GetNumAtoms()
  nH = k2-k

  nP = countAtom(mol,15)
  nF = countAtom(mol,9)
  nTe = countAtom(mol,52)
  nSn = countAtom(mol,50)
  nSb = countAtom(mol,51)
  nSe = countAtom(mol,34)
  nAs = countAtom(mol,33)
  nGe = countAtom(mol,32)
  nSi = countAtom(mol,14)

  # definition: McGowan's Characteristic volume (Vx) p 872 volume descriptors
  # cavity term in linaer solvation energy relationship
  Sumx = 0.0871*nH +0.1635*nC + 0.1439*nN + nO*0.1243 + nP*0.2487 + nS*0.2291 + nF*0.1048 + nCl*0.2095 + nBr*0.2621 + nI*0.3453
  Sumx = Sumx + 0.2683*nSi + 0.1832*nB + 0.3102*nGe +0.2942*nAs + 0.2781*nSe + 0.3935*nSn + 0.3774*nSb + 0.3614*nTe
  
  # number of bonds is defined by this equation:
  B = nH + nC + nN + nO + nP + nS + nF + nCl + nBr + nI + nSi + nB + nGe + nAs + nSe + nSn + nSb + nTe - 1 + R

  # McGowan Volume definition caution already "/ 100"
  Vx = Sumx-0.0656*B

  desc = add_pol_vector(mol)
  # remove the sulfones in nS
  desc[12]=desc[12]-desc[11]
  # remove the Nitro in nN
  desc[9]=desc[9]-desc[8]
  Polcoef = array([10.152, 8.765, 5.702, 3.833, 16.557, 24.123, 38.506, 10.488, 6.335, 4.307, 15.726, 22.366, 11.173])
  Polcoef0 = -1.529
  # polarizability
  Pol = sum(desc*Polcoef)+ Polcoef0 + 3.391*nH


  # MR molecular to molar: x 6.02214179*10^23
  MR=4/3*math.pi*Pol


  # lorentz-lorenz equation molar refractivity p 586
  # n = refractive index
  # MR (unit m3/mol) = (n^2-1)/(n^2+2)* MW / density = (de -1) / ( de+2 ) * MolVol = 4/3 * Pi * No * Polarizability
  # molarVol =  MW / density
  #print "toujours"
  #print "enfin"
  #print "nH:" + str(nH) 
  #print "MV:" + str(0.597+0.6823*Sumx)
  #print "V:" + str(Sumx)
  #print "MolMR:" + str(MolMR)
  #E2 = MolMR/10*Sumx/(0.597+0.6823*Sumx) - 2.83192*Sumx + 0.52553
  # for the moment I need to find how to have MV(Vx???)
  #MV =Sumx
  #E1 = MolMR/10*Sumx/MV - 2.83192*Sumx + 0.52553
  #print "E MR:" + str(E2)
  #print "E MolMR:" + str(E1)
  #print str(smi)
  
  # return a dict store of the elements for json
  result = {}
  result['Pol']=Pol
  result['MRx']=MR
  result['Vx']=Vx
  result['MolMR']=MolMR
  #result['MR']=d['MR']
  #print json.dumps(result)
  return result


def add_cat_ODT(mol):
  NB_fonction =5
  fonction_name = []
  fonction_smarts = [] 
  fonction_name.append('Mercaptan')
  fonction_smarts.append('[SX2H1]')
  fonction_name.append('CARBOXYLIC ACID')
  fonction_smarts.append('[CX3](=O)[OX2H1]')
  fonction_name.append('ALDEHYDE')
  fonction_smarts.append('[$([CX3H][#6]),$([CX3H2])]=[OX1]')
  fonction_name.append('ESTER')
  fonction_smarts.append('[CX3;$([R0][#6]),$([H1R0])](=[OX1])[OX2][#6;!$(C=[O,N,S])]')
  fonction_name.append('UNSATURATED')
  fonction_smarts.append('[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]=[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]')
  # vector of zeros
  desc = zeros(5)
  fonction_count2 = []
  p = Chem.MolFromSmarts('[!r]')
  i= 0
  while i < NB_fonction:
    func = Chem.MolFromSmarts(fonction_smarts[i])
    maps = mol.GetSubstructMatches(func)
    k=len(maps)
    if k>0:
      desc[i]=1
    i +=1
  return desc




###### ABRAHAMS & ODT

def add_Afrags(mol):
##  M = cdk.readstring("smi", str(smi))
  NB_fonction =51
  fonction_name = []
  fonction_smarts = []
  fonction_name.append('Alipahtic -OH 1')
  fonction_smarts.append('[C][OX2H]')
  fonction_name.append('phenol -OH 2')
  fonction_smarts.append('[c][OX2H]')
  fonction_name.append('Ali -NH2 3')
  fonction_smarts.append('[C][NX3;H2]')
  fonction_name.append('Aro -NH2 4')
  fonction_smarts.append('[c][NX3;H2;!$(NC=O)]')
  fonction_name.append('NH 5')
  fonction_smarts.append('[C][NX3;H1;!R][C]')
  fonction_name.append('NH 6')
  fonction_smarts.append('[C][NX3;H1;R][C]')
  fonction_name.append('Aniline NH 7')
  fonction_smarts.append('[c][NX3;H1;!$(NC=O)][C]')
  fonction_name.append('NH pyrrole 8')
  fonction_smarts.append('[c][nX3;H1][c]')
  fonction_name.append('Acide 9')
  fonction_smarts.append('[CX3](=O)[OX1H0-,OX2H1]')
  fonction_name.append('Amide I 10')
  fonction_smarts.append('[CX3](=[OX1])[NX3H2]')
  fonction_name.append('Amine II Ali 11')
  fonction_smarts.append('[CX3](=[OX1])[NX3;H1][C]')
  fonction_name.append('Amine II Aro 12')
  fonction_smarts.append('[CX3](=[OX1])[NX3;H1][c]')
  fonction_name.append('Thiamide 13')
  fonction_smarts.append('[$([SX4](=[OX1])(=[OX1])([!O])[NH,NH2,NH3+]),$([SX4+2]([OX1-])([OX1-])([!O])[NH,NH2,NH3+])]')
  fonction_name.append('urea 1 14')
  fonction_smarts.append('[NX3;H1]C(=[OX1])[NX3;H1]')
  fonction_name.append('urea 2 15')
  fonction_smarts.append('[NX3;H0]C(=[OX1])[NX3;H1]')
  fonction_name.append('carbamate 16')
  fonction_smarts.append('[NX3;H1]C(=[OX1])O')
  fonction_name.append('guanidine 17')
  fonction_smarts.append('[NX3;H1]C(=N)[NX3;H0]') 
  fonction_name.append('Alkyne 18')
  fonction_smarts.append('[C]#[CH]')
  fonction_name.append('Phosphoric acid 19')
  fonction_smarts.append('P[OH,O-]')
  fonction_name.append('RRCHX 20')
  fonction_smarts.append('[CH][F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)]')
  fonction_name.append('RCHX2 21')
  fonction_smarts.append('[CH]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)]')
  fonction_name.append('diacid 22')
  fonction_smarts.append('[CX4]([CX3](=O)[OX1H0-,OX2H1])[CX4][CX3](=O)[OX1H0-,OX2H1]')
  fonction_name.append('acid special 23')
  fonction_smarts.append('[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[CX3](=O)[OX1H0-,OX2H1]')
  fonction_name.append('alcohol special 24')
  fonction_smarts.append('[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[OH]')
  fonction_name.append('special special 25')
  fonction_smarts.append('[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[CX4][OH]')
  fonction_name.append('pyrazole type N 26')
  fonction_smarts.append('[nX3;H1]:n')
  fonction_name.append('imidazole type N 27')
  fonction_smarts.append('[nX3;H1]:c:n')
  fonction_name.append('H-bond 1 59')
  fonction_smarts.append('[OX2;H1]CC[O,N]')
  fonction_name.append('H-bond 2 60')
  fonction_smarts.append('[OX2;H1]C[C,N]=[O,S]')
  fonction_name.append('H-bond 3 61')
  fonction_smarts.append('[OX2;H1]c1ccccc1[O,NX3]')
  fonction_name.append('H-bond 4 62')
  fonction_smarts.append('[OX2;H1]c1ccccc1C=[O,S]')
  fonction_name.append('H-bond 5 63')
  fonction_smarts.append('[OX2;H1]c1ccccc1[$([NX3](=O)=O),$([NX3+](=O)[O-])]')
  fonction_name.append('H-bond 6 64')
  fonction_smarts.append('[NH,NH2,NH3+]CC[O,N]')
  fonction_name.append('H-bond 7 65')
  fonction_smarts.append('[NH,NH2,NH3+]c1ccccc1[O,N]')
  fonction_name.append('H-bond 8 66')
  fonction_smarts.append('[NH,NH2,NH3+]c1ccccc1[C,N]=[O,S]')
  fonction_name.append('H-bond 9 67')
  fonction_smarts.append('[OX2H]c1ccccc1[Cl,Br,I]')
  fonction_name.append('H-bond 10 37  37')
  fonction_smarts.append('[OX1]=[C,c]~[C,c]C[OH]')
  fonction_name.append('8-OH quinoline 38')
  fonction_smarts.append('[OH]c1cccc2cccnc12')
  fonction_name.append('3-X phenol 39')
  fonction_smarts.append('[OH]c1cc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])ccc1')
  fonction_name.append('4-X phenol 40')
  fonction_smarts.append('[OH]c1ccc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cc1')
  fonction_name.append('3-X aniline 41 ')
  fonction_smarts.append('[NH,NH2,NH3+]c1cc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])ccc1')
  fonction_name.append('4-X aniline 42')
  fonction_smarts.append('[NH,NH2,NH3+]c1ccc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cc1')
  fonction_name.append('3 X benzoic acid 43')
  fonction_smarts.append('[CX3](=O)([OX1H0-,OX2H1])c1cc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])ccc1')
  fonction_name.append('4 X benzoic acid 44')
  fonction_smarts.append('[CX3](=O)([OX1H0-,OX2H1])c1ccc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cc1')
  fonction_name.append('2,6 dialkyl phenol 45')
  fonction_smarts.append('[OH]c1c([CX4])cccc1[CX4]')
  fonction_name.append('2,6 dialkyl aniline 46')
  fonction_smarts.append('[NH,NH2,NH3+]c1c([CX4])cccc1[CX4]')
  fonction_name.append('2 CX phenol 47')
  fonction_smarts.append('[OH]c1c(C[F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cccc1')
  fonction_name.append('3 COOH phenol 48')
  fonction_smarts.append('[OH]c1cc([CX3](=O)[OX1H0-,OX2H1])ccc1')
  fonction_name.append('4 COOH phenol 49')
  fonction_smarts.append('[OH]c1ccc([CX3](=O)[OX1H0-,OX2H1])cc1')
  fonction_name.append('3 CO phenol 50')
  fonction_smarts.append('[OH]c1cc([$([CH](=O)),$(C(=O)C)])ccc1')
  fonction_name.append('4 CO phenol 51')
  fonction_smarts.append('[OH]c1ccc([$([CH](=O)),$(C(=O)C)])cc1')
  desc = zeros(51)
  i= 0
  while i < NB_fonction:
     func = Chem.MolFromSmarts(fonction_smarts[i])
     maps = mol.GetSubstructMatches(func)
     k =len(maps)
     desc[i]=k
     i +=1
  return desc


def add_Bfrags(mol):
  #M = cdk.readstring("smi", str(smi))
  NB_fonction =81
  fonction_name = []
  fonction_smarts = []
  fonction_name.append('CH3 sp3 1')
  fonction_smarts.append('[CX4H3]')
  fonction_name.append('CH2 sp3 2')
  fonction_smarts.append('[CX4H2]')
  fonction_name.append('CH sp3 3 ')
  fonction_smarts.append('[CX4H1]')
  fonction_name.append('-C- 4')
  fonction_smarts.append('[CX4H0]')
  fonction_name.append('CH2= 5')
  fonction_smarts.append('*=[CX3H2]')
  fonction_name.append('=CH- 6')
  fonction_smarts.append('[$(*=[CX3H1]),$([cX3H1](a)a)]')
  fonction_name.append('=C- 7')
  fonction_smarts.append('[$(*=[CX3H0]),$([cX3H0](a)(a)A)]')
  fonction_name.append('C fused aromatic 8')
  fonction_smarts.append('c(a)(a)a')
  fonction_name.append('C sp tb 9')
  fonction_smarts.append('*#C')                                                                                                                                                                                                               
  fonction_name.append('Amine 1 sp3 alip 10')
  fonction_smarts.append('[C][NX3;H2]')
  fonction_name.append('Amine 1 sp3 aro 11')
  fonction_smarts.append('[c][NX3;H2]')
  fonction_name.append('Amine 2 sp3 alip 12')
  fonction_smarts.append('[C][NX3;H1][C]')
  fonction_name.append('Amine 2 sp3 aro 13')
  fonction_smarts.append('[c][NX3;H1]')
  fonction_name.append('HN pyrrole 14')
  fonction_smarts.append('[c][nX3;H1][c]')
  fonction_name.append('Amine 3 sp3 alip 15')
  fonction_smarts.append('[C][NX3;H0](C)[C]')
  fonction_name.append('Amine 3 sp3 aro 16')
  fonction_smarts.append('[c][NX3;H0](C)[C]')
  fonction_name.append('RN pyrrole 17')
  fonction_smarts.append('[c][nX3;H0][c]')
  fonction_name.append('N sp2 alip 18')
  fonction_smarts.append('*=[Nv3;!R]')
  fonction_name.append('N sp2 cyclic 19')
  fonction_smarts.append('*=[Nv3;R]')
  fonction_name.append('Pyridine 20')
  fonction_smarts.append('[nX2H0,nX3H1+](a)a')
  fonction_name.append('Nitrile ali 21')
  fonction_smarts.append('N#C[A;!#1]')
  fonction_name.append('Nitrile aro 22')
  fonction_smarts.append('N#C[a;!#1]')
  fonction_name.append('NITRO ali 23')
  fonction_smarts.append('[$([A;!#1][NX3](=O)=O),$([A;!#1][NX3+](=O)[O-])]')
  fonction_name.append('NITRO aro 24')
  fonction_smarts.append('[$([a;!#1][NX3](=O)=O),$([a;!#1][NX3+](=O)[O-])]')
  fonction_name.append('NITRATE 25')
  fonction_smarts.append('[$([NX3](=[OX1])(=[OX1])O),$([NX3+]([OX1-])(=[OX1])O)]')
  fonction_name.append('any -OH 26')
  fonction_smarts.append('[OH]')
  fonction_name.append('O sp3 noncyclic 27')
  fonction_smarts.append('[OX2;H0;!R]')
  fonction_name.append('O sp3 cyclic 28')
  fonction_smarts.append('[OX2;H0;R]')
  fonction_name.append('-O- aromatic 29')
  fonction_smarts.append('[oX2](a)a')
  fonction_name.append('=O sp2 30')
  fonction_smarts.append('*=O')
  fonction_name.append('-S- sp3 31')
  fonction_smarts.append('[SX2](*)*')
  fonction_name.append('-S- aromatic 32')
  fonction_smarts.append('[sX2](a)a')
  fonction_name.append('S= sp2 33')
  fonction_smarts.append('*=[SX1]')
  fonction_name.append('--S= Sp2 34')
  fonction_smarts.append('*=[SX3]')
  fonction_name.append('sulfonate 35')
  fonction_smarts.append('[$([#16X4](=[OX1])(=[OX1])([!#8])[OX2H0]),$([#16X4+2]([OX1-])([OX1-])([!#8])[OX2H0])]')
  fonction_name.append('soufre 36')
  fonction_smarts.append('[S,s]')                                                                                                                                                                                     
  fonction_name.append('phosphore 37')
  fonction_smarts.append('[P,p]')
  fonction_name.append('F alipathic 38')
  fonction_smarts.append('FA')
  fonction_name.append('F aromatic 39')
  fonction_smarts.append('Fa')
  fonction_name.append('Cl 40')
  fonction_smarts.append('Cl')
  fonction_name.append('Br 41 ')
  fonction_smarts.append('Br')
  fonction_name.append('I 42')
  fonction_smarts.append('I')
  fonction_name.append('Non cycli ester 43')
  fonction_smarts.append('[CX3;!R](=[OX1])[OX2H0]')
  fonction_name.append('Lactone 44')
  fonction_smarts.append('[CX3;R](=[OX1])[OX2H0;R]')
  fonction_name.append('Phosphate 45')
  fonction_smarts.append('P(=[OX1])(O)(O)O')
  fonction_name.append('carbonate 46')
  fonction_smarts.append('[CX3](=[OX1])([OX2H0])[OX2H0]')
  fonction_name.append('carboxylique acid 47')
  fonction_smarts.append('[CX3](=O)[OX1H0-,OX2H1]')
  fonction_name.append('Amide Aromatic 48')
  fonction_smarts.append('nC=[OX1]')
  fonction_name.append('Amide non cyclic aliphatic 49')
  fonction_smarts.append('[N;!R]C=[OX1]')
  fonction_name.append('Lactam 50')
  fonction_smarts.append('[N;R][C;R]=[OX1]')
  fonction_name.append('sulfonamide 51')
  fonction_smarts.append('[$([SX4](=[OX1])(=[OX1])([!O])[NX3]),$([SX4+2]([OX1-])([OX1-])([!O])[NX3])]')
  fonction_name.append('urea 52')
  fonction_smarts.append('NC(=[OX1])N')
  fonction_name.append('carbamate 53')
  fonction_smarts.append('[NX3,NX4+][CX3](=[OX1])[OX2,OX1-]')
  fonction_name.append('imide 54')
  fonction_smarts.append('[CX3](=[OX1])[NX3][CX3](=[OX1])')
  fonction_name.append('quinone 55')
  fonction_smarts.append('C1(=[OX1])C=CC(=[OX1])C=C1')
  fonction_name.append('-CX2- 56')
  fonction_smarts.append('[$([CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])]')
  fonction_name.append('--XCX-CX-- 57')
  fonction_smarts.append('[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[CX4][F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)]')
  fonction_name.append('Steroide 58')
  fonction_smarts.append('*1~*2~*(~*3~*(~*~*~*~*3)~*1)~*~*~*1~*2~*~*~*1')                                                                                                                                                                             
  fonction_name.append('H-bond 1 59')
  fonction_smarts.append('[OX2;H1]CC[O,N]')
  fonction_name.append('H-bond 2 60')
  fonction_smarts.append('[OX2;H1]C[C,N]=[O,S]')
  fonction_name.append('H-bond 3 61')
  fonction_smarts.append('[OX2;H1]c1ccccc1[O,NX3]')
  fonction_name.append('H-bond 4 62')
  fonction_smarts.append('[OX2;H1]c1ccccc1C=[O,S]')
  fonction_name.append('H-bond 5 63')
  fonction_smarts.append('[OX2;H1]c1ccccc1[$([NX3](=O)=O),$([NX3+](=O)[O-])]')
  fonction_name.append('H-bond 6 64')
  fonction_smarts.append('[NH,NH2,NH3+]CC[O,N]')
  fonction_name.append('H-bond 7 65')
  fonction_smarts.append('[NH,NH2,NH3+]c1ccccc1[O,N]')
  fonction_name.append('H-bond 8 66')
  fonction_smarts.append('[NH,NH2,NH3+]c1ccccc1[C,N]=[O,S]')
  fonction_name.append('H-bond 9 67')
  fonction_smarts.append('[OX2H]c1ccccc1[Cl,Br,I]')
  fonction_name.append('1,2 diol 68')
  fonction_smarts.append('[CX4]([OH])[CX4][OH]')
  fonction_name.append('pyrazine interaction 69')
  fonction_smarts.append('n:n')
  fonction_name.append('isoxazole interaction 70')
  fonction_smarts.append('o:n')
  fonction_name.append('pyrimidine interaction 71')
  fonction_smarts.append('n:c:n')
  fonction_name.append('oxazole interaction 72')
  fonction_smarts.append('o:c:n')
  fonction_name.append('pyrazine interaction 73')
  fonction_smarts.append('n:c:c:n')
  fonction_name.append('ortho intercation 74')
  fonction_smarts.append('[F,Cl,Br,I,N,O,S]-c:c-[F,Cl,Br,I,N,O,S]')
  fonction_name.append('meta intercation 75')
  fonction_smarts.append('[F,Cl,Br,I,N,O,S]-c:c:c-[F,Cl,Br,I,N,O,S]')
  fonction_name.append('para interaction 76')
  fonction_smarts.append('[F,Cl,Br,I,N,O,S]-c:c:c:c-[F,Cl,Br,I,N,O,S]')
  fonction_name.append('phosphamide 77')
  fonction_smarts.append('P(=[OX1])N')
  fonction_name.append('2 amino pyridine 78')
  fonction_smarts.append('Nc:n')
  fonction_name.append('benzyl alcohol 79')
  fonction_smarts.append('[$(cC[OH]);!$(c[CX3](=O)[OX1H0-,OX2H1])]')
  fonction_name.append('N-oxide 80')
  fonction_smarts.append('[$([#7+][OX1-]),$([#7v5]=[OX1]);!$([#7](~[O])~[O]);!$([#7]=[#7])]')
  fonction_name.append('1,2-dimethoxy 81')
  fonction_smarts.append('[OX2]-c:c-[OX2]')

  desc = zeros(81)
  i= 0
  while i < NB_fonction:
    func = Chem.MolFromSmarts(fonction_smarts[i])
    maps = mol.GetSubstructMatches(func)
    k =len(maps)
    desc[i]=k
    i +=1
  return desc

def compute_ODT(mol):
	# Linear
    Ecoef = array([-0.104,0,0.089,0.187,-0.045,0.068,0.18,0.3,0.04,0.085,0.163,0.138,0.192,-0.03,0.22,0.346,0.083,0.117,0.121,0.046,0,0,0.2,0.21,0,0.061,0.014,0.013,-0.125,-0.041,0.33,0.116,0.364,0.413,0,0.465,0.295,-0.18,-0.23,0.023,0.196,0.533,-0.113,0,-0.1,0,-0.192,0.221,0,0.061,-0.111,-0.11,0,0,0,-0.017,0.012,0.285,0.029,0,-0.069,0,0,0,0,0,-0.1,-0.043,0.092,-0.113,0,0.052,0,0,0,0,-0.08,0.185,0,0,0])
    Ecoef0 = 0.248
    Scoef = array([-0.075,  0,  0.036,  0.071,  -0.085, 0.05, 0.101,  0.121,  0.034,  0.175,  0.383,  0.265,  0.311,  0.221,  0.323,  0.295,  0.265,  0.125,  0.254,  0.223,  0.694,  0.39, 0,  -0.231, -0.476, 0.247,  0.185,  0.185,  0,  0.37, 0.189,  0,  0.618,  1.065,  -0.505, 0.643,  0.703,  -0.042, 0,  0.082,  0.161,  0.198,  -0.225, 0.36, -0.24,  -0.19,  -0.412, -0.076, 0.175,  -0.1, -0.569, -0.553, -0.588, -0.51,  -0.411, -0.05,  0,  1.029,  -0.067, -0.095, -0.237, -0.344, -0.276, -0.102, 0,  -0.14,  -0.12,  0.052,  0.024,  0.047,  -0.04,  0.087,  -0.051, -0.043, -0.038, 0,  -0.452, 0.098,  0,  0.434,  0.38])
    Scoef0 = 0.277
    BHcoef = array([0.007,0,0.011,0.037,0.019,0.011,0,0.019,0.028,0.481,0.275,0.541,0.415,0.316,0.653,0.321,0.392,0.2,0.596,0.321,0.242,0.103,-0.476,-0.525,-0.204,0.307,0.211,0.331,0.047,0.334,0.168,0.043,0.071,0.448,-0.188,0,1.183,-0.036,0,0,-0.011,0,-0.206,-0.214,-0.394,-0.267,-0.308,-0.095,-0.287,-0.231,-0.446,-0.076,-0.252,-0.148,-0.051,-0.014,0.013,0.267,0,-0.068,-0.079,-0.387,-0.126,0,-0.059,-0.045,-0.13,0,-0.132,-0.157,-0.098,-0.17,-0.089,0.031,-0.035,-0.023,-0.668,-0.042,0.131,-0.408,-0.216])
    BHcoef0 = 0.071
    #BOcoef = array([0,0,0.02,0.047,0.024,0.012,0,0.018,0.032,0.486,0.326,0.543,0.426,0.267,0.655,0.338,0.338,0.202,0.589,0.3,0.245,0.093,-0.595,-0.533,-0.202,0.311,0.226,0.33,0.06,0.339,0.175,0.083,0.069,0.319,-0.19,0,1.189,-0.033,0,0,0,0,-0.223,-0.169,-0.408,-0.298,-0.312,-0.038,-0.292,-0.242,-0.443,-0.054,-0.251,-0.149,-0.05,-0.016,0.01,0.218,0,-0.09,-0.122,-0.403,-0.12,0,-0.027,-0.069,-0.13,-0.018,-0.094,-0.141,-0.113,-0.184,-0.073,0.025,-0.033,-0.025,-0.668,-0.057,0.129,-0.405,-0.218])
    #BOcoef0 = 0.064
    Lcoef = array([0.321,0.499,0.449,0.443,0.244,0.469,0.624,0.744,0.332,0.781,0.949,0.568,0.912,1.25,0.4,0.869,0.794,-0.235,-0.24,0.574,0.757,0.732,0.278,0.347,0,0.672,0.36,0.359,0.057,0.495,1.258,0.848,0.954,2.196,0,0.554,2.051,-0.143,-0.147,0.669,1.097,1.59,-0.39,0.406,-0.483,0,-0.369,0,0.603,0.583,0,0,0,0,0,-0.111,0.054,0.488,-0.072,-0.337,0,-0.303,-0.364,0.062,0,0.169,-0.4,0.1,-0.179,0,0.042,0.209,-0.058,-0.081,-0.026,0,0,0.149,-0.145,0,0])
    Lcoef0 = 0.13
    Acoef = array([0.345,0.543,0.177,0.247,0.087,0.321,0.194,0.371,0.243,0.275,0.281,-0.091,0.356,-0.165,-0.119,-0.105,0.17,0.082,0.493,0.019,0.05,-0.362,0.118,0.1,0.051,0.194,0.042,-0.089,-0.161,-0.251,-0.418,-0.45,-0.155,0,-0.093,-0.11,-0.601,-0.475,0.119,0.176,0.08,0.084,0.085,0.055,-0.162,-0.181,0.195,-0.203,0.096,0.185,0.203])
    Acoef0 = 0.003
    # 81 fragments for BH, BO, L, S, E
    odesc = add_Bfrags(mol)
    # correction to remove FRAGMENT 36 IF (31 OR 32 OR 33 OR 34 OR 35) NOT NULL
    if odesc[30]+odesc[31]+odesc[32]+odesc[33]+odesc[34] >0:
      odesc[35]=0

    # 51 fragments for A
    o2desc = add_Afrags(mol)
    E = sum(odesc*Ecoef)+ Ecoef0
    S = sum(odesc*Scoef)+Scoef0
    #BO = sum(odesc*BOcoef)+ BOcoef0
    B = sum(odesc*BHcoef)+ BHcoef0
    L = sum(odesc*Lcoef) + Lcoef0
    A = sum(o2desc*Acoef) + Acoef0
    cat = add_cat_ODT(mol)
    M = cat[0]
    AC = cat[1]
    AL = cat[2]
    UE = cat[3]*cat[4]
    minuslogODT = -1.826+ 0.882*E+0.408*S+0.999*A+2.196*B+0.578*L+4.065*M+1.805*AL+1.424*AC+1.290*UE

    # return a dict store of the elements for json
    result = {}
    result['E']=E
    result['S']=S
    result['B']=B
    result['L']=L
    result['A']=A
    result['minuslogODT']=minuslogODT
    result['M']=M
    result['AC']=AC
    result['AL']=AL
    result['UE']=UE
    #print json.dumps(result)
    return result

def writemol(mol):
  string = "\n"
  if mol != None:
      dict1 = add_V_E(mol)
      dict2 = compute_ODT(mol)
  z = dict1.copy()
  z.update(dict2)
  return z

def main(argv):
    smile = ''
    inputfile = ''
    descriptor = ''
    res = {}

    try:
      opts, args = getopt.getopt(argv,"hs:d:i:",["smile=","descriptors=","inputfile="])
    except getopt.GetoptError:
      print ('SomeDescritpors.py options can be : -s <smilevalue> -d <descriptor ie rdkit, angle, tgg> -i <inputfile>')
      sys.exit(2)
    for opt, arg in opts:
      if opt == '-h':
         print ('SomeDescritpors.py -s <smilevalue> -d <descriptor> -i <inputfile>')
         sys.exit()
      elif opt in ("-s", "--smile"):
         smile = arg
      elif opt in ("-d", "--descriptors"):
         descriptor = arg
      elif opt in ("-i", "--inputfile"):
         inputfile = arg

    Descs=descriptor.split(',')

    if len(smile)>1:
      mol = Chem.MolFromSmiles(str(smile))

    # switch for sdf file or smile
    if len(inputfile)>1:
      suppl = Chem.SDMolSupplier(inputfile)
      mols = [x for x in suppl]
      mol = mols[0]
      smile = Chem.MolToSmiles(mol,True)

    if 'RDKIT' in Descs or 'rdkit' in Descs:
      res.update(rdkitcalcdesc(mol))

    if 'ANGLE' in Descs or 'angle' in Descs:
      res.update(angle_HH(mol))

    if 'TGG' in Descs or 'tgg' in Descs:
      res3 = writemol(mol)
      res.update(res3)

    print (json.dumps(res, ensure_ascii=False))

if __name__=="__main__":
   main(sys.argv[1:])

