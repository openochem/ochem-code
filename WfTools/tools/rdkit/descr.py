import numpy as np
import re

def nans (dim):
    x = np.empty(dim, np.float32)
    x.fill(np.nan)
    return x

from rdkit import Chem
import rdkit.Chem.rdMolDescriptors as rdMD
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MACCSkeys
from rdkit.Chem import AllChem
from rdkit import Avalon
from rdkit.Avalon.pyAvalonTools import *
###

from  rdkit.Chem.rdMolDescriptors import CalcAsphericity, CalcChi0n, CalcChi0v, CalcChi1n, CalcChi1v, CalcChi2n, CalcChi2v, CalcChi3n, CalcChi3v, CalcChi4n, CalcChi4v, CalcChiNn, CalcChiNv, CalcEccentricity, CalcExactMolWt, CalcHallKierAlpha, CalcFractionCSP3, CalcInertialShapeFactor, CalcKappa1, CalcKappa2, CalcKappa3, CalcLabuteASA, CalcNPR1, CalcNPR2, CalcNumAliphaticCarbocycles, CalcNumAliphaticHeterocycles, CalcNumAliphaticRings, CalcNumAmideBonds, CalcNumAromaticCarbocycles, CalcNumAromaticHeterocycles, CalcNumAromaticRings, CalcNumAtomStereoCenters, CalcNumBridgeheadAtoms, CalcNumHBA, CalcNumHBD, CalcNumHeteroatoms, CalcNumHeterocycles, CalcNumLipinskiHBA, CalcNumLipinskiHBD, CalcNumRings, CalcNumRotatableBonds, CalcNumSaturatedCarbocycles, CalcNumSaturatedHeterocycles, CalcNumSaturatedRings, CalcNumSpiroAtoms, CalcNumUnspecifiedAtomStereoCenters, CalcPBF, CalcPMI1, CalcPMI2, CalcPMI3, CalcRadiusOfGyration, CalcSpherocityIndex, CalcTPSA
from rdkit.Chem.Crippen import MolMR, MolLogP # two valoes for CalcCrippenDescriptors

scalar_fns = [CalcAsphericity, CalcChi0n, CalcChi0v, CalcChi1n, CalcChi1v, CalcChi2n, CalcChi2v, CalcChi3n, CalcChi3v, CalcChi4n, CalcChi4v, CalcEccentricity, CalcExactMolWt, CalcHallKierAlpha, CalcFractionCSP3, CalcInertialShapeFactor, CalcKappa1, CalcKappa2, CalcKappa3, CalcLabuteASA, CalcNPR1, CalcNPR2, CalcNumAliphaticCarbocycles, CalcNumAliphaticHeterocycles, CalcNumAliphaticRings, CalcNumAmideBonds, CalcNumAromaticCarbocycles, CalcNumAromaticHeterocycles, CalcNumAromaticRings, CalcNumAtomStereoCenters, CalcNumBridgeheadAtoms, CalcNumHBA, CalcNumHBD, CalcNumHeteroatoms, CalcNumHeterocycles, CalcNumLipinskiHBA, CalcNumLipinskiHBD, CalcNumRings, CalcNumRotatableBonds, CalcNumSaturatedCarbocycles, CalcNumSaturatedHeterocycles, CalcNumSaturatedRings, CalcNumSpiroAtoms, CalcNumUnspecifiedAtomStereoCenters, CalcPBF, CalcPMI1, CalcPMI2, CalcPMI3, CalcRadiusOfGyration, CalcSpherocityIndex, CalcTPSA, MolLogP, MolMR]

scalar_names = ['CalcAsphericity', 'Chi0n', 'Chi0v', 'Chi1n', 'Chi1v', 'Chi2n', 'Chi2v', 'Chi3n', 'Chi3v', 'Chi4n', 'Chi4v', 'Eccentricity', 'ExactMolWt', 'HallKierAlpha', 'FractionCSP3', 'InertialShapeFactor', 'Kappa1', 'Kappa2', 'Kappa3', 'LabuteASA', 'NPR1', 'NPR2', 'NumAliphaticCarbocycles', 'NumAliphaticHeterocycles', 'NumAliphaticRings', 'NumAmideBonds', 'NumAromaticCarbocycles', 'NumAromaticHeterocycles', 'NumAromaticRings', 'NumAtomStereoCenters', 'NumBridgeheadAtoms', 'NumHBA', 'NumHBD', 'NumHeteroatoms', 'NumHeterocycles', 'NumLipinskiHBA', 'NumLipinskiHBD', 'NumRings', 'NumRotatableBonds', 'NumSaturatedCarbocycles', 'NumSaturatedHeterocycles', 'NumSaturatedRings', 'NumSpiroAtoms', 'NumUnspecifiedAtomStereoCenters', 'PBF', 'PMI1', 'PMI2', 'PMI3', 'RadiusOfGyration', 'SpherocityIndex', 'TPSA', 'MolLogP', 'MolMR']

# scalar_names = list(map(lambda fn: fn.func_name[4:], scalar_fns))[:-2] + ['MolLogP', 'MolMR']

def scalars (m):
    res = []
    for fn in scalar_fns:
        try:
            res.append(fn(m))
        except:
            res.append(np.nan)
    return res

###

from rdkit.Chem.EState.EState import MaxEStateIndex, MinEStateIndex, MaxAbsEStateIndex, MinAbsEStateIndex
from rdkit.Chem.EState.EState_VSA import EState_VSA1, EState_VSA2, EState_VSA3, EState_VSA4, EState_VSA5, EState_VSA6, EState_VSA7, EState_VSA8, EState_VSA9, EState_VSA10, EState_VSA11
from rdkit.Chem.MolSurf import PEOE_VSA1, PEOE_VSA2, PEOE_VSA3, PEOE_VSA4, PEOE_VSA5, PEOE_VSA6, PEOE_VSA7, PEOE_VSA8, PEOE_VSA9, PEOE_VSA10, PEOE_VSA11, PEOE_VSA12, PEOE_VSA13,  PEOE_VSA14
from rdkit.Chem.MolSurf import SMR_VSA1, SMR_VSA2, SMR_VSA3, SMR_VSA4, SMR_VSA5, SMR_VSA6, SMR_VSA7, SMR_VSA8, SMR_VSA9, SMR_VSA10
from rdkit.Chem.MolSurf import SlogP_VSA1, SlogP_VSA2, SlogP_VSA3, SlogP_VSA4, SlogP_VSA5, SlogP_VSA6, SlogP_VSA7, SlogP_VSA8, SlogP_VSA9, SlogP_VSA10, SlogP_VSA11, SlogP_VSA12
from rdkit.Chem.Lipinski import FractionCSP3, HeavyAtomCount,  NHOHCount, NHOHSmarts, NOCount, NumHAcceptors, NumHDonors
from rdkit.Chem.QED import qed as QED
from rdkit.Chem.GraphDescriptors import BalabanJ, BertzCT, Ipc

secondary_scalar_fns = [FractionCSP3, HeavyAtomCount,  NHOHCount, NOCount, NumHAcceptors, NumHDonors,
                        MaxEStateIndex, MinEStateIndex, MaxAbsEStateIndex, MinAbsEStateIndex,
                        BalabanJ, BertzCT, Ipc,
                        QED,
                        EState_VSA1, EState_VSA2, EState_VSA3, EState_VSA4, EState_VSA5, EState_VSA6, EState_VSA7, EState_VSA8, EState_VSA9, EState_VSA10, EState_VSA11,
                        PEOE_VSA1, PEOE_VSA2, PEOE_VSA3, PEOE_VSA4, PEOE_VSA5, PEOE_VSA6, PEOE_VSA7, PEOE_VSA8, PEOE_VSA9, PEOE_VSA10, PEOE_VSA11, PEOE_VSA12, PEOE_VSA13,  PEOE_VSA14,
                        SMR_VSA1, SMR_VSA2, SMR_VSA3, SMR_VSA4, SMR_VSA5, SMR_VSA6, SMR_VSA7, SMR_VSA8, SMR_VSA9, SMR_VSA10,
                        SlogP_VSA1, SlogP_VSA2, SlogP_VSA3, SlogP_VSA4, SlogP_VSA5, SlogP_VSA6, SlogP_VSA7, SlogP_VSA8, SlogP_VSA9, SlogP_VSA10, SlogP_VSA11, SlogP_VSA12]

secondary_scalar_names = ['FractionCSP3', 'HeavyAtomCount',  'NHOHCount', 'NOCount', 'NumHAcceptors', 'NumHDonors',
                        'MaxEStateIndex', 'MinEStateIndex', 'MaxAbsEStateIndex', 'MinAbsEStateIndex',
                        'BalabanJ', 'BertzCT', 'Ipc',
                        'QED',
                        'EState_VSA1', 'EState_VSA2', 'EState_VSA3', 'EState_VSA4', 'EState_VSA5', 'EState_VSA6', 'EState_VSA7', 'EState_VSA8', 'EState_VSA9', 'EState_VSA10', 'EState_VSA11',
                        'PEOE_VSA1', 'PEOE_VSA2', 'PEOE_VSA3', 'PEOE_VSA4', 'PEOE_VSA5', 'PEOE_VSA6', 'PEOE_VSA7', 'PEOE_VSA8', 'PEOE_VSA9', 'PEOE_VSA10', 'PEOE_VSA11', 'PEOE_VSA12', 'PEOE_VSA13',  'PEOE_VSA14',
                        'SMR_VSA1', 'SMR_VSA2', 'SMR_VSA3', 'SMR_VSA4', 'SMR_VSA5', 'SMR_VSA6', 'SMR_VSA7', 'SMR_VSA8', 'SMR_VSA9', 'SMR_VSA10',
                        'SlogP_VSA1', 'SlogP_VSA2', 'SlogP_VSA3', 'SlogP_VSA4', 'SlogP_VSA5', 'SlogP_VSA6', 'SlogP_VSA7', 'SlogP_VSA8', 'SlogP_VSA9', 'SlogP_VSA10', 'SlogP_VSA11', 'SlogP_VSA12']
    
#secondary_scalar_names = map(lambda fn: fn.func_name if fn.func_name != "<lambda>" else fn.func_doc.split("(")[0].split('Calc')[-1], secondary_scalar_fns)
#secondary_scalar_names = list(map(lambda x: x.replace(' ', '_'), secondary_scalar_names))

def scalars_secondary (m):
    res = []
    for fn in secondary_scalar_fns:
        try:
            res.append(fn(m))
        except:
            res.append(np.nan)
    return res


from rdkit.Chem.Fragments import fns as Fragments

fragments_scalar_names = list(map(lambda x: "fragment_"+x[0], Fragments))

def scalars_fragments (m):
    res = []
    for name, fn in Fragments:
        try:
            res.append(fn(m))
        except:
            res.append(np.nan)
    return res

### This does not work in singularity image either :(
###
#
#from rdkit.Chem import FragmentCatalog
#fparams = FragmentCatalog.FragCatParams(1,6, 'FunctionalGroups.txt')
#fcat=FragmentCatalog.FragCatalog(fparams)
#fpgen = FragmentCatalog.FragFPGenerator()
#
#def scalars_fragments (m):
#    return [x for x in fpgen.GetFPForMol(m, fcat)]

#def scalars_fragments (m):
#    return []


### processing sparse vector
def sparse_process (fn, task, mols, name):
    keys = []
    results = []
    for m in mols:
        try:
            res = fn(m).GetNonzeroElements()
        except:
            print("Except!")
            res = {}
        results.append(res)
        keys.append(list(res.keys()))

    uniq_keys = set()
    for k in keys:
        uniq_keys |= set(k)
    uniq_keys = list(uniq_keys)

    def fingerprint (r):
        fp = dict(zip(uniq_keys,[0]*len(uniq_keys)))
        for k in r.keys():
            fp[k] = r[k]
        return list(fp.values())

    if name in ['atom_pairs', 'bpf', 'btf']:
        names = list(map(lambda id: name+"_"+re.sub("['\ ]", "", "{}".format(Pairs.ExplainPairScore(id))), uniq_keys))
    elif name == 'torsions':
        names = list(map(lambda id: name+"_"+re.sub("['\ ]", "", "{}".format(Torsions.ExplainPathScore(id))), uniq_keys))
        
    else:
        names = list(map(lambda x: name+"_"+str(x), uniq_keys))
    fps = []
    for r in results:
        fps.append(fingerprint(r))
    return fps, names

def morgan_counts (task, mols):
    return  sparse_process (lambda m: AllChem.GetMorganFingerprint(m, task['morgan_radius'], useFeatures=task['morgan_fcfp']), task, mols, 'morgan_counts')

def avalon_counts (task,mols):
    print("calculate avalon_counts")
    d = sparse_process (lambda m: GetAvalonCountFP(m, nBits=task['avalon_nbits']), task, mols, 'avalon_counts')
    print("finished")
    return d

import rdkit.Chem.AtomPairs.Pairs as Pairs
def atom_pairs (task, mols):
    return sparse_process (Pairs.GetAtomPairFingerprint, task, mols, 'atom_pairs')

import rdkit.Chem.AtomPairs.Sheridan as Sheridan
def btf (task, mols):
    return sparse_process (Sheridan.GetBTFingerprint, task, mols, 'btf')

def bpf (task, mols):
    return sparse_process (Sheridan.GetBPFingerprint, task, mols, 'bpf')

import rdkit.Chem.AtomPairs.Torsions as Torsions

def torsions (task, mols):
    return sparse_process (Torsions.GetTopologicalTorsionFingerprint, task, mols, 'torsions')
###

from sascore import calculateScore as _sascore
def sascore (m):
    return [_sascore(m)]

sascore_names = ['sascore']

###

def topological (m, nbits):
    def _topological (m):
        return FingerprintMols.FingerprintMol(m, minPath=1, maxPath=7, fpSize=nbits, bitsPerHash=2, useHs=False, tgtDensity=0.3, minSize=nbits, maxSize=nbits)
    
    return [x for x in _topological(m)]

###

def make_descr_funcs (task):
    descr_funcs = { 'autocorr2d' : rdMD.CalcAUTOCORR2D,
                    'autocorr3d' : rdMD.CalcAUTOCORR3D,
                    'getaway' : rdMD.CalcGETAWAY,
                    'morse' : rdMD.CalcMORSE,
                    'rdf' : rdMD.CalcRDF,
                    'whim' : lambda m: rdMD.CalcWHIM(m, thresh=task['whim_thresh']),
                    #'maccs' : MACCSkeys.GenMACCSKeys,
                    'maccs' : AllChem.GetMACCSKeysFingerprint,
                    'topological' : lambda m: topological(m, task['topological_nbits']),
                    'morgan' : lambda m: [x for x in AllChem.GetMorganFingerprintAsBitVect(m, task['morgan_radius'], nBits=task['morgan_nbits'], useFeatures=task['morgan_fcfp'])],
                    'avalon' : lambda m: [x for x in GetAvalonFP(m, nBits=task['avalon_nbits'])],
                    'rdkitdef' : lambda m: [x for x in AllChem.RDKFingerprint(m)],
                    'sascore' : sascore,
                    'scalars' : scalars,
                    'scalars_secondary' : scalars_secondary,
                    'scalars_fragments' : scalars_fragments,
                    
                    'avalon_counts' : avalon_counts,
                    'morgan_counts' : morgan_counts,
                    'atom_pairs' : atom_pairs,
                    'btf' : btf,
                    'bpf' : bpf,
                    'torsions' : torsions
                        }

    task['descr_funcs'] = descr_funcs

name_lists = {'sascore' : sascore_names,
              'scalars' : scalar_names,
              'scalars_secondary' : secondary_scalar_names,
              'scalars_fragments' : fragments_scalar_names}

sparse_funcs = ['morgan_counts', 'avalon_counts', 'atom_pairs', 'btf', 'bpf', 'torsions']

def run_descr_sparse (task, mols, descr_name):
    return task['descr_funcs'][descr_name](task, mols)
    

def run_descr (task, mols, descr_name):
    res = []
    dim = 0
    for m in mols:
        try:
            m = removeRadicals(m)
            row = task['descr_funcs'][descr_name](m)
            dim = len(row)
            res.append(row)
        except:
            res = nans(dim)
#        except Exception as e:
#            res.append("error: "+e.__str__())
    if dim > 0 and len(res) > 0:
        names = name_lists.get(descr_name)
        if names is None:
            names = []
            for i in range(dim):
                names.append(descr_name+str(i))

        for i in range(len(res)):
            if len(res[i]) != dim:
                res[i] = nans(dime)
        return res,names
    else:
        return None,None

def removeRadicals(mol):
    try:
        for a in mol.GetAtoms():
            a.SetNumRadicalElectrons(0)
    except:
        pass
    return mol

def read_mols (task):
    try:
        mols = Chem.SDMolSupplier(task['input_file'],sanitize=False)
        mm = []
        for m in mols:
            try:
                Chem.SanitizeMol(m,Chem.SanitizeFlags.SANITIZE_FINDRADICALS|Chem.SanitizeFlags.SANITIZE_SETAROMATICITY|Chem.SanitizeFlags.SANITIZE_SETCONJUGATION|Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION|Chem.SanitizeFlags.SANITIZE_SYMMRINGS,catchErrors=True)
            except:
                pass
            mm.append(m)
    except:
        print ("Bad input file: "+task['input_file'])
    return mm

def run (task):
    make_descr_funcs(task)
    mols = read_mols(task)
    results = []
    for descr_name in list(task['descr_funcs'].keys()):
        if task[descr_name]:
            if descr_name in sparse_funcs:
                res, names = run_descr_sparse(task, mols, descr_name)
            else:
                res,names = run_descr(task, mols, descr_name)
            if res is not None:
                results.append([res, names])
#    header = flatten(map(lambda x: x[1], results))
#    body = flatten(map(lambda x: x[0], results))
    with open(task['output_file'], 'w') as fh:
        for names in list(map(lambda x: x[1], results)):
            for name in names:
                fh.write(name+" ")
        fh.write("\n")
        bodies = list(map(lambda x: x[0], results))
        for i in range(len(mols)):
            for res in bodies:
                for val in res[i]:
                    fh.write(str(float(val))+" ")
            fh.write("\n")
    return results
