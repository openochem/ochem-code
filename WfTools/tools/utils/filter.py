import os
import csv
import sys
import re
import numpy as np
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MolFromSmiles,MolToSmiles
import multiprocessing
import time
import pathlib
import random

global isomeric
global methodtype
global sanitize
global dir
global mg

default_dict = {
    '#': 1,
    '(': 2,
    ')': 3,
    '+': 4,
    '-': 5,
    '/': 6,
    '1': 7,
    '2': 8,
    '3': 9,
    '4': 10,
    '5': 11,
    '6': 12,
    '7': 13,
    '8': 14,
    '=': 15,
    'C': 16,
    'F': 17,
    'H': 18,
    'I': 19,
    'N': 20,
    'O': 21,
    'P': 22,
    'S': 23,
    '[': 24,
    '\\': 25,
    ']': 26,
    '_': 27,
    'c': 28,
    'Cl': 29,
    'Br': 30,
    'n': 31,
    'o': 32,
    's': 33
}

allowedSymbols = "#%()+-./0123456789=ABCFGHIKLMNOPRSTUVZ[\]abcdefghilnorstu"

def RS(smi, center):
    # We have this automatically # of center = # of setAtomMapNum <=> FindMolChiralCenters!
    # we consider only real chiral centers (having :digit) previously assigned
    smi = re.sub('(@{1,2})(H?[+-]?:[0-9]+)',r'^\2', smi) # first step put all to "^" ie "R" 
    # if we have a S center we change it ^ to _
    for p in center:
        if p[1] == 'S':            
            t0 = p[0]+1
            t = '(\^)(H?[+-]?:'+str(t0)+')'
            smi = re.sub(t,r'_\2',smi)
    return re.sub(':[0-9]+','',smi)


def SMILESToChiSMI(smi):
    if "@" not in smi:
        return smi
    m = Chem.MolFromSmiles(smi)
    
    for i in range(2):
        if i == 1:
            p = multiprocessing.Process(target=Chem.rdCIPLabeler.AssignCIPLabels, args=(m))
            p.start()
            p.join(5)

            if p.is_alive():
                print ("still running... let's kill it...")
                p.terminate()
                p.join()
            
        center = Chem.FindMolChiralCenters(m)
    
        idxcenter=[chiidx for (chiidx, chivalue) in center]

        for idx,at in enumerate(m.GetAtoms()):
            if idx in idxcenter:
                at.SetAtomMapNum(idx+1)
        
        if m is  None:
            raise ValueError('Failed to produce molecule from ' + sm)

        smi1 = Chem.MolToSmiles(m, doRandom = False, canonical = False)
    
        chi = RS(smi1, center)
    
        if "@" not in chi:
            return chi
    
    raise Exception("conversion failed: " + chi)

def smiles_to_seq (smiles, char_dict,length):
    smiles_len = len(smiles)
    if smiles_len >= length*1.1:
        raise ValueError('too long SMILES')
    seq = [0]
    keys = char_dict.keys()
    i = 0
    while i < smiles_len:
        # Skip all spaces
        if smiles[i:i + 1] == ' ':
            i = i + 1
        # For 'Cl', 'Br', etc.
        elif smiles[i:i + 2] in keys:
            seq.append(char_dict[smiles[i:i + 2]])
            i = i + 2
        elif smiles[i:i + 1] in keys:
            seq.append(char_dict[smiles[i:i + 1]])
            i = i + 1
        else:
            raise ValueError('character not found in dict')

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def canonize_reaction (mix,num):
    return '>'.join([canonize_mixture(sm,num) for sm in mix.split('>')])

#def canonize_cgr(mix,num):
#    return ' '.join([canonize_reaction(sm,num) for sm in mix.split(' ')])

#def canonize_reactants (mix,num):
#    return '~'.join([canonize_cgr(sm,num) for sm in mix.split('~')])

def canonize_mixture (mix,num):
    smiles = [canonize_smile(sm,num) for sm in mix.split('.')]
    if num > 0 and len(smiles)>1:
        random.shuffle(smiles)
    return '.'.join(smiles)

def canonize_record (mix,num):
    n = 0
    for m in mix.split(","):
        if n == 0:
            mol= canonize_reaction(m,num)
            n += 1
        else:
            mol += ","+canonize_reaction(m,num) # num -- only first is canonical; 0 -- all are canonical ones
    return mol

def canonize_smile (sm,num):
    if len(sm)<1:
        return ""
    
    if "error" in sm:
        raise Exception("error in smiles " + sm)
    
    m = Chem.MolFromSmiles(sm, sanitize = sanitize)

    if (methodtype == "EAGCNG" or methodtype == "PYTORCH" or methodtype == "ATTFP") and len(m.GetBonds())<1:
        raise ValueError('no bonds')
    
    if m.GetNumAtoms() <= 1:
        return sm
    
    if m is  None:
        raise ValueError('failed to produce molecule from ' + sm)
    
    if sanitize:
        AllChem.GetMorganFingerprintAsBitVect(m, 1, nBits=10)
        
    if num == 0 and sanitize:
        smi = Chem.MolToSmiles(m, canonical=True, isomericSmiles=isomeric)
    else:
        smi = Chem.MolToSmiles(m, canonical = False, isomericSmiles=isomeric, doRandom = True)
    
    if methodtype is not None and "TEXTCNN" in methodtype:
        length = methodtype.replace("TEXTCNN","")
        if len(length) == 0:
            length = 100000
        else:
            length = int(length)
        smiles_to_seq (smi, default_dict,length)

    if methodtype == "ATTFP":
        gen_data(smi)

    if methodtype == "DIMENET":
        if "#" in smi:
            raise ValueError('triple bond # is not supported.')

    #if methodtype == "CNF2":
    #    if len(smi)<2:
    #         raise ValueError('should have at least two characters')

    #if fixstereo:
    #    smi = SMILESToChiSMI(smi)

    if methodtype == "DLCA":
        for c in list(smi):
            if c not in allowedSymbols: 
                raise ValueError('not allowed symbol ' + c)

    if methodtype == 'HAMNET':
        sys.path.append(os.path.join(os.path.dirname(__file__), dir))
        from data.encode import encode_smiles
        encode_smiles([smi])


    if  'SMILESX' in methodtype or  'CHEMNLP' in methodtype:
        SmilesOK(smi, 99999)

    if  'PAINN' in methodtype:
        try: mg
        except NameError: 
            from kgcnn.mol.enocder import MolecularGraphRDKit
            mg = MolecularGraphRDKit(add_hydrogen=True, make_conformers=True)
        mg.from_smiles(smi, sanitize=True)
        return smi
        
    if "KERAS" in methodtype or "KGCNN" in methodtype:
        SmilesOK(smi, 99999)
        try:
            mol = MolFromSmiles(smi,sanitize = sanitize)
            if mol.GetNumAtoms() < 2:
                raise ValueError('one atom only')
        except:
            raise ValueError('no atoms in: ' + smi)
            
    return smi

def SmilesOK(smi, maxlen):
    try:
        mol = MolToSmiles(MolFromSmiles(smi,sanitize = sanitize))
    except:
        raise ValueError('SMILES cannot be parsed: ' + smi)
    if len(smi) >= maxlen: # 10% of margin
            raise ValueError('too large molecule')

def gen_data(smiles):
        sys.path.append(os.path.join(os.path.dirname(__file__), dir))
        from AttentiveFP.getFeatures import gen_descriptor_data
        res=gen_descriptor_data(smiles.split("\n"))
        if len(res) == 0:
            raise ValueError('AttentiveFP Features failed')

def str_to_bool(s):
    if s == 'True':
        return True
    elif s == 'False':
        return False
    else:
        raise ValueError

def filter_smiles (infile, outfile, augment):
    with open(infile, 'r') as fp:
        with open(outfile, 'w') as out_fh:
            writer = csv.writer(out_fh, quoting=csv.QUOTE_MINIMAL)
            line = "start"
            while line:
                line = fp.readline().strip()
                iniline = line
                if len(line) == 0:
                    continue                
                for i in range(augment):
                    try:
                        s = canonize_record(line,i)
                        writer.writerow([iniline, s])
                    except Exception as e:
                        writer.writerow([iniline, "ERROR " + str(e)])

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', '-i', help='input file', required=True)
    parser.add_argument('--outfile', '-o', help='output file', required=True)
    parser.add_argument('--augment', '-a', help='augment', required=True)
    parser.add_argument('--isomeric', '-s', help='isomeric', required=True)
    parser.add_argument('--nosanitize', '-n', help='nosanitize', required=False) # by default it is used
    parser.add_argument('--methodtype', '-t', help='Method specific corrections', required=False)
    #parser.add_argument('--fixstereo', '-f', help='Fix Stereo', required=True)
    args = parser.parse_args()
    isomeric = str_to_bool(args.isomeric)
    methodtype = args.methodtype
    if methodtype is not None:
        methodtype = methodtype.upper()
    else:
        methodtype = "NONE"

    sanitize = args.nosanitize is None
    #fixstereo = str_to_bool(args.fixstereo)
    dir = pathlib.Path(args.infile).parent.absolute()

    filter_smiles(args.infile, args.outfile,int(args.augment))
