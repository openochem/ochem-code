from rdkit import Chem
import numpy as np
import pandas as pd

import osmordred as rd

# Define descriptor computation function

def CalcOsmordred(mol, names = False, mynames=[]):

    " expended descriptors with more features and fixed InformationContent in cpp"
    v = 2
    doExEstate = True

    results = []
    descriptor_names = []

    # Define descriptors with names
    descriptors = [
        ("ABCIndex", rd.CalcABCIndex),
        ("AcidBase", rd.CalcAcidBase),
        ("AdjacencyMatrix", lambda mol: rd.CalcAdjacencyMatrix(mol, v)),
        ("Aromatic", rd.CalcAromatic),
        ("AtomCount", lambda mol: rd.CalcAtomCount(mol, v)),
        ("Autocorrelation", rd.CalcAutocorrelation),
        ("BCUT", rd.CalcBCUT),
        ("BalabanJ", rd.CalcBalabanJ),
        ("BaryszMatrix", rd.CalcBaryszMatrix),
        ("BertzCT", rd.CalcBertzCT),
        ("BondCount", rd.CalcBondCount),
        ("RNCGRPCG", rd.CalcRNCGRPCG),
        ("CarbonTypes", lambda mol: rd.CalcCarbonTypes(mol, v)),
        ("Chi", rd.CalcChi),
        ("Constitutional", rd.CalcConstitutional),
        ("DetourMatrix", rd.CalcDetourMatrix),
        ("DistanceMatrix", lambda mol: rd.CalcDistanceMatrix(mol, v)),
        ("EState", lambda mol: rd.CalcEState(mol, doExEstate)),
        ("EccentricConnectivityIndex", rd.CalcEccentricConnectivityIndex),
        ("ExtendedTopochemicalAtom", rd.CalcExtendedTopochemicalAtom),
        ("FragmentComplexity", rd.CalcFragmentComplexity),
        ("Framework", rd.CalcFramework),
        ("HydrogenBond", rd.CalcHydrogenBond),
    ]


    descriptors.append(("LogS", rd.CalcLogS))
    descriptors.append(("InformationContentv2", lambda mol: rd.CalcInformationContent(mol, 5)))

    additional_descriptors = [
        ("KappaShapeIndex", rd.CalcKappaShapeIndex),
        ("Lipinski", rd.CalcLipinski),
        ("McGowanVolume", rd.CalcMcGowanVolume),
        ("MoeType", rd.CalcMoeType),
        ("MolecularDistanceEdge", rd.CalcMolecularDistanceEdge),
        ("MolecularId", rd.CalcMolecularId),
        ("PathCount", rd.CalcPathCount),
        ("Polarizability", rd.CalcPolarizability),
        ("RingCount", rd.CalcRingCount),
        ("RotatableBond", rd.CalcRotatableBond),
        ("SLogP", rd.CalcSLogP),
        ("TopoPSA", rd.CalcTopoPSA),
        ("TopologicalCharge", rd.CalcTopologicalCharge),
        ("TopologicalIndex", rd.CalcTopologicalIndex),
        ("VdwVolumeABC", rd.CalcVdwVolumeABC),
        ("VertexAdjacencyInformation", rd.CalcVertexAdjacencyInformation),
        ("WalkCount", rd.CalcWalkCount),
        ("Weight", rd.CalcWeight),
        ("WienerIndex", rd.CalcWienerIndex),
        ("ZagrebIndex", rd.CalcZagrebIndex),
    ]

    descriptors.extend(additional_descriptors)

    extended_descriptors = [
        ("Pol", rd.CalcPol),
        ("MR", rd.CalcMR),
        ("Flexibility", rd.CalcFlexibility),
        ("Schultz", rd.CalcSchultz),
        ("AlphaKappaShapeIndex", rd.CalcAlphaKappaShapeIndex),
        ("HEState", rd.CalcHEState),
        ("BEState", rd.CalcBEState),
        ("Abrahams", rd.CalcAbrahams),
        ("ANMat", rd.CalcANMat),
        #("ASMat", rd.CalcASMat),
        ("AZMat", rd.CalcAZMat),
        #("DSMat", rd.CalcDSMat),
        ("DN2Mat", rd.CalcDN2Mat),
        ("Frags", rd.CalcFrags),
        ("AddFeatures", rd.CalcAddFeatures),
    ]
    descriptors.extend(extended_descriptors)
    for name, func in descriptors:
        value = func(mol)
        value = np.atleast_1d(np.array(value))
        results.append(value)
        if names:
            descriptor_names.extend([f"{name}_{i+1}" for i in range(len(value))])
    if names:
        return np.concatenate(results), descriptor_names
    return np.concatenate(results)


if __name__ == "__main__":
    _, mynames = CalcOsmordred(Chem.MolFromSmiles('C=C'),names=True)

    results = []
    
    with Chem.SDMolSupplier('datain') as suppl:
        for mol in suppl:
            if mol is None:
                results.append(['error: reading the molecule'])
                continue
            try:
                results.append(CalcOsmordred(mol,names=False))
            except Exception as e:
                results.append(['error:' + str(e)])
                continue
        
    # Convert to DataFrame and save to file
    df_results = pd.DataFrame(results, columns=mynames)
    
    #df_results = df_results.drop(['ASMat_4', 'ASMat_9', 'DSMat_4', 'DSMat_19', 'ANMat_14', 'Chi_41', 'Chi_49','AZMat_10'], axis=1)
    df_results = df_results.drop(['ANMat_14', 'Chi_41', 'Chi_49', 'AZMat_10'], axis=1)

    print(f"Finished processing. Results shape: {df_results.shape}")
    df_results.to_csv('dataout')
    print("youppy saved!")
