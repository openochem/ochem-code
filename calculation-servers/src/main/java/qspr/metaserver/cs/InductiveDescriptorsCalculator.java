
/* Copyright (C) 2022 BIGCHEM GmbH <info@bigchem.de>
 *
 * Contact: info@bigchem.de
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License (AGPL)
 * as published by the Free Software Foundation; either version 3.0
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the Affero GNU General Public License for more details.
 *
 * You should have received a copy of the Affero GNU Lesser General Public License
 * along with this program; If not, see <https://www.gnu.org/licenses/>. 
 */

package qspr.metaserver.cs;

import static java.lang.String.format;
import static java.util.Arrays.fill;

/**
 * Calculation of 'inductive' descriptors.
 * 
 * Based on an MOE script and several publications by Artem Cherkasov:
 * 1) A. Cherkasov, Inductive Electronegativity Scale. Iterative Calculation of Inductive Partial Charges, J.Chem.Inf.Comput.Sci. 43: 2039-2047, 2003.
 * 2) A. Cherkasov, D. Sprous, R. Chen, Three-Dimensional Correlation Analysis - A Novel Approach to the Quantification of Substituent Effects, J.Phys.Chem. 107: 9695-9704, 2003.
 * 3) A. Cherkasov, 'Inductive' Descriptors: 10 Successful Years in QSAR, Curr.Comput.AidedDrug.Des. 1: 21-42, 2005.
 * The calculated descriptors are defined in reference 3). 8 more descriptors were introduced due to discrepancies between published formulas and values computed by the MOE script.
 * The descriptors are basically summary statistics of per-atom descriptors like charge, electronegativity, etc.
 * 
 * @author mrupp
 */
class InductiveDescriptorsCalculator
{
	// Original method constants.
	protected static final float MAX_INTERACTION_DISTANCE = 30;  /** Maximum distance (in Angstrom) up to which interactions are considered. */
	protected static final float ELECTRONEGATIVITY_TERMINATION_CONDITION = 0.005f;  /** Termination condition (smallest delta for which iteration continues) for electronegativity/partial charges computation. */

	/** 
	 * Computed descriptors.
	 * 
	 * These 50 descriptors were introduced in 
	 * A. Cherkasov, Curr.Comput.AidedDrug.Des. 1: 21-42, 2005.
	 * Based partially on MOE SVL code by Artem Cherkasov.
	 * 8 more descriptors were introduced due to discrepancies between published formulas and values computed by the script.
	 * 
	 * The arrays are per-atom descriptors, the variables are per-molecule descriptors.
	 */
	static class Result
	{
		//AtomsCount
		int atomsCount;

		// Electronegativity-related.
		float[] aElectronegativity = null;
		float mEoEqualized  = 0;
		float mAverageEoPos = 0;
		float mAverageEoNeg = 0;

		// Hardness-related.
		float[] aHardness = null;
		float mGlobalHardness      = 0;
		float mSumHardness         = 0;
		float mSumPosHardness      = 0;
		float mSumNegHardness      = 0;
		float mAverageHardness     = 0;
		float mAveragePosHardness  = 0;
		float mAverageNegHardness  = 0;
		float mSmallestPosHardness = 0;
		float mSmallestNegHardness = 0;
		float mLargestPosHardness  = 0;
		float mLargestNegHardness  = 0;
		float mHardnessOfMostPos   = 0;
		float mHardnessOfMostNeg   = 0;

		// Softness-related.
		float[] aSoftness = null;
		float mGlobalSoftness      = 0;
		float mTotalPosSoftness    = 0;
		float mTotalNegSoftness    = 0;
		float mAverageSoftness     = 0;
		float mAveragePosSoftness  = 0;
		float mAverageNegSoftness  = 0;
		float mSmallestPosSoftness = 0;
		float mSmallestNegSoftness = 0;
		float mLargestPosSoftness  = 0;
		float mLargestNegSoftness  = 0;
		float mSoftnessOfMostPos   = 0;
		float mSoftnessOfMostNeg   = 0;

		// Charge-related.
		float[] aPartialCharge   = null;
		float mTotalCharge       = 0;
		float mTotalChargeFormal = 0;
		float mAveragePosCharge  = 0;
		float mAverageNegCharge  = 0;
		float mMostPosCharge     = 0;
		float mMostNegCharge     = 0;

		// Sigma-related (inductive parameter)
		float[] aSigmaMolI       = null;
		float[] aSigmaIMol       = null;
		float mTotalSigmaMolI    = 0;
		float mTotalAbsSigmaMolI = 0;
		float mMostPosSigmaMolI  = 0;
		float mMostNegSigmaMolI  = 0;
		float mMostPosSigmaIMol  = 0;
		float mMostNegSigmaIMol  = 0;
		float mSumPosSigmaMolI   = 0;
		float mSumNegSigmaMolI   = 0;

		// Rs-based (steric parameter)
		float[] aRsMolI       = null;
		float[] aRsIMol       = null;
		float mLargestRsMolI  = 0;
		float mSmallestRsMolI = 0;
		float mLargestRsIMol  = 0;
		float mSmallestRsIMol = 0;
		float mMostPosRsMolI  = 0;
		float mMostNegRsMolI  = 0;
		float mMostPosRsIMol  = 0;
		float mMostNegRsIMol  = 0;

		// These four descriptors originally were buggy implementations of mMostPosRsMolI, mMostNegRsMolI, mMostPosRsIMol, mMostNegRsIMol.
		// in Artem Cherkasov's SVL script. They are provided as additional descriptors after consultation with Artem.		
		float mLargestPosRsMolI  = 0;  
		float mSmallestNegRsMolI = 0;
		float mLargestPosRsIMol  = 0;
		float mSmallestNegRsIMol = 0;

		/** Initializes data structure with zeros. */
		Result(final int numAtoms)
		{
			aElectronegativity = new float[numAtoms]; fill(aElectronegativity, 0);
			aHardness          = new float[numAtoms]; fill(aHardness         , 0);
			aSoftness          = new float[numAtoms]; fill(aSoftness         , 0);
			aPartialCharge     = new float[numAtoms]; fill(aPartialCharge    , 0);
			aSigmaMolI         = new float[numAtoms]; fill(aSigmaMolI        , 0);
			aSigmaIMol         = new float[numAtoms]; fill(aSigmaIMol        , 0);
			aRsMolI            = new float[numAtoms]; fill(aRsMolI           , 0);
			aRsIMol            = new float[numAtoms]; fill(aRsIMol           , 0);
		}
	}

	/**
	 * Holds information parsed from a molecule in SDF format and derived quantities.
	 * 
	 * Derived quantities include the number of lone electron pairs and atom types.
	 */
	protected static class MoleculeData
	{
		final int numAtoms;         // Number of atoms.
		final int numBonds;         // Number of covalent bonds.
		final Element[] elems;      // element types of atoms (index is order in sdf).
		final AtomType[] types;     // Atom types (index is order in sdf).
		final float[][] coords;     // 3d coordinates of atoms (fist index is order in sdf, second is dimension).
		final int[] bondsA;         // First atom index of bond.
		final int[] bondsB;         // Second atom index of bond.
		final int[] bondsType;      // 1 = single, 2 = double, 3 = triple, 4 = aromatic.
		final int[][] connectivity; // Same interpretation as bondsType.
		final int[] formalCharges;  // Formal charges.
		final int[] lonePairs;      // Number of lone pairs = number of electrons in lone pairs / 2.

		/** 
		 * Initialies from a molecule in SDF format.
		 * 
		 * If an error occurs, an exception is thrown.
		 * The SDF format is defined in the MDL/Symyx document "CTFile Formats" from November 2007.
		 * 
		 * @param sdf Input molecule in SDF format, as a string.
		 * @throws java.io.IOException In case of errors in the SDF string.
		 * @throws ComputationFailureException In case of an error in processing the compound.
		 */
		MoleculeData(final String sdf) throws java.io.IOException, ComputationFailureException
		{
			//
			// Parse data from SDF string.
			final java.io.BufferedReader reader = new java.io.BufferedReader(new java.io.StringReader(sdf));
			reader.readLine(); reader.readLine(); reader.readLine();  // Skip header lines.

			// Parse counts line.
			String line = reader.readLine();
			numAtoms = Integer.parseInt(line.substring(0, 3).trim());
			numBonds = Integer.parseInt(line.substring(3, 7).trim());

			if(numAtoms > InductiveDescriptorsServer.MAXATOMS)
				throw new ComputationFailureException("number of atoms more then " + InductiveDescriptorsServer.MAXATOMS, sdf);
			// Parse atom information.
			elems         = new Element[numAtoms];
			coords        = new float[numAtoms][3];
			formalCharges = new int[numAtoms];

			try{
				for(int i = 0; i < numAtoms; ++i) 
				{
					final String s = reader.readLine();
					coords[i][0] = Float.parseFloat(s.substring( 0,10).trim());
					coords[i][1] = Float.parseFloat(s.substring(10,20).trim());
					coords[i][2] = Float.parseFloat(s.substring(20,30).trim());
					elems[i] = Element.fromString(s.substring(31,34).trim());
					switch(Integer.parseInt(s.substring(36,39).trim()))
					{
					case 0: formalCharges[i] =  0; break;
					case 1: formalCharges[i] = +3; break;
					case 2: formalCharges[i] = +2; break;
					case 3: formalCharges[i] = +1; break;
					case 4: throw new ComputationFailureException(format("Atom typing failed - Doublet radical charge on atom %d not handled", i+1));
					case 5: formalCharges[i] = -1; break;
					case 6: formalCharges[i] = -2; break;
					case 7: formalCharges[i] = -3; break;
					default: throw new ComputationFailureException(String.format("Atom typing failed - Unknown charge constant %d in SDF file", Integer.parseInt(s.substring(36,39).trim())), sdf);
					}
				}}catch(NumberFormatException e){
					throw new ComputationFailureException("Cannot understand format of this molecule ",e.getLocalizedMessage());
				}

			// Parse bond information.
			bondsA = new int[numBonds];
			bondsB = new int[numBonds];
			bondsType = new int[numBonds];
			for(int i = 0; i < numBonds; ++i)try
			{
				final String s = reader.readLine();
				bondsA[i]    = Integer.parseInt(s.substring(0,3).trim()) - 1;
				bondsB[i]    = Integer.parseInt(s.substring(3,6).trim()) - 1;
				bondsType[i] = Integer.parseInt(s.substring(7,9).trim());
				if(bondsType[i] == 4) throw new ComputationFailureException(format("Atom typing failed - Aromatic bond %d-%d not handled", bondsA[i]+1, bondsB[i]+1));  // TODO: How to handle aromatic bonds?
				if(bondsType[i] >  4) throw new ComputationFailureException(format("Atom typing failed - Ambivalent bond %d-%d type %d", bondsA[i]+1, bondsB[i]+1, bondsType[i]));
			}catch(NumberFormatException e){
				throw new ComputationFailureException("Error in calculation of bonds indices  " + e.getMessage());
			}
			catch(NullPointerException e){
				throw new ComputationFailureException("Error in calculation of bonds indices: null pointer exception");
			}

			connectivity = new int[numAtoms][numAtoms];
			for(int i = 0; i < numAtoms; ++i) fill(connectivity[i], 0);
			for(int i = 0; i < numBonds; ++i) connectivity[bondsA[i]][bondsB[i]] = connectivity[bondsB[i]][bondsA[i]] = bondsType[i];

			// Handle one-atom compounds.
			if(numAtoms == 1) throw new ComputationFailureException("'Inductive' descriptors are not defined for one-atom records.");

			//
			// Type atoms
			lonePairs = new int[numAtoms];      fill(lonePairs, Integer.MIN_VALUE);
			types     = new AtomType[numAtoms]; fill(types, null);
			final AtomType.Hybridization[] hybridizationStates = new AtomType.Hybridization[]{AtomType.Hybridization.NONE, AtomType.Hybridization.SP, AtomType.Hybridization.SP2, AtomType.Hybridization.SP3, AtomType.Hybridization.SP3D, AtomType.Hybridization.SP3D2};
			for(int ind = 0; ind < numAtoms; ++ind)
			{
				final int numNeighbors   = numNeighbors(ind);
				final int totalBondOrder = totalBondOrder(ind);

				// Determine the number of valence electrons
				int valence = -1;
				switch(elems[ind].periodicTableGroup())
				{
				case  1: valence = 1; break;
				case  2: valence = 2; break;
				case 13: valence = 3; break;
				case 14: valence = 4; break;
				case 15: valence = 5; break;
				case 16: valence = 6; break;
				case 17: valence = 7; break;
				case 18: if(elems[ind] == InductiveDescriptorsCalculator.Element.He) valence = 2; else valence = 8; break;
				case  3:
				case  4:
				case  5:
				case  6:
				case  7: 
				case  8:
				case  9:
				case 10:
				case 11:
				case 12:
					switch(elems[ind].periodicTablePeriod())
					{
					case 4:
					case 5:
						valence = elems[ind].periodicTableGroup(); break;
					case 6:
					case 7:  // TODO: How to define valency for lanthanoids and actinoids?
						if( (elems[ind].atomicNumber() >= 57 && elems[ind].atomicNumber() <=  71) || (elems[ind].atomicNumber() >= 89 && elems[ind].atomicNumber() <= 103) )
							valence = -999;  // This is obviously wrong, but will simply lead to UNKNOWN hybridization.
						valence = 14 + elems[ind].periodicTableGroup(); 
						break;
					default: 
						throw new ComputationFailureException(String.format("Atom typing failed - could not determine number of valence electrons for %s", elems[ind].abbrev()), sdf);
					}
					break;
				default: 
					throw new ComputationFailureException(String.format("Atom typing failed - could not determine number of valence electrons for %s", elems[ind].abbrev()), sdf);
				}

				// Determine number of lone pairs
				final int lonePairElectrons = valence - formalCharges[ind] - totalBondOrder;
				lonePairs[ind] = lonePairElectrons / 2;
				final boolean hybridizationFailure = (lonePairElectrons % 2 != 0) || (lonePairs[ind] < 0);  // Instead of aborting the calculation, we resort to default values.				

				// Sanity check
				if(!hybridizationFailure && formalCharges[ind] != valence - 2*lonePairs[ind] - totalBondOrder)  // Classic formula. Should be the same as valence - fullShell + totalBondOrder.
					throw new ComputationFailureException(String.format("Atom typing failed - internal consistency check failed for atom %d (%s)", ind+1, elems[ind].abbrev()), sdf);

				// Hybridization (see http://en.wikipedia.org/wiki/Orbital_hybridisation#Hybridisation_and_molecule_shape)
				AtomType.Hybridization hybridization = AtomType.Hybridization.UNKNOWN;
				if(numNeighbors == 0) hybridization = AtomType.Hybridization.ION; 
				else if(!hybridizationFailure)  // In case of failure to determine hybridization state, use UNKNOWN, leading to default values.
				{
					final int index = numNeighbors + lonePairs[ind] - 1;
					if(index < 0 || index >= hybridizationStates.length) hybridization = AtomType.Hybridization.UNKNOWN;
					else hybridization = hybridizationStates[index];
				}

				// Atom type.
				types[ind] = AtomType.atomType(elems[ind], hybridization, formalCharges[ind]);
			}
		}

		/** Returns the number of neighbors. */
		final int numNeighbors(final int i)
		{
			int result = 0;
			for(int j = 0; j < numAtoms; ++j) if(connectivity[i][j] != 0) ++result;
			return result;
		}

		/** Returns the sum of the bond orders of the given atom. */
		final int totalBondOrder(final int ind)
		{
			int result = 0;
			for(int j = 0; j < bondsA.length; ++j)
				if(bondsA[j] == ind || bondsB[j] == ind) result += bondsType[j];
			return result;
		}

		/** Returns the index of the j-th bond of atom i. */
		final int bond(final int i, final int j) throws ComputationFailureException
		{
			int k = 0;
			for(int l = 0; l < numBonds; ++l)
				if(bondsA[l] == i || bondsB[l] == i)
					if(++k == j) return l;
			throw new ComputationFailureException(String.format("Atom typing failed: Could not find %d-th neighbor of atom %d.", j, i));
		}

		/** Returns the Euclidean distance between the 3d coordinates of two atoms. */
		final float distance(final int i, final int j)
		{
			return (float) Math.sqrt(Math.pow(coords[i][0]-coords[j][0], 2) + Math.pow(coords[i][1]-coords[j][1], 2) + Math.pow(coords[i][2]-coords[j][2], 2));
		}
	}

	/**
	 * Calculates the descriptors for a given molecule in SDF format.
	 * 
	 * The input molecule must have optimized geometry and explicit hydrogens.
	 * 
	 * @param sdf Input molecule (must have optimized structure and explicit hydrogens) in SDF format as a string.
	 * @return Object containing the calculated descriptors.
	 */
	static Result calculate(final String sdf) throws java.io.IOException, ComputationFailureException
	{
		// Parse the molecule and calculate atom types.
		final MoleculeData md = new MoleculeData(sdf);
		Result r = new Result(md.numAtoms);

		//
		// Compute per-atom descriptors.

		// Set up parts depending on atomic radii, distances, preset electronegativities.
		final float[][] contrib = new float[md.numAtoms][md.numAtoms];  // Stores the (Ri^2 + Rj^2)/r_ij^2 terms. These occur in many formulas.
		final float[] inductiveCharge = new float[md.numAtoms];
		for(int i = 0; i < md.numAtoms; ++i)
		{
			// Electronegativity
			r.aElectronegativity[i] = md.types[i].inductiveElectronegativity();

			// Softness
			for(int j = 0; j < md.numAtoms; ++j)
			{
				if(j == i) continue;
				final float distance = md.distance(i, j); 
				final float squaredDistance = (float) Math.pow(distance, 2); 
				if(distance >= MAX_INTERACTION_DISTANCE) continue;

				contrib[i][j] = (md.types[i].squaredCovalentRadius() + md.types[j].squaredCovalentRadius()) / squaredDistance;
				final float eoDiff = (md.types[j].inductiveElectronegativity() - md.types[i].inductiveElectronegativity());  // Difference \chi_j - \chi_i in electronegativity.

				r.aSoftness[i]      += 2 * contrib[i][j];  // Equation 21.
				r.aPartialCharge[i] += eoDiff * contrib[i][j];  // Equation 17.
				r.aRsMolI[i]        += md.types[j].squaredCovalentRadius() / squaredDistance;
				r.aRsIMol[i]        += md.types[i].squaredCovalentRadius() / squaredDistance;
				r.aSigmaMolI[i]     += +eoDiff * md.types[j].squaredCovalentRadius() / squaredDistance;
				r.aSigmaIMol[i]     += -eoDiff * md.types[i].squaredCovalentRadius() / squaredDistance;
			}

			// Hardness
			r.aHardness[i] = 2 * 5.8139f / r.aSoftness[i];

			// Charge
			inductiveCharge[i] = r.aPartialCharge[i];
			r.aPartialCharge[i] = 0.172f * r.aPartialCharge[i] + md.formalCharges[i]; 
		}

		// Iterative quantities.
		float[] eoEqualized = new float[md.numAtoms];
		float eoGradient = 999;  // L1 metric of average and equalized electronegativity.
		int num = 10000;
		while(eoGradient > ELECTRONEGATIVITY_TERMINATION_CONDITION)
		{
			if(num-- == 0) throw new ComputationFailureException("Number of iterations > 10000 - cannot reach convergence");
			
			for(int i = 0; i < md.numAtoms; ++i)
			{
				eoEqualized[i] = md.types[i].inductiveElectronegativity() + inductiveCharge[i] / r.aSoftness[i];
			}

			final float eoEqualizedAverage = sum(eoEqualized) / eoEqualized.length;
			eoGradient = sumAbsDiff(eoEqualized, eoEqualizedAverage);

			for(int i = 0; i < eoEqualized.length; ++i)
				inductiveCharge[i] += sumDiffTimes(eoEqualized, eoEqualized[i], contrib[i]);
		}

		//
		// Compute per-molecule descriptors (summary statistics).
		r.mSmallestPosHardness = Float.POSITIVE_INFINITY; r.mSmallestPosSoftness  = Float.POSITIVE_INFINITY;
		r.mSmallestNegHardness = Float.POSITIVE_INFINITY; r.mSmallestNegSoftness  = Float.POSITIVE_INFINITY;
		r.mLargestPosHardness  = Float.NEGATIVE_INFINITY; r.mLargestPosSoftness   = Float.NEGATIVE_INFINITY;
		r.mLargestNegHardness  = Float.NEGATIVE_INFINITY; r.mLargestNegSoftness   = Float.NEGATIVE_INFINITY;
		r.mMostNegCharge       = Float.POSITIVE_INFINITY; r.mMostPosCharge        = Float.NEGATIVE_INFINITY;
		r.mMostNegSigmaMolI    = Float.POSITIVE_INFINITY; r.mMostPosSigmaMolI     = Float.NEGATIVE_INFINITY;
		r.mMostNegSigmaIMol    = Float.POSITIVE_INFINITY; r.mMostPosSigmaIMol     = Float.NEGATIVE_INFINITY;
		r.mLargestRsMolI       = Float.NEGATIVE_INFINITY; r.mLargestRsIMol        = Float.NEGATIVE_INFINITY;
		r.mSmallestRsMolI      = Float.POSITIVE_INFINITY; r.mSmallestRsIMol       = Float.POSITIVE_INFINITY;
		r.mMostNegRsMolI       = Float.POSITIVE_INFINITY; r.mMostPosRsMolI        = Float.NEGATIVE_INFINITY;
		r.mMostNegRsIMol       = Float.POSITIVE_INFINITY; r.mMostPosRsIMol        = Float.NEGATIVE_INFINITY;
		r.mLargestPosRsMolI    = Float.NEGATIVE_INFINITY; r.mSmallestNegRsMolI    = Float.POSITIVE_INFINITY;
		r.mLargestPosRsIMol    = Float.NEGATIVE_INFINITY; r.mSmallestNegRsIMol    = Float.POSITIVE_INFINITY;

		int numNegCharge = 0; int numPosCharge = 0;
		for(int i = 0; i < md.numAtoms; ++i)
		{
			final float softness  = r.aSoftness[i];
			final float hardness  = 1 / softness;
			final float charge    = r.aPartialCharge[i];
			final float sigmaMolI = r.aSigmaMolI[i];
			final float sigmaIMol = r.aSigmaIMol[i];
			final float rsMolI    = r.aRsMolI[i];
			final float rsIMol    = r.aRsIMol[i];

			// Unconditional
			r.mGlobalSoftness    += softness;
			r.mSumHardness       += r.aHardness[i]; // With constant.
			r.mTotalCharge       += Math.abs(charge);
			r.mTotalChargeFormal += md.formalCharges[i];
			r.mTotalSigmaMolI    += sigmaMolI;
			r.mTotalAbsSigmaMolI += Math.abs(sigmaMolI);

			// Negative partial charge
			if(charge < 0)
			{
				++numNegCharge;
				if(charge < r.mMostNegCharge)
				{
					r.mMostNegCharge     = charge;
					r.mHardnessOfMostNeg = hardness; // Without constant.
					r.mSoftnessOfMostNeg = softness;
					r.mMostNegRsMolI     = rsMolI;
					r.mMostNegRsIMol     = rsIMol;
				}

				r.mAverageEoNeg += md.types[i].inductiveElectronegativity();

				r.mSumNegHardness     += hardness;  // Without the constant 2*5.8139.
				r.mAverageNegHardness += hardness;  // Without the constant 2*5.8139.
				if(hardness < r.mSmallestNegHardness) r.mSmallestNegHardness = hardness;
				if(hardness > r.mLargestNegHardness ) r.mLargestNegHardness  = hardness;

				r.mTotalNegSoftness   += softness;
				r.mAverageNegSoftness += softness;
				if(softness < r.mSmallestNegSoftness) r.mSmallestNegSoftness = softness;
				if(softness > r.mLargestNegSoftness ) r.mLargestNegSoftness  = softness;

				r.mAverageNegCharge += charge;  // Artem Cherkasov's code has this descriptor with positive sign.

				if(rsMolI < r.mSmallestNegRsMolI) r.mSmallestNegRsMolI = rsMolI;
				if(rsIMol < r.mSmallestNegRsIMol) r.mSmallestNegRsIMol = rsIMol;
			}
			// Non-negative partial charge
			else 
			{
				++numPosCharge;
				if(charge > r.mMostPosCharge)
				{
					r.mMostPosCharge     = charge;
					r.mHardnessOfMostPos = hardness; // Without constant.
					r.mSoftnessOfMostPos = softness;
					r.mMostPosRsMolI     = rsMolI;
					r.mMostPosRsIMol     = rsIMol;
				}

				r.mAverageEoPos += md.types[i].inductiveElectronegativity();

				r.mSumPosHardness += 1 / softness;  // Without the constant 2*5.8139.
				r.mAveragePosHardness += 1 / softness;  // Without the constant 2*5.8139.
				if(hardness < r.mSmallestPosHardness) r.mSmallestPosHardness = hardness;
				if(hardness > r.mLargestPosHardness ) r.mLargestPosHardness  = hardness;

				r.mTotalPosSoftness   += softness;
				r.mAveragePosSoftness += softness;
				if(softness < r.mSmallestPosSoftness) r.mSmallestPosSoftness = softness;
				if(softness > r.mLargestPosSoftness ) r.mLargestPosSoftness  = softness;

				r.mAveragePosCharge += charge;

				if(rsMolI > r.mLargestPosRsMolI) r.mLargestPosRsMolI = rsMolI;
				if(rsIMol > r.mLargestPosRsIMol) r.mLargestPosRsIMol = rsIMol;
			}

			// Sigma mol -> i
			if(sigmaMolI < 0) 
			{ 
				if(sigmaMolI < r.mMostNegSigmaMolI) r.mMostNegSigmaMolI = sigmaMolI;
				r.mSumNegSigmaMolI += sigmaMolI;
			}
			else 
			{
				if(sigmaMolI > r.mMostPosSigmaMolI) r.mMostPosSigmaMolI = sigmaMolI;
				r.mSumPosSigmaMolI += sigmaMolI;
			}

			// Sigma i -> mol
			if(sigmaIMol < 0) { if(sigmaIMol < r.mMostNegSigmaIMol) r.mMostNegSigmaIMol = sigmaIMol; }
			else if(sigmaIMol > r.mMostPosSigmaIMol) r.mMostPosSigmaIMol = sigmaIMol;

			// Rs
			if(rsMolI > r.mLargestRsMolI ) r.mLargestRsMolI  = rsMolI;
			if(rsMolI < r.mSmallestRsMolI) r.mSmallestRsMolI = rsMolI;
			if(rsIMol > r.mLargestRsIMol ) r.mLargestRsIMol  = rsIMol;
			if(rsIMol < r.mSmallestRsIMol) r.mSmallestRsIMol = rsIMol;
		}

		r.mEoEqualized   = eoEqualized[0];
		r.mAverageEoNeg /= numNegCharge;
		r.mAverageEoPos /= numPosCharge;

		r.mGlobalHardness      = 2 * 5.8139f / r.mGlobalSoftness;
		r.mAverageHardness     = (2 * 5.8139f / r.mGlobalSoftness) / md.numAtoms;
		r.mAveragePosHardness /= numPosCharge;
		r.mAverageNegHardness /= numNegCharge;

		r.mAverageSoftness     = r.mGlobalSoftness / md.numAtoms;
		r.mAveragePosSoftness /= numPosCharge;
		r.mAverageNegSoftness /= numNegCharge;

		r.mAverageNegCharge /= numNegCharge;
		r.mAveragePosCharge /= numPosCharge;

		// Check whether any of the computed values is infinity or not-a-number.
		final float check = 
				sum(r.aElectronegativity) + r.mEoEqualized       + r.mAverageEoPos      + r.mAverageEoNeg      +
				sum(r.aHardness         ) + r.mGlobalHardness    + r.mSumHardness       + r.mSumPosHardness    + r.mSumNegHardness   + r.mAverageHardness    + r.mAveragePosHardness + r.mAverageNegHardness  + r.mSmallestPosHardness + r.mSmallestNegHardness + r.mLargestPosHardness + r.mLargestNegHardness + r.mHardnessOfMostPos + r.mHardnessOfMostNeg +
				sum(r.aSoftness         ) + r.mGlobalSoftness    + r.mTotalPosSoftness  + r.mTotalNegSoftness  + r.mAverageSoftness  + r.mAveragePosSoftness + r.mAverageNegSoftness + r.mSmallestPosSoftness + r.mSmallestNegSoftness + r.mLargestPosSoftness  + r.mLargestNegSoftness + r.mSoftnessOfMostPos  + r.mSoftnessOfMostNeg +
				sum(r.aPartialCharge    ) + r.mTotalCharge       + r.mTotalChargeFormal + r.mAveragePosCharge  + r.mAverageNegCharge + r.mMostPosCharge      + r.mMostNegCharge      +
				sum(r.aSigmaMolI        ) + sum(r.aSigmaIMol)    + r.mTotalSigmaMolI    + r.mTotalAbsSigmaMolI + r.mMostPosSigmaMolI + r.mMostNegSigmaMolI   + r.mMostPosSigmaIMol   + r.mMostNegSigmaIMol    + r.mSumPosSigmaMolI     + r.mSumNegSigmaMolI  +
				sum(r.aRsMolI           ) + sum(r.aRsIMol   )    + r.mLargestRsMolI     + r.mSmallestRsMolI    + r.mLargestRsIMol    + r.mSmallestRsIMol     + r.mMostPosRsMolI      + r.mMostNegRsMolI       + r.mMostPosRsIMol       + r.mMostNegRsIMol +
				r.mLargestPosRsMolI       + r.mSmallestNegRsMolI + r.mLargestPosRsIMol  + r.mSmallestNegRsIMol;			
		if(Float.isInfinite(check) || Float.isNaN(check)) throw new ComputationFailureException(String.format("Descriptor computation failed - at least one descriptor has infinite value or is not a number", sdf));

		return r;
	}	

	/** 
	 * Indicates that descriptors could not be computed.
	 * 
	 *  The SDF that caused the problem can be passed along with the descriptive message.
	 *  In that case, the SDF might be logged for debugging.
	 */
	static class ComputationFailureException extends Exception
	{
		private static final long serialVersionUID = 7153380463239646859L;
		final String sdf;

		ComputationFailureException(final String msg) { super(msg); sdf = null; }

		ComputationFailureException(final String msg, final String sdf) { super(msg); this.sdf = sdf; }
	}

	////////////////////////
	//  Helper functions  //
	////////////////////////

	// Sum of a vector.
	private static final float sum(final float[] v) 
	{
		float result = 0;
		for(final float x : v) result += x;
		return result;
	}

	// Sum of absolute differences between a vector and a scalar.
	private static final float sumAbsDiff(final float[] a, final float b)
	{
		float result = 0;
		for(int i = 0; i < a.length; ++i) result += Math.abs(a[i] - b);
		return result;
	}

	// Sum of differences between a vector and a scalar, multiplied component-wise with another vector.
	private static final float sumDiffTimes(final float[] v, final float s, final float[] w)
	{
		float result = 0;
		for(int i = 0; i < v.length; ++i)
			result += (v[i] - s) * w[i];
		return result;
	}

	/**
	 * Periodic table of the elements.
	 * 
	 * Elements with temporary names (three letter abbreviations) are not included.
	 * Provides abbreviation, periodic table group, and periodic table period.
	 */
	private static enum Element {
		H   (  1, "H" , "Hydrogen"     , 1,  1),
		He  (  2, "He", "Helium"       , 1, 18),
		Li  (  3, "Li", "Lithium"      , 2,  1),
		Be  (  4, "Be", "Beryllium"    , 2,  2),
		B   (  5, "B" , "Boron"        , 2, 13),
		C   (  6, "C" , "Carbon"       , 2, 14),
		N   (  7, "N" , "Nitrogen"     , 2, 15),
		O   (  8, "O" , "Oxygen"       , 2, 16),
		F   (  9, "F" , "Fluorine"     , 2, 17),
		Ne  ( 10, "Ne", "Neon"         , 2, 18),
		Na  ( 11, "Na", "Sodium"       , 3,  1),
		Mg  ( 12, "Mg", "Magnesium"    , 3,  2),
		Al  ( 13, "Al", "Aluminimum"   , 3, 13),
		Si  ( 14, "Si", "Silicon"      , 3, 14),
		P   ( 15, "P" , "Phosphorus"   , 3, 15),
		S   ( 16, "S" , "Sulfur"       , 3, 16),
		Cl  ( 17, "Cl", "Chlorine"     , 3, 17),
		Ar  ( 18, "Ar", "Argon"        , 3, 18),
		K   ( 19, "K" , "Potassium"    , 4,  1),
		Ca  ( 20, "Ca", "Calcium"      , 4,  2),
		Sc  ( 21, "Sc", "Scandium"     , 4,  3),
		Ti  ( 22, "Ti", "Titanium"     , 4,  4),
		V   ( 23, "V" , "Vanadium"     , 4,  5),
		Cr  ( 24, "Cr", "Chromium"     , 4,  6),
		Mn  ( 25, "Mn", "Manganese"    , 4,  7),
		Fe  ( 26, "Fe", "Iron"         , 4,  8),
		Co  ( 27, "Co", "Cobalt"       , 4,  9),
		Ni  ( 28, "Ni", "Nickel"       , 4, 10),
		Cu  ( 29, "Cu", "Copper"       , 4, 11),
		Zn  ( 30, "Zn", "Zinc"         , 4, 12),
		Ga  ( 31, "Ga", "Gallium"      , 4, 13),
		Ge  ( 32, "Ge", "Germanium"    , 4, 14),
		As  ( 33, "As", "Arsenic"      , 4, 15),
		Se  ( 34, "Se", "Selenium"     , 4, 16),
		Br  ( 35, "Br", "Bromine"      , 4, 17),
		Kr  ( 36, "Kr", "Krypton"      , 4, 18),
		Rb  ( 37, "Rb", "Rubidium"     , 5,  1),
		Sr  ( 38, "Sr", "Strontium"    , 5,  2),
		Y   ( 39, "Y" , "Yttrium"      , 5,  3),
		Zr  ( 40, "Zr", "Zirconium"    , 5,  4),
		Nb  ( 41, "Nb", "Niobium"      , 5,  5),
		Mo  ( 42, "Mo", "Molybdenum"   , 5,  6),
		Tc  ( 43, "Tc", "Technetium"   , 5,  7),
		Ru  ( 44, "Ru", "Ruthenium"    , 5,  8),
		Rh  ( 45, "Rh", "Rhodium"      , 5,  9),
		Pd  ( 46, "Pd", "Palladium"    , 5, 10),
		Ag  ( 47, "Ag", "Silver"       , 5, 11),
		Cd  ( 48, "Cd", "Cadmium"      , 5, 12),
		In  ( 49, "In", "Indium"       , 5, 13),
		Sn  ( 50, "Sn", "Tin"          , 5, 14),
		Sb  ( 51, "Sb", "Antimony"     , 5, 15),
		Te  ( 52, "Te", "Tellurium"    , 5, 16),
		I   ( 53, "I" , "Iodine"       , 5, 17),
		Xe  ( 54, "Xe", "Xenon"        , 5, 18),
		Cs  ( 55, "Cs", "Cesium"       , 6,  1),
		Ba  ( 56, "Ba", "Barium"       , 6,  2),
		La  ( 57, "La", "Lanthanum"    , 6,  3),
		Ce  ( 58, "Ce", "Cerium"       , 6,  3),
		Pr  ( 59, "Pr", "Praseodymium" , 6,  3),
		Nd  ( 60, "Nd", "Neodymium"    , 6,  3),
		Pm  ( 61, "Pm", "Promethium"   , 6,  3),
		Sm  ( 62, "Sm", "Samarium"     , 6,  3),
		Eu  ( 63, "Eu", "Europium"     , 6,  3),
		Gd  ( 64, "Gd", "Gadolinium"   , 6,  3),
		Tb  ( 65, "Tb", "Terbium"      , 6,  3),
		Dy  ( 66, "Dy", "Dysprosium"   , 6,  3),
		Ho  ( 67, "Ho", "Holmium"      , 6,  3),
		Er  ( 68, "Er", "Erbium"       , 6,  3),
		Tm  ( 69, "Tm", "Thulium"      , 6,  3),
		Yb  ( 70, "Yb", "Ytterbium"    , 6,  3),
		Lu  ( 71, "Lu", "Lutetium"     , 6,  3),
		Hf  ( 72, "Hf", "Hafnium"      , 6,  4),
		Ta  ( 73, "Ta", "Tantalum"     , 6,  5),
		W   ( 74, "W" , "Tungsten"     , 6,  6),
		Re  ( 75, "Re", "Rhenium"      , 6,  7),
		Os  ( 76, "Os", "Osmium"       , 6,  8),
		Ir  ( 77, "Ir", "Iridium"      , 6,  9),
		Pt  ( 78, "Pt", "Platinum"     , 6, 10),
		Au  ( 79, "Au", "Gold"         , 6, 11),
		Hg  ( 80, "Hg", "Mercury"      , 6, 12),
		Tl  ( 81, "Tl", "Thalium"      , 6, 13),
		Pb  ( 82, "Pb", "Lead"         , 6, 14),
		Bi  ( 83, "Bi", "Bismuth"      , 6, 15),
		Po  ( 84, "Po", "Polonium"     , 6, 16),
		At  ( 85, "At", "Astatine"     , 6, 17),
		Rn  ( 86, "Rn", "Radon"        , 6, 18),
		Fr  ( 87, "Fr", "Francium"     , 7,  1), 
		Ra  ( 88, "Ra", "Radium"       , 7,  2),
		Ac  ( 89, "Ac", "Actinium"     , 7,  3),
		Th  ( 90, "Th", "Thorium"      , 7,  3),
		Pa  ( 91, "Pa", "Protactinium" , 7,  3),
		U   ( 92, "U" , "Uranium"      , 7,  3),
		Np  ( 93, "Np", "Neptunium"    , 7,  3),
		Pu  ( 94, "Pu", "Plutonium"    , 7,  3),
		Am  ( 95, "Am", "Americium"    , 7,  3),
		Cm  ( 96, "Cm", "Curium"       , 7,  3),
		Bk  ( 97, "Bk", "Berkelium"    , 7,  3),
		Cf  ( 98, "Cf", "Californium"  , 7,  3),
		Es  ( 99, "Es", "Einsteinium"  , 7,  3),
		Fm  (100, "Fm", "Fermium"      , 7,  3),
		Md  (101, "Md", "Mendelevium"  , 7,  3),
		No  (102, "No", "Nobelium"     , 7,  3),
		Lr  (103, "Lr", "Lawrencium"   , 7,  3),
		Rf  (104, "Rf", "Rutherfordium", 7,  4),
		Db  (105, "Db", "Dubnium"      , 7,  5),
		Sg  (106, "Sg", "Seaborgium"   , 7,  6),
		Bh  (107, "Bh", "Bohrium"      , 7,  7),
		Hs  (108, "Hs", "Hassium"      , 7,  8),
		Mt  (109, "Mt", "Meitnerium"   , 7,  9),
		Ds  (110, "Ds", "Darmstadtium" , 7, 10),
		Rg  (111, "Rg", "Roentgenium"  , 7, 11),
		Cn  (112, "Cn", "Copernicium"  , 7, 12);

		/** Returns the atomic number of the element type. */
		int atomicNumber() { return atomicNumber; }

		/** Returns the official one to three characters abbreviation for the element type. */
		String abbrev() { return abbrev; }

		/** Returns the period of the element in the periodic table (IUPAC). */
		int periodicTablePeriod() { return ptPeriod; }

		/** Returns the group of the element in the periodic table (IUPAC). */
		int periodicTableGroup() { return ptGroup; }

		/**
		 * Converts from a textual representation to the enumeration. 
		 * 
		 * @param abbrev One to three characters abbreviation or full english name of the element.
		 * @return The enumeration type for the element.
		 */
		static Element fromString(String abbrev) throws ComputationFailureException {
			abbrev = abbrev.trim();
			for (Element element : Element.values()) if(element.abbrev.equals(abbrev)) return element;
			for (Element element : Element.values()) if(element.fullName.equals(abbrev)) return element;
			throw new ComputationFailureException(String.format("Unknown element abbreviation '%s' encountered in construction of Element from String.", abbrev));
		}

		// Internal representation.

		private final int atomicNumber; // Atomic number ("Ordnungszahl"), the number of protons.
		private final String abbrev;    // The two-character abbreviation of the element. 
		private final String fullName;  // Complete english name of the element.
		private final int ptPeriod;     // Periodic table (IUPAC) period.
		private final int ptGroup;      // Periodic table (IUPAC) group.

		Element(final int atomicNumber, final String abbrev, final String fullName, final int ptPeriod, final int ptGroup) 
		{ 
			this.atomicNumber = atomicNumber; 
			this.abbrev       = abbrev; 
			this.fullName     = fullName; 
			this.ptPeriod     = ptPeriod;
			this.ptGroup      = ptGroup;
		}
	}

	/**
	 * Atom types.
	 * 
	 * If nothing else is specified and the type was not contained in the original MOE script,
	 * electronegativity values are Pauling electronegativities (source: wikipedia), and, covalent bond radii and ion radii were also taken from wikipedia.
	 *
	 * The table is sorted by atomic number.
	 * Specific types come first. These are the types from the original script, as well as additional types added later on.
	 * Element, hybridization state, and formal charge must match.	    
	 * Then come fallback types These types are partial matches designed as fallbacks. If a more specific type matches first, the specific type will be taken.
	 * Entries should be in order of unspecificity (specific types with charges like ION, SP3, etc., then IONCA, then CA).
	 * 
	 * Table entries are element type, hybridization state, formal charge (in unit charges), electronegativity, squared covalent radius (in squared Angstrom).
	 * 
	 * Updated table received by Artem Cherkasov 2010-05-12.
	 * Whole periodic table added by Matthias Rupp 2010-06-09.
	 */
	private static enum AtomType {	
		H         ( Element.H , Hybridization.NONE   ,  0, 2.10f, 0.0900f), // Original entry from Artem Cherkasov's 'inductive' descriptors script
		//He_ca     ( Element.He, Hybridization.CA     ,  0, ?.??f, ?.????f), // Not defined as helium normally does not form covalent bonds (however, this is possible in principle, e.g., in intense laser fields) 
		Li_ionca  ( Element.Li, Hybridization.IONCA  ,  0, 0.98f, 0.5329f), // Values from wikipedia and webelements.
		Li_ca     ( Element.Li, Hybridization.CA     ,  0, 0.98f, 2.1025f),
		Be_ionca  ( Element.Be, Hybridization.IONCA  ,  0, 0.98f, 0.2025f),
		Be_ca     ( Element.Be, Hybridization.CA     ,  0, 1.57f, 1.1025f),
		B_sp2     ( Element.B , Hybridization.SP2    ,  0, 2.08f, 0.6561f),
		B_sp3_n   ( Element.B , Hybridization.SP3    , -1, 2.08f, 0.6561f), // Copied values.
		B_ion_p3  ( Element.B , Hybridization.ION    , +3, 2.04f, 0.0729f),
		B_ca      ( Element.B , Hybridization.CA     ,  0, 2.04f, 0.7225f),
		C_sp      ( Element.C , Hybridization.SP     ,  0, 3.15f, 0.3600f), 
		C_sp2_neg ( Element.C , Hybridization.SP2    , -1, 2.20f, 0.5929f), 
		C_sp2     ( Element.C , Hybridization.SP2    ,  0, 2.31f, 0.4489f), 
		C_sp2_pos ( Element.C , Hybridization.SP2    , +1, 2.31f, 0.4489f), 
		C_sp3     ( Element.C , Hybridization.SP3    ,  0, 2.20f, 0.5929f), 
		C_sp3_pos ( Element.C , Hybridization.SP3    , +1, 2.31f, 0.4489f),
		C_ca      ( Element.C , Hybridization.CA     ,  0, 2.55f, 0.4900f),
		N_sp      ( Element.N , Hybridization.SP     ,  0, 4.76f, 0.3025f), 
		N_sp_pos  ( Element.N , Hybridization.SP     , +1, 2.59f, 0.4900f), 
		N_sp2_neg ( Element.N , Hybridization.SP2    , -1, 2.49f, 0.4900f), // Value corrected 2010-05-12, from 2.59 to 2.49.
		N_sp2     ( Element.N , Hybridization.SP2    ,  0, 2.59f, 0.4900f),
		N_sp2_pos ( Element.N , Hybridization.SP2    , +1, 3.40f, 0.4900f),
		N_sp3_neg ( Element.N , Hybridization.SP3    , -1, 2.49f, 0.4900f),
		N_sp3     ( Element.N , Hybridization.SP3    ,  0, 2.59f, 0.4900f),
		N_sp3_pos ( Element.N , Hybridization.SP3    , +1, 3.39f, 0.4900f),
		N_ion_n3  ( Element.N , Hybridization.ION    , -3, 3.04f, 2.1316f),
		N_ion_p3  ( Element.N , Hybridization.ION    , +3, 3.04f, 0.0256f),
		N_ion_p5  ( Element.N , Hybridization.ION    , +5, 3.04f, 0.0169f),
		N_ca      ( Element.N , Hybridization.CA     ,  0, 3.04f, 0.4225f),
		O_sp2_neg ( Element.O , Hybridization.SP2    , -1, 1.80f, 0.4356f),
		O_sp2     ( Element.O , Hybridization.SP2    ,  0, 4.90f, 0.3844f),
		O_sp2_pos ( Element.O , Hybridization.SP2    , +1, 4.90f, 0.3844f), // Copied values.
		O_sp3_neg ( Element.O , Hybridization.SP3    , -1, 1.80f, 0.4356f),
		O_sp3     ( Element.O , Hybridization.SP3    ,  0, 3.20f, 0.4356f),
		O_ca      ( Element.O , Hybridization.CA     ,  0, 3.44f, 0.3600f),
		F_sp3     ( Element.F , Hybridization.SP3    ,  0, 4.00f, 0.4096f),
		F_ion_n1  ( Element.F , Hybridization.ION    , -1, 3.98f, 1.7689f),
		F_ion_p7  ( Element.F , Hybridization.ION    , +7, 3.98f, 0.0064f),
		F_ionca   ( Element.F , Hybridization.IONCA  ,  0, 4.00f, 0.4096f),
		F_ca      ( Element.F , Hybridization.CA     ,  0, 3.98f, 0.2500f),
		//Ne_ca     ( Element.Ne, Hybridization.CA     ,  0, ?.??f, 0.3364f),  // Helium is not known to form covalent bonds. Ions such as NeH+, HeNe+ might exist.
		Na_ion_p1 ( Element.Na, Hybridization.ION    , +1, 0.93f, 1.0404f),
		Na_ionca  ( Element.Na, Hybridization.IONCA  ,  0, 2.20f, 2.0000f),
		Na_ca     ( Element.Na, Hybridization.CA     ,  0, 0.93f, 3.2400f),
		Mg_sp3_pos( Element.Mg, Hybridization.SP3    , +2, 2.20f, 2.0000f),
		Mg_ion_p2 ( Element.Mg, Hybridization.ION    , +2, 1.31f, 0.5184f),
		Mg_ionca  ( Element.Mg, Hybridization.IONCA  ,  0, 1.31f, 0.5184f),
		Mg_ca     ( Element.Mg, Hybridization.CA     ,  0, 1.31f, 2.2500f),
		Al_ion_p3 ( Element.Al, Hybridization.ION    ,  0, 1.61f, 0.2852f),
		Al_ionca  ( Element.Al, Hybridization.IONCA  ,  0, 1.61f, 1.5625f),
		Al_ca     ( Element.Al, Hybridization.CA     ,  0, 1.61f, 1.5625f),
		Si_sp3    ( Element.Si, Hybridization.SP3    ,  0, 1.99f, 1.2321f),
		Si_sp3d   ( Element.Si, Hybridization.SP3D   ,  0, 1.99f, 1.2321f),
		Si_ion_p4 ( Element.Si, Hybridization.ION    , +4, 1.90f, 0.1600f),
		Si_ca     ( Element.Si, Hybridization.CA     ,  0, 1.90f, 1.2100f),
		P_sp3     ( Element.P , Hybridization.SP3    ,  0, 2.20f, 1.2100f),
		P_sp3_pos ( Element.P , Hybridization.SP3    , +1, 2.20f, 1.2100f),
		P_dsp3    ( Element.P , Hybridization.SP3D   ,  0, 2.26f, 1.2100f),
		P_sp3d2   ( Element.P , Hybridization.SP3D2  , -1, 2.20f, 1.2100f), // Copied values.
		P_ion_p3  ( Element.P , Hybridization.ION    , +3, 2.19f, 0.1936f),
		P_ion_p5  ( Element.P , Hybridization.ION    , +5, 2.19f, 0.1444f),
		P_ca      ( Element.P , Hybridization.CA     ,  0, 2.19f, 1.0000f),
		S_sp2_neg ( Element.S , Hybridization.SP2    , -1, 2.80f, 0.8836f),
		S_sp2     ( Element.S , Hybridization.SP2    ,  0, 2.80f, 0.8836f),
		S_sp2_pos ( Element.S , Hybridization.SP2    , +1, 2.80f, 0.8836f),
		S_sp3_neg ( Element.S , Hybridization.SP3    , -1, 2.80f, 0.8836f),
		S_sp3     ( Element.S , Hybridization.SP3    ,  0, 2.74f, 1.0816f),
		S_sp3_pos ( Element.S , Hybridization.SP3    , +1, 4.79f, 0.8836f),
		S_sp3_ppos( Element.S , Hybridization.SP3    , +2, 4.80f, 0.8836f),
		S_sp3d2   ( Element.S , Hybridization.SP3D2  ,  0, 2.96f, 0.8836f),  // Copied radius, set eo to Sanderson value. ?
		S_ion_n2  ( Element.S , Hybridization.ION    , -2, 2.58f, 3.3856f),
		S_ion_p4  ( Element.S , Hybridization.ION    , +4, 2.58f, 0.1369f),
		S_ion_p6  ( Element.S , Hybridization.ION    , +6, 2.58f, 0.0841f),
		S_ca      ( Element.S , Hybridization.CA     ,  0, 2.58f, 1.0000f),
		Cl_sp3    ( Element.Cl, Hybridization.SP3    ,  0, 3.28f, 0.9801f),
		Cl_ion_n1 ( Element.Cl, Hybridization.ION    , -1, 3.16f, 3.2761f),
		Cl_ion_p7 ( Element.Cl, Hybridization.ION    , +7, 3.16f, 0.0729f),
		Cl_ca     ( Element.Cl, Hybridization.CA     ,  0, 3.16f, 1.0000f),
		//Ar_ca     ( Element.Ar, Hybridization.CA     ,  0, ?.??f, ?.????f),  // Stable compounds (HArF, XePtF6) are known, but no data available.
		K_ion_p1  ( Element.K , Hybridization.ION    , +1, 0.82f, 1.9044f),
		K_ionca   ( Element.K , Hybridization.IONCA  ,  0, 2.20f, 2.0000f),
		K_ca      ( Element.K , Hybridization.CA     ,  0, 0.82f, 4.8400f),
		Ca_sp3_p1 ( Element.Ca, Hybridization.SP3    , +2, 2.20f, 2.0000f),
		Ca_ion_p2 ( Element.Ca, Hybridization.ION    , +2, 1.00f, 1.0000f),
		Ca_ca     ( Element.Ca, Hybridization.CA     ,  0, 1.00f, 3.2400f),
		Sc_ca     ( Element.Sc, Hybridization.CA     ,  0, 1.36f, 2.0736f),
		Ti_ca     ( Element.Ti, Hybridization.CA     ,  0, 1.54f, 1.8496f),
		V_ca      ( Element.V , Hybridization.CA     ,  0, 1.63f, 1.5625f),
		Cr_ca     ( Element.Cr, Hybridization.CA     ,  0, 1.66f, 1.6129f),
		Mn_sp3_p2 ( Element.Mn, Hybridization.SP3    , +2, 2.20f, 2.0000f),
		Mn_ca     ( Element.Mn, Hybridization.CA     ,  0, 1.55f, 1.9321f),
		Fe_sp3_p2 ( Element.Fe, Hybridization.SP3    , +2, 2.20f, 2.0000f),
		Fe_ion_p4 ( Element.Fe, Hybridization.ION    , +4, 1.83f, 0.3423f),
		Fe_ion_p6 ( Element.Fe, Hybridization.ION    , +6, 1.83f, 0.0625f),
		Fe_ca     ( Element.Fe, Hybridization.CA     , +2, 1.83f, 1.5625f),
		Co_sp3    ( Element.Co, Hybridization.SP3    ,  0, 2.20f, 2.0000f),
		Co_ca     ( Element.Co, Hybridization.CA     ,  0, 1.88f, 1.5876f),
		Ni_ca     ( Element.Ni, Hybridization.CA     ,  0, 1.91f, 1.4641f),
		Cu_sp3_p1 ( Element.Cu, Hybridization.SP3    , +1, 2.20f, 2.0000f),
		Cu_ca     ( Element.Cu, Hybridization.CA     ,  0, 1.90f, 1.9044f),
		Zn_sp3_p2 ( Element.Zn, Hybridization.SP3    , +2, 2.20f, 2.0000f),
		Zn_ca     ( Element.Zn, Hybridization.CA     ,  0, 1.65f, 1.7161f),
		Ga_ion_p3 ( Element.Ga, Hybridization.ION    , +3, 1.81f, 0.3844f),
		Ga_ca     ( Element.Ga, Hybridization.CA     ,  0, 1.81f, 1.6900f),
		Ge_ion_p2 ( Element.Ge, Hybridization.ION    , +2, 2.01f, 0.5329f),
		Ge_ion_p4 ( Element.Ge, Hybridization.ION    , +4, 2.01f, 0.2809f),
		Ge_ca     ( Element.Ge, Hybridization.CA     ,  0, 2.01f, 1.5625f), 
		As_sp2    ( Element.As, Hybridization.SP2    ,  0, 2.38f, 1.4641f),
		As_sp3    ( Element.As, Hybridization.SP3    ,  0, 2.38f, 1.4641f),
		As_sp3_p1 ( Element.As, Hybridization.SP3    , +1, 2.69f, 1.4641f),
		As_ca     ( Element.As, Hybridization.CA     ,  0, 2.18f, 1.3225f),
		Se_sp2    ( Element.Se, Hybridization.SP2    ,  0, 2.54f, 1.3689f),  // Copied values from SP3-hybridized Se.
		Se_sp3    ( Element.Se, Hybridization.SP3    ,  0, 2.54f, 1.3689f),
		Se_ion_n2 ( Element.Se, Hybridization.ION    , -2, 2.55f, 3.2904f),
		Se_ion_p4 ( Element.Se, Hybridization.ION    , +4, 2.55f, 0.2500f),
		Se_ion_p6 ( Element.Se, Hybridization.ION    , +6, 2.55f, 0.1764f),
		Se_ca     ( Element.Se, Hybridization.CA     ,  0, 2.55f, 1.3225f),
		Br_sp3    ( Element.Br, Hybridization.SP3    ,  0, 3.13f, 1.2996f),
		Br_ion_n1 ( Element.Br, Hybridization.ION    , -1, 2.96f, 3.8416f),
		Br_ion_p3 ( Element.Br, Hybridization.ION    , +3, 2.96f, 0.3481f),
		Br_ion_p5 ( Element.Br, Hybridization.ION    , +5, 2.96f, 0.0961f),
		Br_ion_p7 ( Element.Br, Hybridization.ION    , +7, 2.96f, 0.1521f),
		Br_ca     ( Element.Br, Hybridization.CA     ,  0, 2.96f, 1.3225f),
		Kr_ca     ( Element.Kr, Hybridization.CA     ,  0, 3.00f, 1.2544f),  // Covalent radius from webelements.com.
		Rb_ion_p1 ( Element.Rb, Hybridization.ION    , +1, 0.82f, 2.3104f),
		Rb_ca     ( Element.Rb, Hybridization.CA     ,  0, 0.82f, 5.5225f),
		Sr_ion_p2 ( Element.Sr, Hybridization.ION    , +2, 0.95f, 1.3924f),
		Sr_ca     ( Element.Sr, Hybridization.CA     ,  0, 0.95f, 4.0000f),
		Y_ca      ( Element.Y , Hybridization.CA     ,  0, 1.22f, 2.6244f),
		Zr_ca     ( Element.Zr, Hybridization.CA     ,  0, 1.33f, 2.1904f),
		Nb_ca     ( Element.Nb, Hybridization.CA     ,  0, 1.60f, 1.8769f),
		Mo_ca     ( Element.Mo, Hybridization.CA     ,  0, 2.16f, 2.1025f),
		Tc_ca     ( Element.Tc, Hybridization.CA     , +3, 1.90f, 2.4336f),
		Ru_ca     ( Element.Ru, Hybridization.CA     ,  0, 2.20f, 1.5876f),
		Rh_ca     ( Element.Rh, Hybridization.CA     ,  0, 2.28f, 1.8225f),
		Pd_ca     ( Element.Pd, Hybridization.CA     ,  0, 2.20f, 1.7161f),
		Ag_ca     ( Element.Ag, Hybridization.CA     ,  0, 1.93f, 2.3409f),
		Cd_ca     ( Element.Cd, Hybridization.CA     ,  0, 1.69f, 2.1904f),
		In_ion_p3 ( Element.In, Hybridization.ION    , +3, 1.78f, 0.6400f),
		In_ca     ( Element.In, Hybridization.CA     ,  0, 1.78f, 2.4025f),
		Sn_sp3    ( Element.Sn, Hybridization.SP3    ,  0, 1.96f, 1.4161f),
		Sn_sp3_n1 ( Element.Sn, Hybridization.SP3    , -1, 1.96f, 1.4161f), // Copied values from Sn_sp3.
		Sn_sp3_n2 ( Element.Sn, Hybridization.SP3    , -2, 1.96f, 1.4161f),
		Sn_ion_n2 ( Element.Sn, Hybridization.ION    , -2, 1.96f, 3.9204f),
		Sn_ion_p4 ( Element.Sn, Hybridization.ION    , +4, 1.96f, 0.2500f),
		Sn_ion_p6 ( Element.Sn, Hybridization.ION    , +6, 1.96f, 0.1764f),
		Sn_ca     ( Element.Sn, Hybridization.CA     ,  0, 1.96f, 2.1025f),
		Sb_sp3_n1 ( Element.Sb, Hybridization.SP3    , -1, 2.05f, 1.9044f),
		Sb_sp3    ( Element.Sb, Hybridization.SP3    ,  0, 1.41f, 1.9881f),
		Sb_sp3_p1 ( Element.Sb, Hybridization.SP3    , +1, 1.41f, 1.9881f),
		Sb_ion_p3 ( Element.Sb, Hybridization.ION    , +3, 2.05f, 0.5776f),
		Sb_ion_p5 ( Element.Sb, Hybridization.ION    , +5, 2.05f, 0.3600f),
		Sb_ca     ( Element.Sb, Hybridization.CA     ,  0, 2.05f, 2.1025f),
		Te_ion_n2 ( Element.Te, Hybridization.ION    , -2, 2.10f, 4.8841f),
		Te_ion_p4 ( Element.Te, Hybridization.ION    , +4, 2.10f, 0.9409f),
		Te_ion_p6 ( Element.Te, Hybridization.ION    , +6, 2.10f, 0.3136f),
		Te_ca     ( Element.Te, Hybridization.CA     ,  0, 2.10f, 1.9600f),	    
		I_sp3     ( Element.I , Hybridization.SP3    ,  0, 2.93f, 1.7689f),
		I_sp3_p   ( Element.I , Hybridization.SP3    , +1, 2.93f, 1.7689f),  // Copied values from I_sp3.
		I_ion_n1  ( Element.I , Hybridization.ION    , -1, 2.66f, 4.8400f),
		I_ion_p5  ( Element.I , Hybridization.ION    , +5, 2.66f, 0.9025f),
		I_ion_p7  ( Element.I , Hybridization.ION    , +7, 2.66f, 0.2809f),
		I_ca      ( Element.I , Hybridization.CA     ,  0, 2.66f, 1.9600f),
		Xe_ionca  ( Element.Xe, Hybridization.IONCA  ,  0, 2.60f, 0.2304f),
		Xe_ca     ( Element.Xe, Hybridization.CA     ,  0, 2.60f, 1.9600f),
		Cs_ion_p1 ( Element.Cs, Hybridization.ION    , +1, 0.79f, 2.7889f),
		Cs_ca     ( Element.Cs, Hybridization.CA     ,  0, 0.79f, 6.7600f),
		Ba_ca     ( Element.Ba, Hybridization.CA     ,  0, 0.89f, 4.6225f),  // Values from wikipedia and webelements.com
		La_ca     ( Element.La, Hybridization.CA     ,  0, 1.10f, 3.8025f),
		Ce_ion_p3 ( Element.Ce, Hybridization.ION    ,  0, 1.12f, 1.0404f),
		Ce_ion_p4 ( Element.Ce, Hybridization.ION    ,  0, 1.12f, 0.7569f),
		Ce_ca     ( Element.Ce, Hybridization.CA     ,  0, 1.12f, 3.4225f),
		Pr_ca     ( Element.Pr, Hybridization.CA     ,  0, 1.13f, 3.4225f),
		Nd_ca     ( Element.Nd, Hybridization.CA     ,  0, 1.14f, 3.4225f),
		Pm_ca     ( Element.Pm, Hybridization.CA     ,  0, 1.13f, 3.4225f),
		Sm_ca     ( Element.Sm, Hybridization.CA     ,  0, 1.17f, 3.4225f),
		Eu_ca     ( Element.Eu, Hybridization.CA     ,  0, 1.20f, 3.4225f),
		Gd_ca     ( Element.Gd, Hybridization.CA     ,  0, 1.20f, 3.2400f),
		Tb_ca     ( Element.Tb, Hybridization.CA     ,  0, 1.10f, 3.0625f),
		Dy_ca     ( Element.Dy, Hybridization.CA     ,  0, 1.22f, 3.0625f),
		Ho_ca     ( Element.Ho, Hybridization.CA     ,  0, 1.23f, 3.0625f),
		Er_ca     ( Element.Er, Hybridization.CA     ,  0, 1.24f, 3.0625f),
		Tm_ca     ( Element.Tm, Hybridization.CA     ,  0, 1.25f, 3.0625f),
		Yb_ca     ( Element.Yb, Hybridization.CA     ,  0, 1.10f, 3.0625f),
		Lu_ca     ( Element.Lu, Hybridization.CA     ,  0, 1.27f, 3.0625f),
		Hf_ca     ( Element.Hf, Hybridization.CA     ,  0, 1.30f, 2.2500f),
		Ta_ca     ( Element.Ta, Hybridization.CA     ,  0, 1.50f, 1.9044f),
		W_ca      ( Element.W , Hybridization.CA     ,  0, 2.36f, 2.1316f),
		Re_ca     ( Element.Re, Hybridization.CA     , +1, 1.90f, 2.5281f),
		Os_ca     ( Element.Os, Hybridization.CA     ,  0, 2.20f, 1.6384f),
		Ir_ca     ( Element.Ir, Hybridization.CA     ,  0, 2.20f, 1.8769f),
		Pt_ca     ( Element.Pt, Hybridization.CA     ,  0, 2.28f, 1.6384f),
		Au_ca     ( Element.Au, Hybridization.CA     ,  0, 2.54f, 2.0736f),
		Hg_ca     ( Element.Hg, Hybridization.CA     ,  0, 2.00f, 2.2201f),
		Tl_ca     ( Element.Tl, Hybridization.CA     ,  0, 1.62f, 2.2201f),	    
		Pb_ion_p2 ( Element.Pb, Hybridization.ION    , +2, 2.33f, 1.4161f),
		Pb_ion_p4 ( Element.Pb, Hybridization.ION    , +4, 2.33f, 0.6006f),
		Pb_ca     ( Element.Pb, Hybridization.CA     ,  0, 2.33f, 3.2400f),
		Bi_ca     ( Element.Bi, Hybridization.CA     ,  0, 2.02f, 2.5600f),
		Po_ion_p4 ( Element.Po, Hybridization.ION    , +4, 2.00f, 0.8836f),
		Po_ion_p6 ( Element.Po, Hybridization.ION    , +6, 2.00f, 0.4489f),
		Po_ca     ( Element.Po, Hybridization.CA     ,  0, 2.00f, 3.6100f),
		At_ion_p7 ( Element.At, Hybridization.ION    , +7, 2.20f, 0.3844f),
		At_ca     ( Element.At, Hybridization.CA     ,  0, 2.20f, 2.1025f),  // Covalent radius from webelements.com.
		//Rn_ca     ( Element.Rn, Hybridization.CA     ,  0, 2.20f, ?.????f),  // Information about covalent radius not available.
		Fr_ion_p1 ( Element.Fr, Hybridization.ION    , +1, 0.70f, 3.2400f),
		//Fr_ca     ( Element.Fr, Hybridization.CA     ,  0, 0.70f, ?.????f),  // Information about covalent radius not available.
		Ra_ion_p2 ( Element.Ra, Hybridization.ION    , +2, 0.90f, 2.1904f),
		Ra_ca     ( Element.Ra, Hybridization.CA     ,  0, 0.90f, 4.6225f),
		Ac_ca     ( Element.Ac, Hybridization.CA     ,  0, 1.10f, 3.8025f),
		Th_ca     ( Element.Th, Hybridization.CA     ,  0, 1.30f, 3.2400f),
		Pa_ca     ( Element.Pa, Hybridization.CA     ,  0, 1.50f, 3.2400f),
		U_ca      ( Element.U , Hybridization.CA     ,  0, 1.38f, 3.0625f),
		Np_ca     ( Element.Np, Hybridization.CA     ,  0, 1.36f, 3.0625f),
		Pu_ca     ( Element.Pu, Hybridization.CA     ,  0, 1.28f, 3.0625f),
		Am_ca     ( Element.Am, Hybridization.CA     ,  0, 1.13f, 3.0625f)
		;

		/** 
		 * Hybridization states for atom typing.
		 * 
		 * IONCA = ion catch all, CA = catch all (fall-back types).
		 */
		static enum Hybridization { NONE, SP, SP2, SP3, SP3D, SP3D2, UNKNOWN, ION, IONCA, CA };  // IONCA = Ion Catch All, CA = Catch All.

		/** Returns the 'inductive' effective electronegativity of an atom type. */
		float inductiveElectronegativity() { return electronegativity; }

		/** Returns the squared covalent radius of an atom type. */
		float squaredCovalentRadius() { return squaredCovalentRadius; }

		/** Returns the appropriate atom type for given element, hybridization, and ionization state. */
		static AtomType atomType(final Element element, final Hybridization hybridization, final int ionization) throws ComputationFailureException 
		{
			for (final AtomType type : AtomType.values()) 
			{
				if(type.element == element && type.hybridization == hybridization && type.ionization == ionization) return type;  // Exact match.
				if(type.element == element && hybridization == Hybridization.ION && type.hybridization == Hybridization.IONCA) return type;  // Fallback: generic ion type.
				if(type.element == element && type.hybridization == Hybridization.CA) return type;  // Fallback: generic element type.
			}
			throw new ComputationFailureException(String.format("Unknown atom type encountered (%s, %s, %d)", element.abbrev(), hybridization.toString(), ionization));
		}

		// Internal representation.
		private final Element element;
		private final Hybridization hybridization;
		private final int ionization;
		private final float electronegativity;  // 'Inductive" effective electronegativity.
		private final float squaredCovalentRadius;  // Squared covalent radius in Angstrom^2.

		AtomType(final Element element, final Hybridization hybridization, final int ionization, final float electronegativity, final float squaredCovalentRadius) 
		{ 
			this.element               = element;
			this.hybridization         = hybridization;
			this.ionization            = ionization;
			this.electronegativity     = electronegativity;
			this.squaredCovalentRadius = squaredCovalentRadius;
		}
	}
}
