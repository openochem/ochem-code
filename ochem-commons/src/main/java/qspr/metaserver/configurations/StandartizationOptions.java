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

package qspr.metaserver.configurations;

import java.io.Serializable;
import java.util.Arrays;
import java.util.List;

import javax.xml.bind.annotation.XmlRootElement;

import qspr.dao.ChemInfEngine;
import qspr.metaserver.util.MoleculeStandartizer.MolFormat;

@XmlRootElement(name = "standartization")
public class StandartizationOptions implements Serializable
{
	private static final long serialVersionUID = 1L;

	public ChemInfEngine desaltWith;
	public ChemInfEngine neutralizeWith;
	public ChemInfEngine standardizeWith;
	public ChemInfEngine cleanStructureWith;
	public ChemInfEngine addExplicitHydrogensWith;
	public ChemInfEngine dearomatizeWith;

	public MolFormat outFormat;

	public StandartizationOptions() {
	}

	public StandartizationOptions(boolean defaultOptions) {
		if(defaultOptions)
			setDefaults();
	}

	public boolean standartizationRequired()
	{
		return desaltWith != null || neutralizeWith != null || standardizeWith != null || cleanStructureWith != null || addExplicitHydrogensWith != null;
	}

	public StandartizationOptions setDefaults() {
		desaltWith = neutralizeWith = standardizeWith = cleanStructureWith = ChemInfEngine.CDK;
		return this;
	}

	public StandartizationOptions setDefaults(ChemInfEngine standardizer) {
		desaltWith = neutralizeWith = standardizeWith = cleanStructureWith = standardizer;
		return this;
	}

	public static List<ChemInfEngine> getStandardizers() {
		return Arrays.asList(ChemInfEngine.values());
	}

	public static void main(String args[]) {

		ChemInfEngine ret;
		ret = ChemInfEngine.valueOf("CDK"); 
		System.out.println("Selected : " + ret);                              
	}

	public ChemInfEngine getDefault() {
		ChemInfEngine one = 
				desaltWith != null ? desaltWith: 
					neutralizeWith != null ? neutralizeWith:
						standardizeWith != null ? standardizeWith:
							cleanStructureWith != null? cleanStructureWith : 
								addExplicitHydrogensWith != null? addExplicitHydrogensWith : dearomatizeWith;
		if (one == null)return ChemInfEngine.CDK; // nothing is defined, CDK is used
		return one;
	}

	public ChemInfEngine synchronise() {
		ChemInfEngine touse = getDefault();
		ChemInfEngine second = desaltWith != null && desaltWith != touse ? desaltWith: 
			neutralizeWith != null && neutralizeWith != touse ? neutralizeWith:
				standardizeWith != null && standardizeWith != touse ? standardizeWith:
					cleanStructureWith != null && cleanStructureWith != touse? cleanStructureWith : 
						addExplicitHydrogensWith != null && addExplicitHydrogensWith != touse? addExplicitHydrogensWith : 
							dearomatizeWith != null? dearomatizeWith: touse;

		if(second != touse) touse = ChemInfEngine.CDK; // conflict, using default, which is CHEMAXON
		
		if (desaltWith != null) desaltWith = touse;
		if (neutralizeWith != null) neutralizeWith = touse;
		if (standardizeWith != null) standardizeWith = touse;
		if (cleanStructureWith != null) cleanStructureWith = touse;
		if (addExplicitHydrogensWith != null) addExplicitHydrogensWith = touse;
		if (dearomatizeWith != null) dearomatizeWith = touse;

		return touse;
	}

}
