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

package qspr.metaserver.util;

import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeMap;

import com.eadmet.utils.NumericalValueStandardizer;

import qspr.dao.Various;

// Mixture Related Code
public class MixtureAttachment implements Serializable 
{
	private static final long serialVersionUID = 1L;

	public TreeMap<String,Double> fractions = new TreeMap<String,Double>();
	public List<String> smiles = new ArrayList<String>(); // may store SMILES or Molecular IDS (for faster creation of dataset for solvents) of molecules

	public MixtureAttachment(){
	}

	public MixtureAttachment(String sdf) throws IOException {

		String mols[]= Various.molecule.splitOrderedByChargeAndSize(sdf);
		for(String s:mols) {
			String inchi = Various.molecule.getInChIKeyNoStero(s);
			double val = 0;
			if(fractions.containsKey(inchi)) val = fractions.get(inchi);
			fractions.put(inchi, val + 1./mols.length);
		}
	}

	public void validateMixtureAttachment() throws IOException {
		double sum = 0;
		for(Entry<String,Double> e:fractions.entrySet())
			sum += e.getValue();
		if(Math.abs(sum - 1.)>0.0001)throw new IOException("sum of molar fractions is not 1: sum=" + sum);
	}

	public double[] getMolarFractions() {
		double val[] = new double[fractions.size()];
		int i = 0;
		for(Entry<String,Double> e:fractions.entrySet())
			val[i++] = e.getValue();
		return val; 
	}
	
	public String toString() {
		String str="";
		int i =0;
		for(Entry<String,Double> e:fractions.entrySet())
			str += (str.length() == 0? "":" ") + smiles.get(i++) + ";" + NumericalValueStandardizer.getSignificantDigits(e.getValue());
		return str;	
	}

	public String smiles() {
		String str="";
		for(String smile:smiles)
			str+= (str.length() == 0? "":".") + smile;
		return str;	
	}
}
