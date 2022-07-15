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

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.eadmet.exceptions.UserFriendlyException;

import qspr.dao.Various;
import qspr.metaserver.configurations.DescriptorsConfiguration.MixturesProcessing;
import qspr.metaserver.util.MixtureAttachment;
import qspr.workflow.utils.QSPRConstants;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;

abstract public class MixtureDescriptorsProcessor 
{
	public static MixtureDescriptorsProcessor getInstance(MixturesProcessing processing) throws IOException
	{
		switch (processing)
		{
		case SUMONLY:
			return new WAverageOnlyProcessor();
		case FRACTION:
			return new FractionMixtureProcessor();
		case CONCATENATE:
			return new ConcatenateMixtureProcessor();
		case WSUM:
			return new WSumMixtureProcessor();
		case WSUMDIFF:
			return new WSumDiffMixtureProcessor();
		default:
			throw new IOException("Mixture processor is not defined");
		}
	}

	public DataTable processMixtureDescriptors(DataTable dtMolecules, DataTable dtDescriptors)
	{
		Map<String, Integer> componentsToRows = new HashMap<String, Integer>();

		for(int i=0;i<dtDescriptors.getRowsSize();i++)
			try { // looking only for single molecules
				String split[] = Various.molecule.splitOrderedByChargeAndSize(dtMolecules.getSDF(i));
				if(split.length>1)continue;
				String inchi =(String)dtMolecules.getRow(i).getAttachment(QSPRConstants.INCHIKEYS);
				if(inchi == null)
					throw new UserFriendlyException("Inchi should not be null");
				componentsToRows.put(inchi, i);
			}catch(Exception e) {
			}

		DataTable newDescriptors=dtDescriptors.getEmptyCopy();
		newDescriptors.columns = new  ArrayList<String>();

		for(int i=0;i<dtDescriptors.getRowsSize();i++)try{
			newDescriptors.addRow();

			MixtureAttachment ma = (MixtureAttachment)dtMolecules.getRow(i).getAttachment(QSPRConstants.MIXTURE_ATTACHMENT);

			if(ma == null) ma = new MixtureAttachment((String)dtMolecules.getValue(i, 0)); // create on the fly :)
			else
				ma.validateMixtureAttachment();

			AbstractDataRow rows[];

			if(this instanceof FractionMixtureProcessor) {
				rows = new AbstractDataRow[1];
				rows[0] =  dtDescriptors.getRow(i);  // whole molecule descriptors are used
			}
			else {
				rows =  new AbstractDataRow[ma.fractions.size()]; // fractions
				int j = 0;
				for(String key:ma.fractions.keySet()){
					if(!componentsToRows.containsKey(key))throw new IOException("No descriptors were found for component with inchies: " + key);
					int row = componentsToRows.get(key);
					dtDescriptors.filterNaN(row);
					rows[j] = dtDescriptors.getRow(row);
					if(rows[j].isError()) throw new IOException(rows[j].detailedStatus + " for " + Various.molecule.convertToCanonicalSMILES(dtMolecules.getSDF(row)));
					j++;
				}
			}

			processMixtureCompound(newDescriptors, dtDescriptors.getColumns(),rows,normaliseFraction(ma.getMolarFractions()));
			newDescriptors.getCurrentRow().addAttachment(QSPRConstants.MIXTURE_ATTACHMENT, ma); // attachment will be further required for correct validation
		}catch(Exception e) {
			newDescriptors.getCurrentRow().setError(e.getMessage());
		}

		return newDescriptors;
	}

	/**
	 * Processing and creation of descriptors for a single row
	 * @param newDes
	 * @param list
	 * @param rows
	 * @param molarFactions
	 */


	abstract void processMixtureCompound(DataTable newDes, List<String> list, AbstractDataRow[] rows, double[] molarFactions);

	double [] normaliseFraction(double [] val) throws IOException {
		double sum =0;
		for(double v:val) {
			if(v < 0)throw new IOException("one of fraction was less than 0 < " + v);
			sum += v;
		}

		if(sum == 0)
			for(int i=0;i<val.length;i++)
				val[i] = 1./val.length;
		else
			if(sum>1.001) throw new IOException("sum of fraction >1: " + sum);
			else
				if(sum < 0.999)  throw new IOException("sum of fraction <1: " + sum);
				else
					val[0] = val[0] + (1-sum);
		return val;
	}

	void addFractions(DataTable newDes, double []fractions, List<String> descriptors, AbstractDataRow[] rows) {
		for (int i=0; i<fractions.length; i++)
			newDes.setValue("fraction_"+i,fractions[i]);

		if(descriptors != null)
			for (int i=0; i<descriptors.size(); i++)
				newDes.setValue(descriptors.get(i), rows[0].getValue(i)); // only one row is assumed and added
	}
}
class FractionMixtureProcessor extends MixtureDescriptorsProcessor
{
	@Override
	void processMixtureCompound(DataTable newDes, List<String> descriptors, AbstractDataRow[] rows, double []fractions) {
		addFractions(newDes, fractions,  descriptors, rows);
	}
}

class ConcatenateMixtureProcessor extends MixtureDescriptorsProcessor
{
	@Override
	void processMixtureCompound(DataTable newDes, List<String> descriptors, AbstractDataRow[] rows, double []fractions) {
		if(rows.length != 2) {
			newDes.getCurrentRow().setError("For concatenate descriptors only binary mixtures are supported");
			return;
		}

		addFractions(newDes, fractions, null, null);
		for (int i=0; i<descriptors.size(); i++){
			for (int j=0;j<rows.length;j++)
				newDes.setValue(descriptors.get(i) +"_" +j, rows[j].getValue(i));
		}
	}
}

class WAverageOnlyProcessor extends MixtureDescriptorsProcessor
{
	String suffix="";

	@Override
	void processMixtureCompound(DataTable newDes, List<String> descriptors, AbstractDataRow[] rows, double []fractions) {
		for (int i=0; i<descriptors.size(); i++){
			double val = 0.; int j =0;
			for (AbstractDataRow row:rows) 
				val += (Double)row.getValue(i)*fractions[j++];
			newDes.setValue(descriptors.get(i) +suffix, val);
		}
	}
}

class WSumMixtureProcessor extends MixtureDescriptorsProcessor
{
	String suffix="_sum";

	@Override
	void processMixtureCompound(DataTable newDes, List<String> descriptors, AbstractDataRow[] rows, double []fractions) {
		addFractions(newDes, fractions, null, null);
		for (int i=0; i<descriptors.size(); i++){
			double val = 0.; int j =0;
			for (AbstractDataRow row:rows) 
				val += (Double)row.getValue(i)*fractions[j++];
			newDes.setValue(descriptors.get(i) +suffix, val);
		}
	}
}

class WSumDiffMixtureProcessor extends MixtureDescriptorsProcessor
{
	String suffix="_sum", suffix1="_diff";

	@Override
	void processMixtureCompound(DataTable newDes, List<String> descriptors, AbstractDataRow[] rows, double []fractions) {
		addFractions(newDes, fractions, null, null);
		if(fractions.length == 2)
			for (int i=0; i<descriptors.size(); i++)
			{
				double val = Math.abs((Double)rows[0].getValue(i)*fractions[0] + (Double)rows[1].getValue(i))*fractions[1];
				newDes.setValue(descriptors.get(i) +suffix, val);
				val = Math.abs((Double)rows[0].getValue(i)*fractions[0] - (Double)rows[1].getValue(i))*fractions[1];
				newDes.setValue(descriptors.get(i) +suffix1, val);
			}
		else
			newDes.getCurrentRow().setError("WSumDiffMixtureProcessor is only defined for bi-componental mixtures but this is a mixture of " + fractions.length + " components");
	}
}

