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

package qspr.metaserver.cs.util;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringWriter;
import java.io.Writer;

//import javafx.util.Pair;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.utils.NumericalValueStandardizer;
import com.eadmet.utils.OCHEMUtils;

public class LIBSVMUtils {

	public static void writeLibSvmTrainingSet(DescriptorsTable dtDescriptors, LabelsTable dtExpValues,  String filename) throws IOException{
		writeLibSvmTrainingSet(dtDescriptors,dtExpValues,null,filename);
	}

	public static void writeLibSvmTrainingSet(DescriptorsTable dtDescriptors, LabelsTable dtExpValues, String singleClassLabel, String filename) throws IOException
	{
		File f = new File(filename);
		if(f.exists())return;

		FileWriter fis = new FileWriter(f);
		BufferedWriter bw = new BufferedWriter(fis,QSPRConstants.FILE_BUFFER);

		for (int mol = 0; mol < dtExpValues.getDataSize(); mol++)
		{
			String label = dtExpValues.getOneValueOnlyString(mol); 

			if (singleClassLabel != null && !singleClassLabel.equals(label)) //Training set for one-class classification holds only one selected class
				continue;

			bw.write(label);
			writeLibSvmDescriptorLine(dtDescriptors, mol, bw);
			bw.write("\n");
		}
		bw.close();

	}

	/*	
	public static void writeMultiClassValues(DescriptorsTable dtDescriptors, LabelsTable dtExpValues, String filename) throws IOException {
		FileWriter fis = new FileWriter(filename);
		BufferedWriter bw = new BufferedWriter(fis,QSPRConstants.FILE_BUFFER);

		for (int mol = 0; mol < dtExpValues.getDataSize(); mol++){

			String values[] = dtExpValues.getScaledValuesString(mol);
			for(int n = 0; n < values.length; n++){
				if(!values[n].equals(LabelsTable.MISSED_VALUES))
					bw.append(values[n]);
				if(n != values.length -1) bw.append(",");
			}
			writeLibSvmDescriptorLine(dtDescriptors, mol, bw);
			bw.write("\n");
		}
		bw.close();
	}
	 */ 

	static String getHash(DescriptorsTable dtDescriptors, int mol) throws IOException {
		if(dtDescriptors.getRawRow(mol).attachments.containsKey(QSPRConstants.EXTERNAL_ID))
			return (String)dtDescriptors.getRawRow(mol).attachments.get(QSPRConstants.EXTERNAL_ID); 
		StringWriter sw = new StringWriter();
		sw.getBuffer().setLength(0);
		writeLibSvmDescriptorLine(dtDescriptors, mol, sw);
		return OCHEMUtils.getMD5(sw.toString());
	}

	/**
	 *  Aggregates values for molecules with similar descriptors
	 *  hasVal may be not unique for large datasets
	 * @param dtDescriptors
	 * @param dtExpValues
	 * @param filename
	 * @param missedValues 
	 * @throws IOException
	 */

	/*
	private static void writeMultiClassValues(DescriptorsTable dtDescriptors, LabelsTable dtExpValues, String filename) throws IOException {

		Map<String,ArrayList<Map.Entry<Integer,String>>> val = new LinkedHashMap<String,ArrayList<Map.Entry<Integer,String>>> ();
		ArrayList<Map.Entry<Integer,String>> vals ;
		int outputs = dtExpValues.getScaledValuesString(0).length;
		String values[], hashVal, value;

		// collection of data according to hash values by descriptors
		for (int mol = 0; mol < dtExpValues.getDataSize(); mol++){

			hashVal = getHash(dtDescriptors,mol);

			if(val.containsKey(hashVal))
				vals = val.get(hashVal);
			else
			{
				vals = new ArrayList<Map.Entry<Integer,String>>();
				val.put(hashVal,vals);
			}

			values = dtExpValues.getScaledValuesString(mol);
			for(int n = 0; n < outputs; n++)
				if(!values[n].equals(LabelsTable.MISSED_VALUES))
					vals.add(new AbstractMap.SimpleEntry<Integer,String>(n,values[n])); // all non missed values are stored for each descriptor combination, e.g. molecule
		}

		System.out.println("Unique " + val.size() + " for total "+dtExpValues.getDataSize());

		// saving all values
		FileWriter fis = new FileWriter(filename);
		BufferedWriter bw = new BufferedWriter(fis,QSPRConstants.FILE_BUFFER);

		FileWriter fi = new FileWriter(filename+".ids");
		BufferedWriter bi= new BufferedWriter(fi,QSPRConstants.FILE_BUFFER);

		StringWriter sw = new StringWriter();

		for (int mol = 0; mol < dtExpValues.getDataSize(); mol++){
			hashVal = getHash(dtDescriptors,mol); 
			if(!val.containsKey(hashVal))continue; // already was saved...
			sw.getBuffer().setLength(0);

			vals = val.get(hashVal);
			for(int n = 0; n < outputs; n++){
				Map.Entry<Integer,String> foundValue = null;
				for(Map.Entry<Integer,String> p:vals){
					if(p.getKey() != n)continue; // we need values for this output
					foundValue = p;
					break;
				}

				if(foundValue != null){
					sw.append(foundValue.getValue());
					vals.remove(foundValue);
				}
				else //
					if( (value = dtExpValues.getImplicitValue(mol, n)) != null)
						sw.append(value); //

				if(n != outputs -1)sw.append(",");
			}
			if(vals.size() == 0)val.remove(hashVal); // no new values to save for this descriptor
			bi.write(""+ mol+ "\t" + dtDescriptors.getRawRow(mol).attachments.get(QSPRConstants.MOLECULE_ID_NO_STEREOCHEM ) + "\t"  + dtDescriptors.getRawRow(mol).attachments.get(QSPRConstants.MOLECULE_ID_STEREOCHEM ) + "\t" + hashVal+"\t"+sw.toString()+"\n");
			bw.write(sw.toString());
			writeLibSvmDescriptorLine(dtDescriptors, mol, bw); 
			bw.write("\n");
		}
		bi.close();
		bw.close();
	}
	 */
	public static void writeLibSvmTestSet(DescriptorsTable dtDescriptors, String filename) throws Exception
	{
		File f = new File(filename);
		if(f.exists())return;

		FileWriter fis = new FileWriter(f);
		BufferedWriter bw = new BufferedWriter(fis,QSPRConstants.FILE_BUFFER);

		for (int mol = 0; mol < dtDescriptors.getDataSize(); mol++)
		{
			bw.write("0");
			writeLibSvmDescriptorLine(dtDescriptors, mol, bw);
			bw.write("\n");
		}
		bw.close();
	}

	/**
	 *  Changed to work with float values to speed up saving of files
	 * @param dtDescriptors
	 * @param conf
	 * @param out
	 * @throws Exception
	 */

	public static void writeLibSvmDescriptorLine(DescriptorsTable dtDescriptors, int mol, Writer out) throws IOException
	{
		double vals[] = dtDescriptors.getScaledValues(mol);

		int i = 1;
		for (double v : vals)
		{
			if (v != 0)
				out.write(" " + i + ":" +  NumericalValueStandardizer.getSignificantDigits(v));
			i++;
		}

		if( mol == 0 && vals[vals.length - 1] == 0) // we need to add zero just to indicate length of data for SVM
			out.write(" " + vals.length + ":0");
	}

}
