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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import qspr.workflow.datatypes.CompactDataRow;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.utils.NumericalValueStandardizer;

public class SvmModel extends AbstractSvmModel{

	private static final long serialVersionUID = 1L;

	DataTable data;
	String COEFF = "COEFF";
	String headers = "";

	/**
	 *  Prepares DataTable for storage of vectors from models
	 *  Allows faster storage of data
	 * @param model
	 * @throws NumberFormatException
	 * @throws IOException
	 */

	private int prepareDataTable(String filename) throws IOException {
		FileReader fis = new FileReader(filename);
		BufferedReader model = new BufferedReader(fis,QSPRConstants.FILE_BUFFER);

		data = new DataTable(true);

		Set<Integer> columns = new TreeSet<Integer>();

		String line;
		boolean dataSection = false;
		int coeff = 0;
		data = new DataTable(true);
		try{
			while ((line = model.readLine()) != null)
			{
				if (!dataSection)
					dataSection = line.equals("SV");
				else
				{
					boolean first = true;
					String[] pieces = line.split("\\s+");
					for (int i=0; i<pieces.length; i++)
					{
						String[] subpieces = pieces[i].split(":");
						if (subpieces.length == 2){ //Already vectors, storing columns
							columns.add(Integer.parseInt(subpieces[0]));
							if(coeff==0)coeff=i; // number of coefficients
							if(first && coeff != i)
								throw new IOException("Different number of coefficients in SVM per SVM vectors is not supported");
							first = false;
						}
					}
				}
			}
		}catch(IOException e){
			throw new IOException(e.getMessage());
		}finally{
			model.close();
			fis.close();
		}

		for(int i=0;i<coeff;i++)
			data.addColumn(COEFF+i);

		List<Integer> sortedList=new ArrayList<Integer>(columns);

		for(Integer i:sortedList)
			data.addColumn(""+i);

		return coeff;
	}


	public void readModelFromFile(String filename) throws IOException
	{

		int numberOfCoefficients = prepareDataTable(filename);

		FileReader fis = new FileReader(filename);
		BufferedReader model = new BufferedReader(fis,QSPRConstants.FILE_BUFFER);

		String line;
		boolean dataSection = false;
		try{
			while ((line = model.readLine()) != null)
			{
				if (!dataSection)
				{
					headers += line+"\n";
					dataSection = line.equals("SV");
				}
				else
				{
					String[] pieces = line.split("\\s+");
					if(pieces.length==0)continue;				
					data.addRow();
					for (int i=0; i<pieces.length; i++)
						if(i>=numberOfCoefficients)
						{
							String[] subpieces = pieces[i].split(":");
							data.setValue(subpieces[0], Float.parseFloat(subpieces[1]));
						}
						else
							data.setValue(COEFF+i,Float.parseFloat(pieces[i]));
				}
			}
		}catch(IOException e){
			throw new IOException(e.getMessage());
		}finally{
			model.close();
			fis.close();
		}
	}

	public void writeModelToFile(String filename) throws Exception
	{
		FileWriter fis = new FileWriter(filename);
		BufferedWriter model = new BufferedWriter(fis,QSPRConstants.FILE_BUFFER);
		try{
			model.write(headers);
			data.reset();
			String col;
			while(data.nextRow()){
				CompactDataRow row = (CompactDataRow) data.getCurrentRow();
				float val[] = row.toArray();
				for(int i = 0; i < val.length;i++)
					if((col=data.getColumn(i)).contains(COEFF))
						model.write(""+val[i]+" ");
					else
						if(val[i] != 0)model.write(col+":"+NumericalValueStandardizer.formattedFloatValuesForWEKA(val[i])+" ");
				model.write("\n");
			}
		}catch(IOException e){
			throw new IOException(e.getMessage());
		}finally{
			model.close();
			fis.close();
		}
	}
}