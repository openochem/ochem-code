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

import qspr.workflow.utils.QSPRConstants;

import com.eadmet.utils.NumericalValueStandardizer;

public class LibSvmModel extends AbstractSvmModel
{
	private static final long serialVersionUID = 1L;

	List<String> headers = new ArrayList<String>();
	float[][] vector_coeficients;
	short[][] indices;
	float[][] vectors;

	public void readModelFromFile(String filename) throws IOException
	{
		FileReader fis = new FileReader(filename);
		BufferedReader model = new BufferedReader(fis,QSPRConstants.FILE_BUFFER);

		String line;
		List<short[]> tmpIndices = new ArrayList<short[]>();
		List<float[]> tmpVectors = new ArrayList<float[]>();
		List<float[]> tmpVectorCoefs = new ArrayList<float[]>();
		boolean dataSection = false;

		try{
			while ((line = model.readLine()) != null)
			{
				if (!dataSection)
				{
					headers.add(line);
					dataSection = line.equals("SV");
				}
				else
				{
					String[] pieces = line.split("\\s+");
					short[] indices = null;
					float[] coefs = null;
					float[] weights = null;

					for (int i=0; i<pieces.length; i++) //Ugly workaround, probably possible to get number of classes from model headers
					{
						String[] subpieces = pieces[i].split(":");
						if (subpieces.length == 2) //We are already vectors
						{
							coefs = new float[i];
							weights = new float[pieces.length-i];
							indices = new short[pieces.length-i];
							break;
						}
					}

					if (coefs == null) //Funny fix for empty support vectors
					{
						coefs = new float[1];
						weights = new float[0];
						indices = new short[0];
					}

					for (int i=0; i<pieces.length; i++)
					{
						if (i < coefs.length)
							coefs[i] = new Float(pieces[0]);
						else
						{
							String[] subpieces = pieces[i].split(":");
							int column = Integer.parseInt(subpieces[0]);
							if(column>Short.MAX_VALUE)throw new IOException("Too large value");
							indices[i-coefs.length] = (short)column;
							weights[i-coefs.length] = new Float(subpieces[1]);
						}
					}
					tmpIndices.add(indices);
					tmpVectors.add(weights);
					tmpVectorCoefs.add(coefs);
				}
			}
		}catch(IOException e){
			throw new IOException(e.getMessage());
		}finally{
			model.close();
		}

		vector_coeficients = new float[tmpVectorCoefs.size()][];
		for (int i=0; i<tmpVectorCoefs.size(); i++)
			vector_coeficients[i] = tmpVectorCoefs.get(i);

		indices = new short[tmpIndices.size()][];
		for (int i=0; i<tmpIndices.size(); i++)
			indices[i] = tmpIndices.get(i);

		vectors = new float[tmpVectors.size()][];
		for (int i=0; i<tmpVectors.size(); i++)
			vectors[i] = tmpVectors.get(i);
	}

	public void writeModelToFile(String filename) throws IOException
	{
		FileWriter fis = new FileWriter(filename);
		BufferedWriter model = new BufferedWriter(fis,QSPRConstants.FILE_BUFFER);

		try{
			for (int i=0; i<headers.size(); i++)
				model.write(headers.get(i)+"\n");

			for (int i=0; i<vector_coeficients.length; i++)
			{
				for (int j=0; j<vector_coeficients[i].length; j++)
				{
					if (j != 0)
						model.write(" ");
					model.write(""+NumericalValueStandardizer.formattedFloatValuesForWEKA(vector_coeficients[i][j]));
				}

				for (int j=0; j<indices[i].length; j++)
				{
					model.write(" "+indices[i][j]+":"+NumericalValueStandardizer.formattedFloatValuesForWEKA(vectors[i][j]));
				}
				model.write("\n");
			}
		}catch(IOException e){
			throw new IOException(e.getMessage());
		}finally{
			model.close();
		}
	}

}

