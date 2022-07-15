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

import java.io.BufferedReader;
import java.io.IOException;

import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.workflow.datatypes.DataTable;

public class MOPAC2016Server extends MOPAC7Server{

	public MOPAC2016Server()
	{
		supportedTaskType = DescriptorsConfiguration.MOPAC2016;
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
		repostSize = 20;
	}

	@Override
	String getKeywords(){
		return PM7 + " " + KEYWORDS + " EF SUPER " + (gpuCard == NO_GPU?"NOGPU":"")+ " DISP THREADS="+THREADS;
	}

	@Override
	String getOutputFile(){
		return FILENAME+".out";
	}

	@Override
	protected void addCOSMO(BufferedReader reader, DataTable dtResult) throws IOException {
		dtResult.setValue("CosmoArea", parseValue("COSMO AREA", reader));
		dtResult.setValue("CosmoVolume", parseValue("COSMO VOLUME", reader));
	}

	@Override
	protected void addSuperDelocalizabilities(BufferedReader reader, DataTable dtResult) throws IOException{
		String line = skipTo("Mulliken electronegativity:", reader);
		double val;
		dtResult.setValue("MullikenElectronegativity",Double.parseDouble(last(line,1)));
		dtResult.setValue("AbsoluteHardness:", val=Double.parseDouble(last(read(reader),1)));
		dtResult.setValue("AbsoluteSoftness:", 0.5/(val+Double.MIN_VALUE));
		dtResult.setValue("SchuurmannMOShiftAlpha", Double.parseDouble(last(read(reader),1)));
		//skipTo("q(r) - Z(r)",reader);
		line = skipTo("Total:",reader);
		dtResult.setValue("TotalElectrophilicDelocalizability", Double.parseDouble(last(line,1)));
		dtResult.setValue("TotalNucleophilicDelocalizability", Double.parseDouble(last(line,2)));
		line = skipTo("Total:",reader);
		dtResult.setValue("TotalSelfPolarizability", Double.parseDouble(last(line,1)));
	}

	@Override
	protected void addAlphaPolar(BufferedReader reader, DataTable dtResult) throws IOException{
		skipTo("COMPONENTS OF ALPHA", reader);
		read(reader, 2);
		double aPol [] = new double[9];
		int i,j,n=0;
		for(i=0;i<3;i++){
			String line = read(reader);
			for(j=3;j>0;j--)
				aPol [n++]=Double.parseDouble(last(line,j));
		}

		double polar = parseValue("ISOTROPIC AVERAGE ALPHA", reader); 

		polar = polar * polar;

		double anisotropy = (aPol[0]*aPol[0] + aPol[4]*aPol[4] + aPol[8]*aPol[8] - 3 * polar)/(6 * polar);

		dtResult.setValue("PolarizabilityAnisotropy", anisotropy);
	}

}
