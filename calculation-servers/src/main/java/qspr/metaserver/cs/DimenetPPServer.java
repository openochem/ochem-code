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

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;

import com.eadmet.exceptions.CriticalException;

import qspr.dao.ChemInfEngine;
import qspr.dao.Various;
import qspr.metaserver.configurations.DIMENETConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.metaserver.cs.util.LabelsTable;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.utils.QSPRConstants;

public class DimenetPPServer extends SmilesOnlyAbstractServer{
	/** Initializes supported task and workflow. */
	public DimenetPPServer()
	{
		supportedTaskType = QSPRConstants.DIMENET;
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);	
		batchApplySize = 10000;
		startPositionResults = 0;
	}

	@Override
	void saveHeaders(BufferedWriter bw, DescriptorsTable dtDescriptors, LabelsTable dtExpValues, ModelAbstractConfiguration conf) throws IOException {
		super.saveHeaders(bw, dtDescriptors, dtExpValues, conf);

		BufferedWriter writer = getAliasedBufferedWriter(dtExpValues != null? "train.sdf":"apply.sdf");
		for(int i = 0; i< dtDescriptors.getDataSize(); i++) {
			AbstractDataRow r = dtDescriptors.getRawData().getRow(i);
			String sdf = (String)r.getAttachment(QSPRConstants.SDF_COLUMN);
			if(sdf == null)continue;
			writer.append(Various.molecule.addPropertyToSDF(sdf, QSPRConstants.SMILES_ATTACHMENT, (String)r.getAttachment(QSPRConstants.SMILES_ATTACHMENT)));
		}
		writer.close();
	}


	@Override
	public boolean isCritical(String message) {
		if(message.contains("Explicit valence for atom"))
			throw new CriticalException("Explicit valence for atom: " + message);
		if(message.contains("'NoneType' object has no attribute"))
			throw new CriticalException("'NoneType' object has no attribute: " + message);
		return super.isCritical(message);
	}
	
	@Override
	protected String getGPUCard(boolean forceCPU) { // always card 0, since environment is used to restrict only one card, which becomes 0
		return "" + (noGPU() || forceCPU ? NO_GPU : 0);
	}

	@Override
	protected String[] getMessages() {
		return new String[] { "Epoch","Train: tensor(","Valid: ["};
	}

	@Override
	void saveConfig(ModelAbstractConfiguration configuration, DescriptorsTable mols, boolean train, boolean forceCPU) throws IOException{

		Various.molecule = Various.getCheminfImpl(ChemInfEngine.CHEMAXON);

		DIMENETConfiguration conf = (DIMENETConfiguration) configuration;

		conf.sanitize = true;

		File f = new File(getAliasedFileName(CFG));
		if(f.exists()) {
			File ff = new File(getAliasedFileName(CFG+".train"));
			if(!ff.exists())
				f.renameTo(ff);
		}

		BufferedWriter writer = getAliasedBufferedWriter(CFG);

		writer.write("[Task]\n"+
				"train_mode = "+  (train?"True":"False")+ "\n"+
				(train?"\ntrain_data_file = " + DATAFILE:"apply_data_file = " + APPLYFILE) + "\n"+
				"result_file = " + PREDICTIONS + "\n" +
				"model_file = "+ MODEL + "\n");

		writer.write("\n\n[Details]\n"
				+   "seed = "+ conf.getSeed() + "\n"
				+ 	"gpu = "+ getGPUCard(forceCPU)   + "\n"
				+   "batch = " + conf.batch  + "\n"
				+   "nbepochs = "  + conf.nbepochs  + "\n"
				+ 	"canonize = False" + "\n"
				+	"early = "  + (conf.early?"True":"False")+ "\n" + "\n"
				+ (conf.isExternal3D()?"external3D = True":"") + "\n"
				);

		writer.write("\n");
		writer.close();
		return;
	}

}
