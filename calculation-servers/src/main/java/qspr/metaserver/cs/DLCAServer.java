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
import java.io.IOException;

import qspr.metaserver.configurations.DLCAConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.workflow.utils.QSPRConstants;


/**
 *  Requires specific handling, since instead of descriptors SDF data are provided
 * @author itetko
 *
 */

public class DLCAServer extends SmilesOnlyAbstractServer
{

	/** Initializes supported task and workflow. */
	public DLCAServer()
	{
		supportedTaskType = QSPRConstants.DLCA;
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);
		batchApplySize = QSPRConstants.MODEL_REPOST_SIZE;
		startPositionResults = 0;
	}

	@Override
	void saveConfig(ModelAbstractConfiguration configuration, DescriptorsTable mols, boolean traningMode, boolean forceCPU) throws IOException{

		DLCAConfiguration conf = (DLCAConfiguration) configuration;	
		BufferedWriter writer = getAliasedBufferedWriter(CFG);
		writer.write("[Task]\n");
		writer.write("\ntrain_mode = " + (traningMode?"True":"False"));
		writer.write("\nmodel_file = " + MODEL);
		writer.write("\ntrain_data_file = " + DATAFILE);
		writer.write("\napply_data_file = " + APPLYFILE);

		writer.write("\n\n[Details]");
		writer.write("\ngpu = "+ getGPUCard(forceCPU));
		writer.write("\nseed = " +  conf.getSeed() );
		writer.write("\nn_epochs = " + conf.epochs);
		writer.write("\nbatch_size = " + conf.batch);
		writer.write("\n");
		writer.close();
	}

}
