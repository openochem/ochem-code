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

package qspr.modelling;

import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;

import qspr.entities.ModelMapping;
import qspr.metaserver.configurations.ConsensusModelConfiguration;
import qspr.metaserver.protocol.Task;
import qspr.workflow.utils.QSPRConstants;

public class ConsensusModelProcessor extends BasicModelProcessor
{

	@Override
	public Serializable getApplierConfiguration(boolean recalculated, boolean forceRecalculateDescriptors) throws Exception
	{
		ConsensusModelConfiguration consConf = (ConsensusModelConfiguration) getModelConfiguration();
		consConf.units = new ArrayList<String>();
		for (ModelMapping mm : model.modelMappings){
			consConf.units.add(mm.unit.getName());
		}
		return consConf; 
	}

	@Override
	public Task createApplierTask() throws Exception
	{
		return new Task(QSPRConstants.CONSENSUS, null, null);
	}

	@Override
	protected Serializable getTeacherConfiguration() throws Exception
	{
		throw new IOException("Consensus model cannot be trained: it can be only created.");
	}

	@Override
	protected String getTaskType()
	{
		return QSPRConstants.CONSENSUS;
	}

}
