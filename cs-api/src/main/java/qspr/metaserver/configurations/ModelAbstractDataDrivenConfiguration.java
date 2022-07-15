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

import javax.xml.bind.annotation.XmlTransient;

import qspr.interfaces.DataDrivenConfiguration;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;

public abstract class ModelAbstractDataDrivenConfiguration extends
ModelAbstractConfiguration implements DataDrivenConfiguration{
	private static final long serialVersionUID = 1L;

	public CompressedObject<WorkflowNodeData> trainingSet = new CompressedObject<WorkflowNodeData>();

	@XmlTransient
	public void setTrainingSet(WorkflowNodeData data) {
		trainingSet.set(data);
	}

	public DataTable getTrainingSetDescriptors(){
		WorkflowNodeData trainingWorkflow = trainingSet.get();
		return trainingWorkflow ==null ? null : trainingWorkflow.ports.get(0);
	}

	public DataTable getTrainingSetValues(){
		WorkflowNodeData trainingWorkflow = trainingSet.get();
		return trainingWorkflow == null ? null : trainingWorkflow.ports.get(1);
	}

	@Override
	public void removeTrainingSet()
	{
		trainingSet.set(null);
	}

	@Override
	public boolean isTrainingConfiguration()
	{
		return getTrainingSetDescriptors() == null && savedmodel == null; // saved Model is for multilearning
	}

	public boolean isModelSaved() {
		return getTrainingSetDescriptors() != null; 	
	}

}
