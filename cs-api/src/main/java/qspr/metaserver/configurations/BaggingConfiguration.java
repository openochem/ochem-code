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

import java.util.BitSet;
import java.util.List;
import java.util.Map;

import javax.xml.bind.annotation.XmlRootElement;

import qspr.workflow.datatypes.WorkflowNodeData;
import qspr.workflow.utils.QSPRConstants;

@XmlRootElement(name = "bagging-configuration")
public class BaggingConfiguration extends ValidationConfiguration
{
	private static final long serialVersionUID = 1L;
	
	public Integer numberInstances;
		
	public Boolean noOverrun;
	
	/*
	 * These fields are required for DataDriven configuration
	 */
	
	public CompressedObject<WorkflowNodeData> trainingSet = null;

	/**
	 * Contains IDs of rows which are NOT in the training set
	 */
	public CompressedObject<List<BitSet>> validationRowNumbers = null;
	/**
	 * Indicates number of repetitions for each row in the training set, which has at least 2 repetitions
	 */
	public CompressedObject<List<Map<Integer,Integer>>> trainingSetRepetions = null;
	
	/*
	 * End of fields required for data driven configuration
	 */
	
	/**
	 * Contains models used for bagging
	 */
	public List<Object> models;
	
	public int getInstances(){
		return numberInstances==null?0:numberInstances;
	}
	
	public String toString()
	{
		return "";
	}

	@Override
	public boolean isTrainingConfiguration()
	{
		return (models == null || models.isEmpty()) && savedmodel == null; // for multilearning!
	}
	
	@Override
	public void cleanBulkyStuff() {
		super.cleanBulkyStuff();
		if (models != null)
			models.clear();
	}

	@Override
	public String getDefaultName() {
		return QSPRConstants.BAGGING;
	}

	@Override
	public int requireMinimumRecords() {
		return 32;
	}
	
	public boolean isModelSaved() {
		return !isTrainingConfiguration();
	}
}
