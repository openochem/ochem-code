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

import java.util.HashMap;
import java.util.Map;

import qspr.dao.Various;
import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.configurations.DescriptorsQNPRConfiguration;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.exceptions.UserFriendlyException;

public class QNPRServer extends DescriptorsAbstractServer
{
	@Override
	public WorkflowNodeData calculateDescriptors(WorkflowNodeData dtInput, DescriptorsAbstractConfiguration receivedConfiguration) throws Exception
	{

		// Test configuration passed, parse it
		if (receivedConfiguration == null || !(receivedConfiguration instanceof DescriptorsQNPRConfiguration))
			throw new UserFriendlyException("Invalid configuration passed, should be instance of QNPRConfiguration");

		DescriptorsQNPRConfiguration conf = (DescriptorsQNPRConfiguration)receivedConfiguration;
		int min_length = conf.minFragmentLength;
		int max_length = conf.maxFragmentLength;

		if (min_length > max_length)
			throw new UserFriendlyException("Invalid configuration passed : Minimum fragment length > Maximum fragment length");


		DataTable dtMolecules = dtInput.ports.get(0);
		DataTable dtResults = new DataTable(true);

		for(int i=0;i<dtMolecules.getRowsSize();i++)
		{			
			String sdf = (String)dtMolecules.getValue(i,0);	
			//sdf = Various.molecule.convertToFormat(sdf, QSPRConstants.SDFAROM_BASIC_WITHH); // is not required
			dtResults.addRow();

			if (sdf == null || sdf.trim().length() == 0) {
				dtResults.getCurrentRow().setError("SDF was not provided.");
				continue;				
			}

			String smiles=null;
			try {
				// Workaround to enable use of coordination bonds with QNPRF

				smiles = Various.molecule.convertToFormatFixMetal(sdf,QSPRConstants.SMILESUniqueNoHAromatic);
				smiles = smiles.replaceAll("[0123456789]+", "\u00A7");	

			} catch (Exception e) {
				dtResults.getCurrentRow().setError(e.getMessage());
				continue;
			}

			Map<String, Integer> descriptors = new HashMap<String, Integer>();

			if (smiles.length() <= min_length){
				dtResults.setValue(smiles,1);
				continue;
			}

			for (int len=min_length; len <= max_length; len++) 	    		
				for (int l=0; l<= smiles.length()-len; l++){
					String newdesc = smiles.substring(l,l+len);
					Integer n = descriptors.get(newdesc);
					if(n == null) n = 0;
					descriptors.put(newdesc, n+1); // increase
				}	

			for (Map.Entry<String, Integer> entry : descriptors.entrySet())
				dtResults.setValue(entry.getKey(),entry.getValue());
		}


		setStatus("Finished "+supportedTaskType);

		return new WorkflowNodeData(dtResults);
	}

	public QNPRServer()
	{
		supportedTaskType = "QNPR";
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
	}

}
