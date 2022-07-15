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

import java.io.Serializable;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.metaserver.configurations.LeverageConfiguration;
import qspr.workflow.utils.QSPRConstants
;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;
import Jama.Matrix;

// The standard DM calculation:
// (Descriptors, Structures) -> (DM values, DM model)
public class LeverageServer extends WorkflowNodeServer
{
	private static transient final Logger logger = LogManager.getLogger(LeverageServer.class);
	
	@Override
	public WorkflowNodeData calculate(WorkflowNodeData dtInput, Serializable configuration) throws Exception
	{
		if(configuration==null)return calculateIt(dtInput,configuration);

		LeverageConfiguration c=(LeverageConfiguration)configuration;
		DataTable dtDescriptors=dtInput.ports.get(0);
		
		addCorrectionDescriptor(dtDescriptors); // as mentioned by Simona

		return apply(dtDescriptors,c.invMatrix);		
	}
	private WorkflowNodeData calculateIt(WorkflowNodeData dtInput, Serializable configuration) throws Exception
		{
		DataTable dtDescriptors=dtInput.ports.get(0);
		
		addCorrectionDescriptor(dtDescriptors); // as mentioned by Simona
		
		int nTrainSet = getTrainingSetSize(dtDescriptors);
		int nDesc = dtDescriptors.getColumnsSize();
		
		out.println("" + dtDescriptors.getRowsSize() + " compounds, " + nTrainSet + " thereof in the training set");
		
		// H = X * (Xt * X)^(-1) * Xt = X * COV^(-1) * Xt
		double[][] cov = new double[nDesc][nDesc];
		
		setStatus("Creating the covariation matrix");
		for (int i = 0; i < nDesc; i++)
			for (int k = 0; k <= i; k++)
			{
				double sum = 0;
				for (int j = 0; j < nTrainSet; j++)
					sum += (Double)dtDescriptors.getValue(j, i) * (Double)dtDescriptors.getValue(j, k);
				cov[i][k] = cov[k][i] = sum;
			}
		
		setStatus("Inverting the covariation matrix");
		Matrix covMatrix = new Matrix(cov);
		logger.info("Determinant of cov matrix: " + covMatrix.det());
		
		Matrix invMatrix = covMatrix.inverse();
		
		double[][] invArray = invMatrix.getArray();
		
		return apply(dtDescriptors,invArray);
	}
	
	WorkflowNodeData apply(DataTable dtDescriptors, double[][] invMatrix) throws Exception
	{
		setStatus("Calculating leverages");
		addCorrectionDescriptor(dtDescriptors);
		
		DataTable dtLeverages = new DataTable(true);
		dtLeverages.id = "LeverageValues";
		dtLeverages.addColumn("Leverage");
		
		int nMols = dtDescriptors.getRowsSize();
		int nDesc = dtDescriptors.getColumnsSize();
		
		for (int mol = 0; mol < nMols; mol++)
		{
			if (mol % 10 == 0)
				setStatus("Calculating leverages: " + mol + " out of " + nMols);
			double sum = 0;
			for (int i = 0; i < nDesc; i++)
				for (int k = 0; k < nDesc; k++)
					sum += (Double)dtDescriptors.getValue(mol, i) * invMatrix[i][k] * (Double)dtDescriptors.getValue(mol, k);
			
			if (sum < 0)
				out.println("WARNING: Negative leverage " + sum  +", something must be wrong");
			dtLeverages.addRow();
			dtLeverages.setValue(sum);
		}
		
		LeverageConfiguration c= new LeverageConfiguration();
		c.invMatrix=invMatrix;
		
		return new WorkflowNodeData(dtLeverages).addPort(new DataTable(c));
	}
	
	// The method adds special unity descriptor (always "1")
	// which is necessary for the correct calculation of Leverage (mentioned by Simona, sounds reasonable)
	private void addCorrectionDescriptor(DataTable dtDescriptors)
	{
		if (dtDescriptors.containsColumn("CORRECTION-1"))
			return;
		dtDescriptors.addColumn("CORRECTION-1");
		for (int i = 0; i < dtDescriptors.getRowsSize(); i++)
			dtDescriptors.setValue(i, dtDescriptors.getColumnsSize() - 1, 1.0d);
	}

	
	private int getTrainingSetSize(DataTable dtDescriptors)
	{
		int res = 0;
		dtDescriptors.reset();
		while (dtDescriptors.nextRow())
			if (dtDescriptors.getCurrentRow().getAttachment(QSPRConstants.VALIDATION) == null)
				res++;
			else
				break;
		return res;
	}
	
	public LeverageServer()
	{
		supportedTaskType = "Leverage";
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
	}
	
}
