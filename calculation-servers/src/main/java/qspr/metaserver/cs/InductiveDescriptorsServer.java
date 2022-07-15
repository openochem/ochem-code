
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

import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.metaserver.configurations.DescriptorsInductiveDescriptorsConfiguration;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;

/**
 * 'Inductive' descriptors calculation server.
 * 
 * Based mainly on publications of Artem Cherkasov.
 * These descriptors require optimized geometries.
 * 
 * @author mrupp
 */
public class InductiveDescriptorsServer extends DescriptorsAbstractServer 
{

	public InductiveDescriptorsServer()
	{
		supportedTaskType = DescriptorsConfiguration.INDUCTIVE;
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
	}

	private static final String[] colNamesMElectronegativity = {"EoEqualized", "AverageEoPos", "AverageEoNeg"};
	private static final String[] colNamesMHardness          = {"GlobalHardness", "SumHardness", "SumPosHardness", "SumNegHardness", "AverageHardness", "AveragePosHardness", "AverageNegHardness", "SmallestPosHardness", "SmallestNegHardness", "LargestPosHardness", "LargestNegHardness", "HardnessOfMostPos", "HardnessOfMostNeg"};
	private static final String[] colNamesMSoftness          = {"GlobalSoftness", "TotalPosSoftness", "TotalNegSoftness", "AverageSoftness", "AveragePosSoftness", "AverageNegSoftness", "SmallestPosSoftness", "SmallestNegSoftness", "LargestPosSoftness", "LargestNegSoftness", "SoftnessOfMostPos", "SoftnessOfMostNeg"};
	private static final String[] colNamesMPartialCharges    = {"TotalCharge", "TotalChargeFormal", "AveragePosCharge", "AverageNegCharge", "MostPosCharge", "MostNegCharge"};
	private static final String[] colNamesMInductiveParam    = {"TotalSigmaMolI", "TotalAbsSigmaMolI", "MostPosSigmaMolI", "MostNegSigmaMolI", "MostPosSigmaIMol", "MostNegSigmaIMol", "SumPosSigmaMolI", "SumNegSigmaMolI"};
	private static final String[] colNamesMStericParam       = {"LargestRsMolI", "SmallestRsMolI", "LargestRsIMol", "SmallestRsIMol", "MostPosRsMolI", "MostNegRsMolI", "MostPosRsIMol", "MostNegRsIMol", "LargestPosRsMolI", "SmallestNegRsMolI", "LargestPosRsIMol", "SmallestNegRsIMol"};
	public static final int MAXATOMS = 500;

	@Override
	public WorkflowNodeData calculateDescriptors(WorkflowNodeData wnData, DescriptorsAbstractConfiguration receivedConfiguration) throws Exception
	{
		DataTable dtMolecules = wnData.ports.get(0);
		DataTable dtResult = new DataTable(true);
		dtResult.id = "descriptors";

		DescriptorsInductiveDescriptorsConfiguration configuration =  (DescriptorsInductiveDescriptorsConfiguration) receivedConfiguration;
		out.println("Computing 'inductive' descriptors...");
		out.println("Selected options: " + configuration.getSummary());

		// Create columns.
		if(configuration.calcMElectronegativity) for(final String s : colNamesMElectronegativity) dtResult.addColumn(s);
		if(configuration.calcMHardness         ) for(final String s : colNamesMHardness         ) dtResult.addColumn(s);
		if(configuration.calcMSoftness         ) for(final String s : colNamesMSoftness         ) dtResult.addColumn(s);
		if(configuration.calcMPartialCharges   ) for(final String s : colNamesMPartialCharges   ) dtResult.addColumn(s);
		if(configuration.calcMInductiveParam   ) for(final String s : colNamesMInductiveParam   ) dtResult.addColumn(s);
		if(configuration.calcMStericParam      ) for(final String s : colNamesMStericParam      ) dtResult.addColumn(s);

		out.println("We will calculate " + dtResult.getColumnsSize() + " inductive descriptors...");

		int okCount = 0;
		dtMolecules.reset();
		while (dtMolecules.nextRow())
		{
			setStatus("Processing molecule " + (dtMolecules.currentRow + 1) + " out of " + dtMolecules.getRowsSize());
			dtResult.addRow();
			//for debugging purpose 
			//			setStatus("Size of columns "+dtMolecules.getColumnsSize()+"\n");
			//			setStatus(String.valueOf(dtMolecules.getCurrentRow().getAttachment(QSPRConstants.MOLECULE_ID)));
			//			if(dtMolecules.getColumnIndex(QSPRConstants.MOLECULE_ID) != -1)
			//				setStatus("molecule mapping id "+ dtMolecules.getValue(QSPRConstants.MOLECULE_ID)+" and the INCHI is "+dtMolecules.getValue(QSPRConstants.INCHI));
			try 
			{
				final String sdf = (String) dtMolecules.getValue();
				final InductiveDescriptorsCalculator.Result calc = InductiveDescriptorsCalculator.calculate(sdf); 

				int columnIndex = 0;

				if(configuration.calcMElectronegativity)
				{
					dtResult.setValue(columnIndex++, (double) calc.mEoEqualized);
					dtResult.setValue(columnIndex++, (double) calc.mAverageEoPos);
					dtResult.setValue(columnIndex++, (double) calc.mAverageEoNeg);
				}

				if(configuration.calcMHardness)
				{
					dtResult.setValue(columnIndex++, (double) calc.mGlobalHardness     );
					dtResult.setValue(columnIndex++, (double) calc.mSumHardness        );
					dtResult.setValue(columnIndex++, (double) calc.mSumPosHardness     );
					dtResult.setValue(columnIndex++, (double) calc.mSumNegHardness     );
					dtResult.setValue(columnIndex++, (double) calc.mAverageHardness    );
					dtResult.setValue(columnIndex++, (double) calc.mAveragePosHardness );
					dtResult.setValue(columnIndex++, (double) calc.mAverageNegHardness );
					dtResult.setValue(columnIndex++, (double) calc.mSmallestPosHardness);
					dtResult.setValue(columnIndex++, (double) calc.mSmallestNegHardness);
					dtResult.setValue(columnIndex++, (double) calc.mLargestPosHardness );
					dtResult.setValue(columnIndex++, (double) calc.mLargestNegHardness );
					dtResult.setValue(columnIndex++, (double) calc.mHardnessOfMostPos  );
					dtResult.setValue(columnIndex++, (double) calc.mHardnessOfMostNeg  );
				}

				if(configuration.calcMSoftness)
				{
					dtResult.setValue(columnIndex++, (double) calc.mGlobalSoftness     );
					dtResult.setValue(columnIndex++, (double) calc.mTotalPosSoftness   );
					dtResult.setValue(columnIndex++, (double) calc.mTotalNegSoftness   );
					dtResult.setValue(columnIndex++, (double) calc.mAverageSoftness    );
					dtResult.setValue(columnIndex++, (double) calc.mAveragePosSoftness );
					dtResult.setValue(columnIndex++, (double) calc.mAverageNegSoftness );
					dtResult.setValue(columnIndex++, (double) calc.mSmallestPosSoftness);
					dtResult.setValue(columnIndex++, (double) calc.mSmallestNegSoftness);
					dtResult.setValue(columnIndex++, (double) calc.mLargestPosSoftness );
					dtResult.setValue(columnIndex++, (double) calc.mLargestNegSoftness );
					dtResult.setValue(columnIndex++, (double) calc.mSoftnessOfMostPos  );
					dtResult.setValue(columnIndex++, (double) calc.mSoftnessOfMostNeg  );
				}

				if(configuration.calcMPartialCharges)
				{
					dtResult.setValue(columnIndex++, (double) calc.mTotalCharge      );
					dtResult.setValue(columnIndex++, (double) calc.mTotalChargeFormal);
					dtResult.setValue(columnIndex++, (double) calc.mAveragePosCharge );
					dtResult.setValue(columnIndex++, (double) calc.mAverageNegCharge );
					dtResult.setValue(columnIndex++, (double) calc.mMostPosCharge    );
					dtResult.setValue(columnIndex++, (double) calc.mMostNegCharge    );
				}

				if(configuration.calcMInductiveParam)
				{
					dtResult.setValue(columnIndex++, (double) calc.mTotalSigmaMolI      );
					dtResult.setValue(columnIndex++, (double) calc.mTotalAbsSigmaMolI   );
					dtResult.setValue(columnIndex++, (double) calc.mMostPosSigmaMolI    );
					dtResult.setValue(columnIndex++, (double) calc.mMostNegSigmaMolI    );
					dtResult.setValue(columnIndex++, (double) calc.mMostPosSigmaIMol    );
					dtResult.setValue(columnIndex++, (double) calc.mMostNegSigmaIMol    );
					dtResult.setValue(columnIndex++, (double) calc.mSumPosSigmaMolI     );
					dtResult.setValue(columnIndex++, (double) calc.mSumNegSigmaMolI     );
				}

				if(configuration.calcMStericParam)
				{
					dtResult.setValue(columnIndex++, (double) calc.mLargestRsMolI    );
					dtResult.setValue(columnIndex++, (double) calc.mSmallestRsMolI   );
					dtResult.setValue(columnIndex++, (double) calc.mLargestRsIMol    );
					dtResult.setValue(columnIndex++, (double) calc.mSmallestRsIMol   );
					dtResult.setValue(columnIndex++, (double) calc.mMostPosRsMolI    );
					dtResult.setValue(columnIndex++, (double) calc.mMostNegRsMolI    );
					dtResult.setValue(columnIndex++, (double) calc.mMostPosRsIMol    );
					dtResult.setValue(columnIndex++, (double) calc.mMostNegRsIMol    );
					dtResult.setValue(columnIndex++, (double) calc.mLargestPosRsMolI );
					dtResult.setValue(columnIndex++, (double) calc.mSmallestNegRsMolI);
					dtResult.setValue(columnIndex++, (double) calc.mLargestPosRsIMol );
					dtResult.setValue(columnIndex++, (double) calc.mSmallestNegRsIMol);
				}

				++okCount;
			}
			catch(InductiveDescriptorsCalculator.ComputationFailureException e) 
			{
				dtResult.getCurrentRow().setError(e.getMessage());
				out.println(String.format("Calculation failure: %s", e.getMessage()));
				if(e.sdf != null) out.println(String.format("SDF:\n%s\n", e.sdf));
				continue;
			}catch(Exception e) {
				dtResult.getCurrentRow().setError(e.getMessage());
				continue;
			}

		}

		out.println("'Inductive' descriptors calculated for " + okCount + " out of " + dtMolecules.getRowsSize() + " molecules.");

		return new WorkflowNodeData(dtResult);
	}

	/*
	public static void main(String[] args) throws JAXBException, Exception
	{
		// Just a test for debugging. Delete it later.
		final DataTable dt = DataTable.fromXml(JAXBContextFactory.get("qspr.workflow.datatypes:qspr.workflow.structure"), "/ews/WorkflowSamples/files/test-samples/datatables/dt_test_inductive_1.xml");

		InductiveDescriptorsConfiguration conf = new InductiveDescriptorsConfiguration();
		conf.calcAElectronegativity = true;
		conf.calcMHardness = false;
		Task task = new Task(DescriptorsConfiguration.INDUCTIVE, conf, new WorkflowNodeData(dt));
		task.setPreferredServer("ibis60-52.EclipseCSServer");

		task = new CalculationClient("Inductive Test").calculateTask(task);
		task.check();
	}
	 */
}
