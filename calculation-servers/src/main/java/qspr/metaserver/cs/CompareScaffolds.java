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

import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.distribution.HypergeometricDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.metaserver.configurations.CompareScaffoldsConfiguration;
import qspr.workflow.utils.QSPRConstants
;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;

import com.eadmet.utils.NumericalValueStandardizer;

/**
 * Calculates the significantly overrepresented Descriptors
 * @author itetko
 *
 */

public class CompareScaffolds extends WorkflowNodeServer
{

	public static final double LIMIT = 0.000001;
	static public final String SAME_SET = "Ignored as duplicate within this set";
	static public final String PREVIOUS_SET = "Ignored as duplicate with the 1st set";
	static public final String pValue = "pValue1";

	private static transient final Logger logger = LogManager.getLogger(CompareScaffolds.class);

	public CompareScaffolds()
	{
		supportedTaskType = QSPRConstants.COMPARE_SCAFFOLDS;
		setInputFlowGroup(0, 1);
		setInputFlowGroup(1, 2);
		setOutputFlowGroup(1, 1);
		setOutputFlowGroup(2, 2);
	}

	@Override
	public WorkflowNodeData calculate(WorkflowNodeData wndInput, Serializable configuration) throws Exception
	{

		CompareScaffoldsConfiguration config = (CompareScaffoldsConfiguration) configuration;

		DataTable dataSet1 = wndInput.ports.get(0);
		DataTable dataSet2 = wndInput.ports.get(1);

		logger.info("Detect duplicates");

		findDuplicates(dataSet1, dataSet2);

		int nSet1 = dataSet1.getRowsNoErrorsSize();
		int nSet2 = dataSet2.getRowsNoErrorsSize();

		DataTable overrepresented = new DataTable();
		List <Double> order = new ArrayList<Double>();
		DataTable descriptors1 = dataSet1.replicateDataTable(0, dataSet1.getRowsSize(), true, true, false);
		DataTable descriptors2 = dataSet2.replicateDataTable(0, dataSet2.getRowsSize(), true, true, false);

		Set<String> descriptors = new HashSet<String>();

		if(dataSet1.getRowsSize()>dataSet2.getRowsSize()){ // to preserver order to speed-up search in a larger table
			descriptors.addAll(dataSet1.getColumns());
			descriptors.addAll(dataSet2.getColumns());
		}else{
			descriptors.addAll(dataSet2.getColumns());
			descriptors.addAll(dataSet1.getColumns());				
		}

		logger.info("Calculate descriptors");

		Map<String,Integer> N1=countDescriptors(dataSet1),N2=countDescriptors(dataSet2);

		logger.info("Detect significant descriptors");

		for (String descriptor : descriptors)
		{
			int n1 = N1.containsKey(descriptor)?N1.get(descriptor):0;
			int n2 = N2.containsKey(descriptor)?N2.get(descriptor):0;

			double pval = nSet2 == 0? 0.5: getPValue(n1, n2, nSet1, nSet2);
			if (hasEffect(pval, config.pvalue))
			{
				addDescriptors(descriptor, descriptors1, dataSet1);
				addDescriptors(descriptor, descriptors2, dataSet2);

				overrepresented.addRow();
				overrepresented.setValue("Descriptor", descriptor);
				overrepresented.setValue("C1", n1);
				overrepresented.setValue("C2", n2);
				overrepresented.setValue(pValue, pval);
				overrepresented.setValue("pValue2", -pval);
				order.add(pval!=0?Math.abs(pval):-n1);
			}
		}

		logger.info("Sort significant dscriptors");

		Collections.sort(order);

		DataTable sorted = overrepresented.getEmptyCopy();

		double previous=Double.MAX_VALUE;

		for(double v: order){
			if(v == previous)continue;
			previous = v;
			for(int i=0;i<overrepresented.getRowsSize();i++)
				if(Math.abs((Double)overrepresented.getValue(i, pValue)) == v)
					sorted.addRow(overrepresented.getRow(i));
		}
		

		WorkflowNodeData w = new WorkflowNodeData(sorted);
		w.addPort(descriptors1);
		w.addPort(descriptors2);

		logger.info("Finishing");

		return w;
	}

	private Map<Integer, Integer> nonRedundantIDs(DataTable dataSet) throws IOException
	{
		Map<Integer, Integer> mols = new HashMap<Integer, Integer>();

		for (int i = 0; i < dataSet.getRowsSize(); i++)
		{
			Integer mid = (Integer) dataSet.getRow(i).getAttachment(QSPRConstants.MOLECULE_ID_NO_STEREOCHEM);
			if (mid == null)
				throw new IOException("MID should be provided");
			if (mols.containsKey(mid))
				dataSet.getRow(i).setError(SAME_SET);
			mols.put(mid, i);
		}
		return mols;
	}

	private void findDuplicates(DataTable dataSet1, DataTable dataSet2) throws IOException
	{

		Map<Integer, Integer> mols1 = nonRedundantIDs(dataSet1);
		Map<Integer, Integer> mols2 = nonRedundantIDs(dataSet2);

		for (Integer mol : mols1.keySet())
			if (mols2.containsKey(mol))
				dataSet2.getRow(mols2.get(mol)).setError(PREVIOUS_SET);

	}

	static private double getPValue(int n1, int n2, int N1, int N2)
	{
		double val;

		HypergeometricDistribution hypergeom = new HypergeometricDistribution(N1 + N2, n1 + n2, N1);
		if (n1 > hypergeom.getNumericalMean())
			val =  hypergeom.upperCumulativeProbability(n1);
		else
			val = -hypergeom.cumulativeProbability(n1);

		return 	Double.parseDouble(NumericalValueStandardizer.getSignificantDigits(val));
	}

	private boolean hasEffect(double pval, double threshold)
	{
		return Math.abs(pval) < threshold;
	}

	private boolean descriptorExist(Double val)
	{
		return Math.abs(val) > LIMIT;
	}

	private void addDescriptors(String descriptor, DataTable descriptors, DataTable dataSet)
	{

		if (!dataSet.containsColumn(descriptor))
			return;

		int column = dataSet.indexOfColumn(descriptor);

		descriptors.setValue(0, descriptor, 0); // just to add new column

		int colnew = descriptors.indexOfColumn(descriptor);

		for (int i = 0; i < dataSet.getRowsSize(); i++)
			if(descriptorExist((Double) dataSet.getValue(i, column)))
				descriptors.setValue(i, colnew, 1);

	}


	private Map<String,Integer> countDescriptors(DataTable dataSet)
	{
		Map<String,Integer> m=new HashMap<String,Integer>();	
		int count[]=new int[dataSet.getColumnsSize()];

		for(int row=0;row<dataSet.getRowsSize();row++)
			if(!dataSet.getRow(row).isError()) // if error, we ignore this row for analysis
				for(int i=0;i<dataSet.getColumnsSize();i++)
					if (descriptorExist((Double) dataSet.getValue(row, i)))
						count[i]++;

		for(int i=0;i<dataSet.getColumnsSize();i++)
			if(count[i]>0)m.put(dataSet.getColumn(i), count[i]);

		return m;
	}

	public static void main(String[] args) throws Exception
	{

		HypergeometricDistribution hypergeom = new HypergeometricDistribution(100, 40, 20);

		System.out.println(" prob=" + hypergeom.upperCumulativeProbability(1));
		System.out.println(" prob=" + hypergeom.upperCumulativeProbability(10));
		System.out.println(" prob=" + hypergeom.upperCumulativeProbability(8));

		double val = getPValue(15, 5, 60, 40);
		System.out.println(" P1=" + val );
		val = getPValue(5, 15, 40, 60);
		System.out.println(" P1=" + val);

	}

}
