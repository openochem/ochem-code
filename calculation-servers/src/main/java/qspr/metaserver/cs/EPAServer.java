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

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.Vector;

import com.eadmet.utils.OCHEMUtils;

import qspr.dao.ChemInfEngine;
import qspr.dao.Various;
import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.metaserver.configurations.DescriptorsEPAConfiguration;
import qspr.workflow.datatypes.DataTable;

public class EPAServer extends CDKDescriptorsServer{

	public static EPAServer working;

	String epa[];


	public EPAServer(){
		supportedTaskType = DescriptorsConfiguration.EPA;
		repostSize = 10;
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
		working = this;
	}

	@SuppressWarnings("unchecked")
	@Override
	protected
	DataTable calculateDescriptors(final DataTable dtMolecules, DescriptorsAbstractConfiguration config, int start,
			int batchSize) throws Exception {

		if(epa == null) {
			if(Various.defaultEngine ==  ChemInfEngine.CDK) { 
				epa = OCHEMUtils.loadJars(getExeFile());
				setStatus("Using libraries " + epa);
			}
		}

		DescriptorsEPAConfiguration conf = (DescriptorsEPAConfiguration) config;

		final DataTable dtResults = getResults();

		if(dtResults.getColumnsSize()==0) {

			Class<?> clazz = Class.forName("gov.epa.ccte.apps.modelingservices.testDescriptors.domain.Descriptors.DescriptorFactory.DescriptorData");
			java.lang.reflect.Method method = clazz.getMethod("getDescriptorNames", boolean.class);
			Vector<String> vals = (Vector<String>)method.invoke(null,conf.requires3D());

			//Vector<String> vals = DescriptorData.getDescriptorNames(conf.requires3D());
			for(String v: vals)
				dtResults.addColumn(v);
			System.out.println(dtResults.getColumnsSize());
		}

		for(int i = start; i < start+batchSize; i++)
		{
			dtResults.addRow();
			final int mol = i; // required to be used in run()

			Thread t = new Thread()
			{
				@Override
				public void run()
				{
					try 
					{
						setStatus("Processing molecule " + (mol + 1) + " out of " + dtMolecules.getRowsSize());


						Class<?> clazz = Class.forName("gov.epa.ccte.apps.modelingservices.testDescriptors.service.TestDescriptorService");
						java.lang.reflect.Method method =  conf.requires3D()? 
								clazz.getMethod("goDescriptors2d_and_3dVector", String.class):
									clazz.getMethod("goDescriptors2dVector", String.class)
									;
						Vector<String> vals = (Vector<String>)method.invoke(null,(String)dtMolecules.getValue(mol, 0));

						//Vector<String> vals = conf.requires3D() ? 
						//		TestDescriptorService.goDescriptors2d_and_3dVector((String)dtMolecules.getValue(mol, 0)):
						//			TestDescriptorService.goDescriptors2dVector((String)dtMolecules.getValue(mol, 0))
						//			;

						if(vals.size() != dtResults.getColumnsSize()) {
							dtResults.getCurrentRow().setError("Number of calculated descriptors = " + vals.size() + "  != expected = " + dtResults.getColumnsSize());
						}else {
							double val = 0;
							for(int i=0; i<vals.size();i++) {
								val += Math.abs(Double.valueOf(vals.get(i)));
								dtResults.setValue(i, vals.get(i));
							}
							if(val == 0)
								dtResults.getCurrentRow().setError("all values are 0, calculations failed");

						}
					}
					catch(Throwable e){ // Runtime exception due to the timeout is here ...
						e.printStackTrace();
						StringWriter errors = new StringWriter();
						e.printStackTrace(new PrintWriter(errors));
						dtResults.getCurrentRow().setError(e.getMessage() == null?errors.toString():e.getMessage());
						setStatus("Molecule " + mol + " failed for calculation " + e.getMessage() + " " + errors.toString());
					}
				}

			};
			t.start();
			t.join(conf.getTimeoutInSeconds()*1000);
		}

		return dtResults;
	}

	public static String getFile(String filename) {
		filename = filename.replace(" ", "_");
		return working.getAliasedFileName(filename);
	}


}
