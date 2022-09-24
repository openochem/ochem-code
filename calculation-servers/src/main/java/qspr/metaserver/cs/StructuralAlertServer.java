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
import java.util.Arrays;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.dao.Various;
import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.metaserver.configurations.DescriptorsStructuralAlertsConfiguration;
import qspr.metaserver.protocol.Task;
import qspr.workflow.utils.QSPRConstants;
import qspr.metaserver.util.SMARTMatcher;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;
import cern.colt.Timer;

public class StructuralAlertServer extends DescriptorsAbstractServer {

	private static transient final Logger logger = LogManager.getLogger(StructuralAlertServer.class);

	public StructuralAlertServer() 
	{
		supportedTaskType = DescriptorsConfiguration.StructuralAlerts;
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
	}

	@Override
	public int getRepostSize(Serializable configuration)
	{
		// Define the reposting size based on the number of alerts
		DescriptorsStructuralAlertsConfiguration conf = (DescriptorsStructuralAlertsConfiguration) configuration;
		if (conf.size() == 0)
			return 0;

		return Math.min(200, 100000 / conf.size())/2;
	}

	@Override
	public WorkflowNodeData calculateDescriptors(WorkflowNodeData wnData, DescriptorsAbstractConfiguration receivedConfiguration) throws Exception 
	{

		if (receivedConfiguration == null || ! (receivedConfiguration instanceof DescriptorsStructuralAlertsConfiguration))
			throw new Exception("Invalid configuration passed, should be instance of StructuralAlertsConfiguration");

		DescriptorsStructuralAlertsConfiguration configuration = (DescriptorsStructuralAlertsConfiguration) receivedConfiguration;
		out.println(configuration.printConfig());

		Timer overall = new Timer(); overall.start();
		out.println("Preparing to run structural alerts computation...");

		// get the molecules from the input datatable
		DataTable dtMolecules = wnData.ports.get(0);

		// create datatable for the results
		DataTable result = new DataTable(true);	

		SMARTMatcher matcher = new SMARTMatcher(configuration.alertPatterns, Various.molecule.engine);
		//List<Molecule> moleculePattern = configuration.getAlertSetAsJChemMolecules(configuration.alertPatterns);
		if (configuration.compactMode)
			for (int i = 0; i < (configuration.size() - 1) / 32 + 1; i++)
				result.addColumn("AlertsBits" + i);
		else
			for (int i = 0; i < configuration.size(); i++)
				result.addColumn("Alert" + i + "_" + configuration.get(i));

		dtMolecules.reset();
		while (dtMolecules.nextRow())
		{
			int[] alertsBits = null;
			if (configuration.compactMode)
				alertsBits = new int[result.getColumnsSize()];
			try {
				setStatus("" + dtMolecules.currentRow + " of " + dtMolecules.getRowsSize() + " done");

				String sdf = preprocessSDF((String) dtMolecules.getValue());
				sdf = Various.molecule.convertToFormat(sdf, QSPRConstants.SDFAROM_GENERAL_WITHH);
				result.addRow();
				int failed = 0, count = 0;
				for (int i = 0; i < configuration.size(); i++) {
					String curPattern = configuration.get(i);
					try 
					{
						if(i == 0) matcher.getMatchCount(sdf, 0);

						if (configuration.compactMode)
						{
							if (matcher.matchPattern(sdf, i)) {
								alertsBits[i / 32] |= (1 << (i%32));
								count++;
							}
						}
						else {
							int number =  matcher.getMatchCount(sdf, i);
							count += number;
							result.setValue("Alert" + i + "_" + curPattern, number);
						}
					} catch (Exception e) 
					{
						failed++;
					}
				}

				if(failed>0 && count == 0) { // molecule was failing 
					if(!configuration.compactMode)
						result.getCurrentRow().setError("Molecule failed for all alerts");
					else
						for(int i=0;i<configuration.size();i++)
							if(configuration.get(i).contains(QSPRConstants.ERROR))
								alertsBits[i / 32] |= (1 << (i%32));
				}

				if (configuration.compactMode)
					for (int i = 0; i < alertsBits.length; i++)
						result.setValue(i, Float.intBitsToFloat(alertsBits[i]));

			} catch (Exception e) {
				result.getCurrentRow().setError(supportedTaskType +":" + e.getMessage());
				e.printStackTrace(out);
				out.println(e.getMessage() + " for molecule " + dtMolecules.currentRow);
			}

		}

		overall.stop();
		out.println("We ran structural alert analysis in " + overall.seconds() + " seconds.");

		result.compact();

		out.println("Used memory: "+usedMemory()/(1024*1024)+"MB");
		return new WorkflowNodeData(result);
	}

	static private String trimRight(String st)
	{
		return st.replaceAll("[\\s|\\u00A0]+$", "");
	}

	static private String preprocessSDF(String sdf)
	{
		String st = trimRight(sdf);
		if (st.charAt(st.length() - 1) != '$')
			return st + "\n$$$$\n";
		else
			return st + "\n";
	}


	public static void main(String[] args) throws Exception
	{



		DescriptorsStructuralAlertsConfiguration conf = new DescriptorsStructuralAlertsConfiguration();
		conf.compactMode = false;
		String alerts[] = {"NOT C=C AND NC","NOT NCCC","CC","CCCCC","C=C AND NC AND CCC","CCF OR NCC","CCF OR NOT CCF","NOT CC OR NOT CCN"};
		conf.setAlerts(Arrays.asList(alerts));

		//		Molecule smart = conf.getAlertSetAsJChemMolecules(conf.alertPatterns).get(0);
		//logger.info(matchPattern(MolImporter.importMol(MolImporter.importMol("CCCCNCCCCCCC=C", "smiles").exportToFormat("sdf"), "sdf"), smart, false));
		StructuralAlertServer server = new StructuralAlertServer();
		//server.toolsDirectory = "/ews/WfTools/tools";

		Task task = new Task();
		String sdf = Various.molecule.convertToFormat("CCCCNCCCCCCC=C", QSPRConstants.SDF);

		DataTable dtMols = new DataTable();
		dtMols.addColumn(QSPRConstants.SDF_COLUMN);
		dtMols.addRow();
		dtMols.setValue(sdf);
		dtMols.addRow();
		dtMols.setValue(sdf);
		task.setData(new WorkflowNodeData(dtMols));


		task.setConfiguration(conf);
		server.calculate(task);
		task.check();

		WorkflowNodeData wnd = WorkflowNodeData.fromTask(task);
		for (int i = 0; i < wnd.ports.get(0).getColumnsSize(); i++)
		{
			logger.info( wnd.ports.get(0).getColumn(i) + "   " + ((Double)wnd.ports.get(0).getValue(0, i)));
		}
	}

	//	public static void main2(String[] args) throws Exception {
	//		try {
	//			MolSearch s = new MolSearch();
	//
	//			// queryMode = true forces string to be imported as SMARTS
	//			// If SMILES import needed, set queryMode = false.
	//			MolHandler mh1 = new MolHandler("[$([*R2]([*R])([*R])([*R]))].[$([*R2]([*R])([*R])([*R]))]", true);
	//
	//			// The query molecule must be aromatized if it uses the
	//			// alternating single/double bonds for the description of
	//			// aromaticity.
	//			mh1.aromatize();
	//			s.setQuery(mh1.getMolecule());
	//
	//			// use Molfile molecule as target
	//			//              BufferedInputStream tis=null;
	//			//              tis = new BufferedInputStream(new FileInputStream("/Volumes/data/work/aData/jcsearchTesting/molsdf.sdf"));
	//			//              MolInputStream tmis = new MolInputStream(tis);
	//			//              MolImporter tmolimp = new MolImporter(tmis);
	//			//Molecule target = tmolimp.read();
	//
	//			List<Molecule> mols = sdfImport("/Volumes/data/work/aData/jcsearchTesting/testAlerts.sdf");
	//
	//			for (Molecule target : mols) {
	//
	//
	//				target.aromatize(true);
	//				s.setTarget(target);
	//
	//				// search all matching substructures and print hits
	//				int[][] hits=null;
	//
	//				hits = s.findAll();
	//				if(hits==null)
	//					logger.info("No hits");
	//				else  {
	//					for(int i=0; i < hits.length; i++) {
	//						System.out.print("Hit " + (i+1) + ":  ");
	//						int[] hit = hits[i];
	//						for(int j=0; j < hit.length; j++) {
	//							System.out.print(hit[j]+" ");
	//						}
	//						logger.info("\n");
	//					}
	//				}//end else
	//			}
	//		} catch (IOException e) {
	//			e.printStackTrace();
	//			System.exit(1);
	//		} catch (Exception e) {
	//			e.printStackTrace();
	//			System.exit(1);
	//		}//end catch
	//		//		      }//end main
	//		//		  }//end searchTest
	//
	//		//		logger.info(matchPattern(test, smart));
	//		//		logger.info(matchPattern(smart, test));
	//	}


}
