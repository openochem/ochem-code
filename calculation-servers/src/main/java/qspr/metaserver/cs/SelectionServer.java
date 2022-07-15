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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.metaserver.configurations.CompressedObject;
import qspr.metaserver.configurations.SelectionConfiguration;
import qspr.metaserver.protocol.DataSize;
import qspr.workflow.utils.QSPRConstants
;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.CompactDataRow;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.StringList;
import qspr.workflow.datatypes.WorkflowNodeData;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.FileUtils;
import com.eadmet.utils.NumericalValueStandardizer;

/*
 * Since UFS reuses a lot of code from my calculation servers, it is extension of my abstract classes
 * Better solution?
 */

public class SelectionServer extends ExecutableServer
{
	private static transient final Logger logger = LogManager.getLogger(SelectionServer.class);
	final static String OCHEM="ochem";
	final static String CFG="cfg";
	final static String DATAFILE = "data";

	public SelectionServer()
	{
		supportedTaskType = "Selection";
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);
	}


	@Override
	public WorkflowNodeData calculate(WorkflowNodeData wndInput,
			Serializable receivedConfiguration) throws Exception 
	{
		DataTable dtDescriptors = wndInput.ports.get(0);

		Set<String> initialDescriptors = new HashSet<>(dtDescriptors.columns);

		SelectionConfiguration selectionConfiguration = (SelectionConfiguration)receivedConfiguration;

		if (dtDescriptors.getColumnsSize() != 0) {  // we have descriptors...

			if (receivedConfiguration == null || !(receivedConfiguration instanceof SelectionConfiguration))
				throw new Exception("Invalid configuration: " + 
						receivedConfiguration + " must be instance of SelectionConfiguration");

			setStatus("Starting filtering process");

			if (!dtDescriptors.isCompactRowFormat()) throw new IOException("only CompactRowFormat is expected here!");

			if(selectionConfiguration.getDescriptorsSize() > 0){  // previously filtered or requested by user
				setStatus("Filtering by list");
				dtDescriptors.keepByList(selectionConfiguration.descriptorAsStrings());  // if there was a request to pre-filer data by a list; could be also a part of the training set configuration
			}

			// filtering block using cfg information
			if (selectionConfiguration.fixedDescriptors == null && // we do not have list of filtered descriptors, which is possible on apply model only
					selectionConfiguration.doFiltering()  && dtDescriptors.getColumnsSize() > 0 && guessTrainingSet(dtDescriptors) > 0) // we can have no descriptors, and this is OK if we use conditions as inputs
				filterIt(selectionConfiguration, dtDescriptors);
			else
				// Marking molecules with large values as errors, nothing else
				filterByMaximumValue(dtDescriptors,selectionConfiguration.maximumValueThreshold);

			setStatus("Rounding descriptors values");
			dtDescriptors.roundValues();
		}

		if (selectionConfiguration!= null && selectionConfiguration.fixedDescriptors != null) // application and we have conditions and we should add them
		{
			// Check if we really have received the correct fixed descriptors
			if (wndInput.ports.size() < 2)
				throw new UserFriendlyException("No obligatory conditions were provided but they are required for " + selectionConfiguration.fixedDescriptors);

			DataTable condTable = wndInput.ports.get(1);
			List<String> conditions = condTable.getColumns();

			setStatus("Updating condition values");

			for(String condition: selectionConfiguration.fixedDescriptors.values){
				if(!conditions.contains(condition))
					throw new UserFriendlyException("Obligatory condition \n: \"" + condition + "\" has not been provided. Has it been renamed?");
				if(dtDescriptors.getColumnIndex(condition) == -1)
					dtDescriptors.addColumn(condition); // adding conditions as descriptor
				int columnDescr = dtDescriptors.getColumnIndex(condition);
				int columnCondit = conditions.indexOf(condition);
				for(int row=0 ; row < dtDescriptors.getRowsSize() ; row++)
					dtDescriptors.setValue(row, columnDescr, condTable.getValue(row, columnCondit));
			}
		}
		else // training set, adding condition if required
			if (wndInput.ports.size() > 1)
			{
				setStatus("Merging with conditions: <" + wndInput.ports.get(1).id + ">");
				DataTable conditionsTable =  wndInput.ports.get(1);

				Set<String> conditionsToDelete = new HashSet<String>(); // removing already used conditions
				for(String col: conditionsTable.columns) {
					String name = (col+QSPRConstants.USED).toLowerCase();
					for(String desc: initialDescriptors)
						if(desc.contains(name)) {
							conditionsToDelete.add(desc);
							conditionsToDelete.add(col);
						}
				}

				dtDescriptors.mergeColumnsWith(conditionsTable);
				dtDescriptors.deleteByList(conditionsToDelete);

				List<String> conditionsToKeep = new ArrayList<String>();
				for(String name:conditionsTable.getColumns())
					if(!conditionsToDelete.contains(name))
						conditionsToKeep.add(name);

				if(conditionsToKeep.size()>0)
					selectionConfiguration.fixedDescriptors = new StringList(conditionsToKeep);
			}

		if (selectionConfiguration != null)
			selectionConfiguration.storeDescriptorAsStrings(dtDescriptors.getColumns()); // all kept descriptors

		setStatus("FILTERED: finished filtering process, number of selected descriptors:"+dtDescriptors.getColumnsSize());
		dtDescriptors.id = "descriptors";

		dtDescriptors.normalize();

		DataSize n = dtDescriptors.getValuesCount();

		setStatus("There are " + n.nonZero + " non zero values out of " + n.all + " sparseness is " + n.sparseness());

		// Filtering NaN and Inf values, in particular for molecules for prediction datasets

		if(selectionConfiguration.useAUTO != null && selectionConfiguration.useAUTO) { // filtering by autoencoder
			if(selectionConfiguration.autoEncoder == null)
				selectionConfiguration.autoEncoder	= new CompressedObject<Object>(createFilter(dtDescriptors));
			applyFilter(dtDescriptors,selectionConfiguration.autoEncoder);
		}

		if (dtDescriptors.getColumnsSize() != 0)  // we have descriptors...
			dtDescriptors.filterNaN();


		return new WorkflowNodeData(dtDescriptors).addPort(new DataTable(selectionConfiguration).setId("selection-configuration"));
	}



	private byte[] createFilter(DataTable dtDescriptors) throws IOException, InterruptedException {
		init();
		BufferedWriter writer = getAliasedBufferedWriter(DATAFILE);
		for (int mol = 0; mol < dtDescriptors.getRowsSize(); mol++){
			for(int i=0;i<dtDescriptors.getColumnsSize();i++){
				if(i==0)writer.append("  "+mol+",  ");
				String val = "" + dtDescriptors.getValue(mol, i);
				writer.append(val+",   ");
			}
			writer.append("\n");
		}
		writer.close();

		try {
		
		String program ="/Users/itetko/Desktop/mul/a.out";

		String[] commands = new String[] { program, DATAFILE, "--SCALE_MUL"};
		executeBinary(commands, DATAFILE+".txt", 0);

		commands[0] = "cp";commands[1]=DATAFILE+".txt"; commands[2]=DATAFILE+".tes";
		executeBinary(commands, DATAFILE+".tes", 0);

		commands[0] = program;commands[1]=DATAFILE+".txt"; commands[2]="--NET_AMISS";
		executeBinary(commands, "data.dbd", 0);

		commands[0] = program;commands[1]=DATAFILE+".txt"; commands[2]="--AAA";
		executeBinary(commands, DATAFILE+".wgt", 0);

		commands[0] = program;commands[1]=DATAFILE+".txt"; commands[2]="--DESCALE_AUTO";
		executeBinary(commands, "reconstruct.txt", 0);

		commands[0] = program;commands[1]="reconstruct.txt"; commands[2]="--BB_2_I";
		executeBinary(commands, "reconstruct.txt.txt", 0);

		commands = new String[] { "tar", "-cf", "model.tar",  "data.wgt", "data.dbd", "statts.txt"};
		executeBinary(commands, "model.tar", 0);
		}catch(Exception ee) {
			ee.printStackTrace();
		}

		return FileUtils.getFileAsBytes(getAliasedFileName("model.tar"));
	}


	private void applyFilter(DataTable dtDescriptors, CompressedObject<Object> autoEncoder) throws IOException, InterruptedException {
		String program ="/Users/itetko/Desktop/mul/a.out";

		String[] commands = new String[] { program, DATAFILE, "--SCALE_TES"};
		executeBinary(commands, DATAFILE+".txt", 0);

		commands[0] = "cp";commands[1]=DATAFILE+".txt"; commands[2]=DATAFILE+".val";
		executeBinary(commands, DATAFILE+".val", 0);

		commands[0] = program;commands[1]=DATAFILE+".val"; commands[2]="--AAA_TEST";
		executeBinary(commands, DATAFILE+".val", 0);

		
	}


	/**
	 * Trying to guess training set size; unsupervised descriptor selection is done using it only
	 * @param dtDescriptors
	 * @return
	 */

	int guessTrainingSet(DataTable dtDescriptors){
		int trainingSetSize = 0;
		dtDescriptors.reset();

		while(dtDescriptors.nextRow() && 
				!Boolean.TRUE.equals(dtDescriptors.getCurrentRow().getAttachment(QSPRConstants.VALIDATION)))
			trainingSetSize++;

		logger.info("Approximate training set size is "+trainingSetSize + " full set size is "+dtDescriptors.getRowsSize());

		return trainingSetSize;
	}


	/**
	 *  Perform actual filtering of descriptors using UFS or not
	 * @param dtDescriptors
	 * @param maxval
	 * @throws InterruptedException 
	 * @throws IOException 
	 */

	private void filterIt(SelectionConfiguration selectionConfiguration, DataTable dtDescriptors) throws IOException, InterruptedException{

		int trainingSetSize = guessTrainingSet(dtDescriptors);

		Set<String> columnsToDelete = new HashSet<String>();

		@SuppressWarnings("unchecked")
		Map<Float,Integer> uniques[]= new HashMap[dtDescriptors.getColumnsSize()];
		for(int i=0;i<uniques.length;i++)
			uniques[i] = new HashMap<Float,Integer>(); 
		if(selectionConfiguration.numDifferentValues > 0 || selectionConfiguration.maximumValueThreshold != Integer.MAX_VALUE ){

			setStatus("Fast filtering of columns with 0 and maximum values");  // pre-filtering by removing columns with all 0 and to large values to speed-up UFS

			for(int row=0;row<dtDescriptors.getRowsSize();row++){ 
				float vals[] = ((CompactDataRow)dtDescriptors.getRow(row)).toArray();
				for(int col=0;col<vals.length;col++){
					if(row<trainingSetSize && uniques[col].size() <= selectionConfiguration.maximumValueThreshold) {
						if(!uniques[col].containsKey(vals[col]))uniques[col].put(vals[col],1);
						else
							uniques[col].put(vals[col],1+uniques[col].get(vals[col]));
					}
					if(Math.abs(vals[col])>selectionConfiguration.maximumValueThreshold || !Float.isFinite(vals[col]) ) // also excluding too large values for the test set; later such errors will be marked as errors
						columnsToDelete.add(dtDescriptors.getColumn(col));
				}
			}
			for(int col=0;col<uniques.length;col++){
				if(uniques[col].size()<=selectionConfiguration.numDifferentValues){
					int sum =0;
					for(Float f:uniques[col].keySet())
						if(uniques[col].get(f)<= (1+uniques.length/2))sum += uniques[col].get(f); // only small counts are considered
					if(sum <= selectionConfiguration.numDifferentValues)
						columnsToDelete.add(dtDescriptors.getColumn(col));
				}
			}

			setStatus("Zero/Maxvalue filtering removed " + columnsToDelete.size()+ " columns."); 

			if (dtDescriptors.getColumnsSize() == columnsToDelete.size())
				throw new UserFriendlyException("All descriptors have identical/constant values and have been filtered out! Model calculation is impossible for " + dtDescriptors.getRowsSize() + " rows");

		}

		if(selectionConfiguration.useAUTO != null && selectionConfiguration.useAUTO) {

		}
		else
			if( (selectionConfiguration.useUFS || (
					selectionConfiguration.correlationThreshold >0 && 
					selectionConfiguration.correlationThreshold < 1)) && (dtDescriptors.getColumnsSize()  - columnsToDelete.size()) > 1){ 

				init();

				createCfgFile(selectionConfiguration,dtDescriptors.getColumnsSize()-columnsToDelete.size());
				saveOnlyTrainingSetDescripors(dtDescriptors,trainingSetSize,columnsToDelete); // 

				String[] commands = new String[] { getExeFile(), CFG, DATAFILE };
				executeBinary(commands, OCHEM, 0);

				Set<Integer> ufsFiltered=readColumns(OCHEM);

				cleanup();

				for(int i=0,n=0;i<dtDescriptors.getColumnsSize();i++)
					if(!columnsToDelete.contains(dtDescriptors.getColumn(i))){ // skipping already filtered columns
						if(ufsFiltered.contains(n))
							columnsToDelete.add(dtDescriptors.getColumn(i));
						n++;
					}

				setStatus("UFS has finished -- removing "+columnsToDelete.size()+" columns");
			}
			else
				setStatus("UFS was not used -- removing "+columnsToDelete.size()+" columns");

		setStatus("Selection filtering removed " + columnsToDelete.size()+ " columns."); 

		if(columnsToDelete.size()>0)
			dtDescriptors.deleteByList(columnsToDelete);

		if (dtDescriptors.getColumnsSize() == 0)
			throw new UserFriendlyException("All descriptors have been filtered out! Model calculation is impossible.");

		selectionConfiguration.setNoFiltering();

	}

	/*
	 * Some descriptors may fail on the stage of application
	 * Molecules with such  descriptors (not all molecules!) should be flagged as errors
	 */

	private void filterByMaximumValue(DataTable dtDescriptors, double maxval) {

		dtDescriptors.reset();

		while(dtDescriptors.nextRow()){
			AbstractDataRow r=dtDescriptors.getCurrentRow();
			for(int i=0;i<r.size();i++){
				if(Math.abs((Double)r.getValue(i))>=maxval){
					r.setError("Descriptor maxvalue error: "+dtDescriptors.getColumn(i)+"  failed, value="+r.getValue(i));
				}	

			}
		}

	}


	private Set<Integer> readColumns(String datafile) throws IOException
	{
		Set<Integer> columns=new HashSet<Integer>();

		BufferedReader file = getAliasedBufferedReader(datafile);
		String line;
		while ((line = file.readLine()) != null)
			if(line.length()>0)
				columns.add(Integer.valueOf(line).intValue());
		file.close();

		return columns;
	}

	private void createCfgFile(SelectionConfiguration selectionConfiguration,int variables) throws IOException{
		logger.info("Writing cfg file");
		String parameters = "FILE=data\nNAMES=1\nVARIABLE_NAMES=1\n"+
				"\nUFS="+(selectionConfiguration.useUFS?1:0)+
				"\nNONZERO="+selectionConfiguration.numDifferentValues+
				"\nCORRELATION="+selectionConfiguration.correlationThreshold+
				"\nSTD="+selectionConfiguration.stdThreshold+
				"\nMAXVALUE="+selectionConfiguration.maximumValueThreshold+
				"\nINPUTS="+variables;
		FileUtils.saveStringToFile(parameters+"\nSTOP\nEND", getAliasedFileName(CFG));
	}

	/**
	 *  Saves only descriptors (the qualitative options are supposed not to come to the filtering) -- to be checked!!!
	 * @param dtDescriptors
	 * @param molecules
	 * @throws IOException
	 */
	private void saveOnlyTrainingSetDescripors(DataTable dtDescriptors,
			int molecules, Set<String> skip) throws IOException {

		// checking that we do not have qualitative descriptors
		for (int col = 0; col < dtDescriptors.getColumnsSize(); col++)
			if (dtDescriptors.getColumnAttachment(dtDescriptors.getColumn(col), QSPRConstants.IS_CONDITION_COLUMN) != null) throw new IOException("Datatable in selection server contains also conditions -- but it should not!");

		List<Integer> rows = new ArrayList<Integer>();
		for (int i=0; i<molecules; i++) 
			rows.add(i);

		int saveMolecules = dtDescriptors.getColumnsSize() -skip.size() ==0? 0: 1024*1024*1024/(dtDescriptors.getColumnsSize() -skip.size()); // limitation for UFS filtering is 8GB (4GB for data and 4GB for UFS) or 1Gb of values 

		if(saveMolecules < molecules){
			setStatus("limiting number of molecules for UFS to "+saveMolecules);
			Collections.shuffle(rows, new Random(QSPRConstants.SEED));  // first we shuffle all data to avoid their biasing 
		}else
			saveMolecules = molecules;

		BufferedWriter writer = getAliasedBufferedWriter(DATAFILE);

		// write names
		writer.append(""+saveMolecules);
		writer.append("\nmols");

		for(int i=0;i<dtDescriptors.getColumnsSize();i++)
			if(!skip.contains(dtDescriptors.getColumn(i)))
				writer.append("\tc" + i);
		writer.append("\n");

		// write descriptors
		dtDescriptors.reset();
		for (int i = 0; i < saveMolecules; i++) {
			dtDescriptors.setCurrentRow(rows.get(i));
			writer.append("mol_" + i);

			float vals[]= ((CompactDataRow)dtDescriptors.getCurrentRow()).toArray();

			for(int col=0; col<dtDescriptors.getColumnsSize();col++)
				if(!skip.contains(dtDescriptors.getColumn(col))){
					writer.append("\t"+(col>=vals.length?0:NumericalValueStandardizer.getSignificantDigits(vals[col]))); // conversion is just to decrease the file size ...
				}
			writer.append("\n");
		}
		writer.append("\n");

		writer.flush();
		writer.close();
	}


	/**
	 * Required to implement non-existing method
	 */
	@Override
	protected DataTable calculateExecutableTask(DataTable dtInpout, Serializable configuration)
			throws Exception {
		return null;
	}

}
