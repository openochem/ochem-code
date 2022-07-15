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
import java.io.PrintWriter;
import java.io.Serializable;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import com.eadmet.utils.OCHEMUtils;

import qspr.dao.Various;
import qspr.metaserver.configurations.StructureOptimisationConfiguration;
import qspr.metaserver.cs.util.NaiveConnectionManager;
import qspr.metaserver.cs.util.SavedMolecule;
import qspr.workflow.utils.QSPRConstants;
import qspr.workflow.utils.SDFProcessor;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;

public abstract class AbstractStructureConversionServer extends SimpleAbstractServer {

	NaiveConnectionManager cml = new NaiveConnectionManager();

	abstract void init() throws IOException, InterruptedException;
	abstract void cleanup() throws IOException, InterruptedException;
	abstract String convertedStructure(String molecule, StructureOptimisationConfiguration config) throws Exception;

	final static int DEFAULT_BATCH = 1000;

	String method;

	@Override
	public void setParam(String name, String value)
	{
		name = name.toUpperCase();
		super.setParam(name, value);
		if(cml != null) {
			switch(name) {
			case "DBURL": cml.dbUrl = value; break;
			case "DBUSERNAME": cml.dbUsername = value; break;
			case "DBPASSWORD": cml.dbPassword = value; break;
			case "BYPASS": cml = null;
			}
		}
	}

	@Override
	DataTable calculateGeneralTask(DataTable dtMolecules, Serializable config) throws Exception
	{
		if( Various.molecule == null)
			Various.molecule = Various.getDefaultCheminfImpl();

		StructureOptimisationConfiguration configuration=(StructureOptimisationConfiguration)config;

		NaiveConnectionManager cmlocal = cml;

		if(cmlocal == null || cmlocal.dbUrl == null || configuration.bypassCacheConfiguration()) {
			cmlocal = null; // only for this round!
			setStatus("caching is not used.");
		}
		else
			Class.forName(QSPRConstants.DEFAULTDATABASEDRIVER).newInstance();

		init();

		DataTable dtResults = new DataTable();
		dtResults.id = "optimized-structures";
		dtResults.addColumn(QSPRConstants.SDF_COLUMN);

		method = configuration.getTaskType();

		Map<Integer,SavedMolecule> mol= new HashMap<Integer,SavedMolecule>();
		Map<String,Integer> stored = new HashMap<String,Integer>();

		int size = dtMolecules.getRowsSize();
		int batchSize = repostSize == 0 ? DEFAULT_BATCH : repostSize/10;

		int skipAll = 0, cachedAll = 0, errors = 0;

		while(dtResults.getRowsSize() < size) { // as long as we have not processed all records!

			int skip = 0;
			mol.clear();

			for(int i = dtResults.getRowsSize(); i< size && mol.size() < batchSize; i++)
				try{
					dtResults.addRow(); // add new result

					AbstractDataRow data = dtMolecules.getRow(i);
					String externalID = (String)data.getAttachment(QSPRConstants.EXTERNAL_ID);
					Integer mapping2 = (Integer)data.getAttachment(QSPRConstants.MOLECULE_ID_STEREOCHEM);

					SavedMolecule molecule = null;
					molecule = new SavedMolecule((String)data.getValue(0), externalID,  mapping2 == null? 0 : mapping2);

					if(molecule.originalStructure == null) {  //empty structure!
						if(molecule.externId == null) // can be allowed if descriptors are detected by the ID; thus just skipping these structures unless they are errors
							dtResults.getCurrentRow().setError("Empty structure was received.");
						skip++;
						continue;
					}

					if(stored.containsKey(molecule.md5)) {  // found from previous analyzes; but NOT in this run (to be cached first)!
						int storedRow = stored.get(molecule.md5);
						dtResults.getCurrentRow().setValue(0, dtResults.getValue(storedRow, 0));
						if(dtResults.getRow(storedRow).isError())
							dtResults.getCurrentRow().setError(dtResults.getRow(storedRow).detailedStatus);
						skip++;
						continue;
					}
					mol.put(i, molecule);
				}catch(IOException e) {
					dtResults.getCurrentRow().setError(e.getLocalizedMessage());
					skip++;
					continue;
				}

			int cached = cmlocal != null? getFromCache(mol.values(), cmlocal, out) : 0;

			if(mol.size() != cached){ // new molecules are available!
				if(cmlocal != null) cmlocal.closeConnection(); // to avoid keeping connections for time when performing the conversion

				// long operation
				for(SavedMolecule molecule:mol.values())
					if (!molecule.cached)
						processMolecule(molecule, configuration); 

				if(cmlocal != null) saveToCache(mol.values(), cmlocal);  // connection is open again at this time
			}

			setStatus("Using " +  cached + " cached " +( mol.size() != cached? "and saving " + (mol.size() - cached) +" new molecules " :"")  + (skip > 0? "skip "+ skip:"") + 
					" starting with  mol n = " + dtResults.getRowsSize() + " out of " + size);

			skipAll += skip;
			cachedAll += cached;

			for(Integer row: mol.keySet()) {
				SavedMolecule molecule = mol.get(row);
				if(molecule.error != null) {
					dtResults.getRow(row).setError(molecule.error); // can be conversion with an error; let anyway
					errors++;
				}
				else
					dtResults.setValue(row, 0, molecule.optimizedStructure);
				stored.put(molecule.md5, row); // will be our new reference for this md5
			}
		}

		setStatus("Finished conversion: skipped " +  skipAll +  " cached " +  cachedAll + " and saved  " + (dtResults.getRowsSize() - cachedAll - skipAll) + " errors = " + errors + " out of " + size);

		if(cmlocal != null) cmlocal.closeConnection();
		cleanup();
		return dtResults;
	}


	void processMolecule(SavedMolecule molecule, StructureOptimisationConfiguration config) throws InterruptedException
	{
		int MAX_TRIALS = 3;
		for(int trial=0; trial<MAX_TRIALS;trial++)
			try
		{
				String sdf = convertedStructure(molecule.processedStructure, config); // converting molecule
				molecule.optimizedStructure = SDFProcessor.standartize(sdf);
				String diff = molecule.compareFormulas();
				if(diff != null)
					throw new IOException("structure of the molecule was changed: " + diff);
				break;
		} catch (Exception e)
		{
			molecule.error = supportedTaskType + ": " + e.getMessage();
			e.printStackTrace(out);
			Thread.sleep(1000);
		}
	}

	int getFromCache(Collection<SavedMolecule> mol, NaiveConnectionManager cm, PrintWriter out) throws IOException, SQLException, InterruptedException
	{
		if(mol.size() == 0)return 0;

		StringBuffer buf = new StringBuffer();

		for(SavedMolecule molecule: mol) 
			buf.append( (buf.length() > 0 ?",":"") + "'"+molecule.md5+"'");

		PreparedStatement ps = cm.getStatement("select optimized_structure, error, md5 from SavedMolecule where method='"+
				method + "' and md5 in (" + buf.toString() + ")", false);

		ResultSet rs = ps.executeQuery();

		int cached = 0;

		while(rs.next()) 
			for(SavedMolecule molecule: mol)
				if(rs.getString(3).equals(molecule.md5)) {
					molecule.optimizedStructure = rs.getBytes(1) == null? null : new String(OCHEMUtils.MySqlCompatibleUncompress(rs.getBytes(1)));
					molecule.error = rs.getString(2);
					molecule.cached =  molecule.error == null || !molecule.error.contains("timed out");
					if(molecule.cached)cached++;
				}

		rs.close();

		return cached;
	}

	void saveToCache(Collection<SavedMolecule> mol, NaiveConnectionManager cm) throws SQLException, InterruptedException, IOException
	{
		PreparedStatement ps = cm.getStatement("insert ignore into SavedMolecule(N, md5, original_structure, processed_structure, optimized_structure, error, method, mapping2, time) "
				+ "values(0, ?, ?, ?, ?, ?, ?, ?, NOW())", true);
		if(ps == null) return;

		for(SavedMolecule molecule: mol) 
		{
			ps.setString(1, molecule.md5);
			ps.setBytes(2, OCHEMUtils.MySqlCompatibleCompress(molecule.originalStructure));
			ps.setBytes(3, OCHEMUtils.MySqlCompatibleCompress(molecule.processedStructure));
			ps.setBytes(4, OCHEMUtils.MySqlCompatibleCompress(molecule.optimizedStructure));
			ps.setString(5, molecule.error);
			ps.setString(6, method);
			ps.setInt(7, molecule.mapping2);
			ps.execute();
		}
	}



}


