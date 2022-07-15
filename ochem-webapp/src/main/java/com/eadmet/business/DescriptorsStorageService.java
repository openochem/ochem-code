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

package com.eadmet.business;

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.StringReader;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;

import qspr.Globals;
import qspr.dao.Repository;
import qspr.entities.Mapping2;
import qspr.entities.Molecule;
import qspr.export.ExportableMolecule;
import qspr.export.ExportableSetConfiguration;
import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.metaserver.configurations.DescriptorsEmptyConfiguration;
import qspr.workflow.utils.SDFProcessor;
import qspr.util.ExportThread;
import qspr.util.MoleculePeer;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.descriptorcache.CacheEntry;
import com.eadmet.descriptorcache.DescriptorConfigEntry;
import com.eadmet.descriptorcache.DescriptorsCache;
import com.eadmet.descriptorcache.DescriptorsRepository;
import com.eadmet.descriptorcache.DescriptorsRepositoryFactory;
import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.parsers.SimpleParser;

public class DescriptorsStorageService 
{
	static final int ENTRIES_TO_STORE = 1000;

	public class CacheUploadResult 
	{
		public Boolean externalIDsUpload;
		public Integer count;
	}


	private List<DescriptorConfigEntry> getSortedConfigs(String user, boolean all) throws Exception
	{
		if (Globals.isGuestUser())
			throw new UserFriendlyException("Guest users cannot work with the descriptors storage. Please, register and login to get this functionality.");

		DescriptorsRepository repository = DescriptorsRepositoryFactory.getRepository();

		List<DescriptorConfigEntry> configs = all ? repository.getAllConfigurations() : repository.getConfigurations(user);

		Collections.sort(configs, comparator);
		return configs;	
	}

	public List<DescriptorConfigEntry> getPublicConfigs() throws Exception
	{
		return getSortedConfigs(null, true);
	}

	public List<DescriptorConfigEntry> getPrivateConfigs() throws Exception
	{
		List<DescriptorConfigEntry> my = getSortedConfigs(Globals.userSession().user.login, false);
		my.addAll(getSortedConfigs(QSPRConstants.PUBLISHER, false));
		return my;
	}

	public void deleteCache(String configId, String user) throws Exception
	{
		DescriptorsRepository repository = DescriptorsRepositoryFactory.getRepository();
		DescriptorConfigEntry config = new DescriptorConfigEntry();
		config.objectID = configId;

		if ((user != null && !user.equals("")) && !user.equals(Globals.userSession().user.login))
			if (!Globals.isSuperUser())
				throw new UserFriendlyException("You do not have sufficient privileges for this operation.");

		if (user == null || user.equals(""))
		{
			if (!Globals.isSuperUser())
				throw new UserFriendlyException("Only super-users can delete public cache.");
			user = null;
		}

		config.setUser(user);
		repository.clearCache(config);
	}

	public CacheUploadResult descriptorsCacheUpload(File f, String descriptorType, String descriptorConfigXML) throws Exception
	{
		if (f == null)
			throw new UserFriendlyException("No file has been provided!");

		FileInputStream fis = new FileInputStream(f);
		CacheUploadResult result = descriptorsCacheUpload(fis, f.getName(), descriptorType, descriptorConfigXML);
		fis.close();
		return result;
	}

	public CacheUploadResult descriptorsCacheUpload(InputStream is, String fileName, String descriptorType, String descriptorConfigXML) throws Exception
	{
		if (Globals.isGuestUser())
			throw new UserFriendlyException("Guest users cannot work with the descriptors storage. Please, register and login to get this functionality.");

		if (descriptorType == null || descriptorType.equals(""))
			throw new UserFriendlyException("Descriptor type not specified!");

		DescriptorsRepository repository = DescriptorsRepositoryFactory.getRepository();

		SimpleParser parser = SimpleParser.getParser(fileName).setSource(is);
		parser.reset();
		List<String> header = parser.sheetColumns.get(0);

		String firstColumn = header.get(0).toUpperCase();
		if (!QSPRConstants.MOLECULE.equalsIgnoreCase(firstColumn) && !QSPRConstants.MOLECULEID.equalsIgnoreCase(firstColumn) 
				&& !QSPRConstants.SMILES_FORMAT.equalsIgnoreCase(firstColumn) && !QSPRConstants.EXTERNALID.equalsIgnoreCase(firstColumn))
			throw new UserFriendlyException("First column should be '"+QSPRConstants.MOLECULE+"' or '" + QSPRConstants.EXTERNALID +"' but it is " + firstColumn+
					". Change the order of columns in your file and repeat upload.");

		boolean externalIDsUpload = QSPRConstants.EXTERNALID.equalsIgnoreCase(firstColumn);

		List<String> descNames = new ArrayList<String>();
		for (int i = 1; i < header.size(); i++)
		{
			String name = header.get(i);
			if (name == null || name.trim().equals(""))
				break;
			descNames.add(name);
		}

		DescriptorsAbstractConfiguration conf = null;


		if (descriptorConfigXML == null || descriptorConfigXML.trim().length()==0)
			conf = new DescriptorsEmptyConfiguration();
		else
			try
		{
				Object oConf = Globals.jaxbContext.createUnmarshaller().unmarshal(new StringReader(descriptorConfigXML));
				if (!(oConf instanceof DescriptorsAbstractConfiguration))
					throw new Exception();
				conf = (DescriptorsAbstractConfiguration) oConf;
		}
		catch (Exception e)
		{
			throw new UserFriendlyException("The entered configuration is not a valid configuration XML. If unsure, leave the configuration field empty.");
		}

		/*
		if(descriptorType.contains(DescriptorsConfiguration.ExternalDescriptors)) {
			conf = new DescriptorsExternalConfiguration(descriptorType.substring(
					DescriptorsConfiguration.ExternalDescriptors.length()+1));
			descriptorType = DescriptorsConfiguration.ExternalDescriptors;
		}
*/
		DescriptorConfigEntry dConf = new DescriptorConfigEntry(conf, descriptorType);

		if(!Globals.userSession().user.isSuperUser())
			dConf.setUser(Globals.userSession().user.login);

		repository.saveConfig(dConf);

		Timestamp now = new Timestamp(Calendar.getInstance().getTimeInMillis());

		DescriptorsCache cache = new DescriptorsCache();

		int total = 0;

		boolean finish = false;
		while (!finish)
		{
			List<CacheEntry> entries = new ArrayList<CacheEntry>();
			int counter = 0;
			while (parser.hasNext() && (counter++ < ENTRIES_TO_STORE))
			{
				total++;
				List<String> values = parser.next();
				String molStr = values.get(0);
				if (molStr == null || molStr.trim().equals(""))
					break;

				CacheEntry ce = new CacheEntry();

				ce.config = dConf;
				ce.user = dConf.user;
				ce.dateCreated = now;

				if (!externalIDsUpload)
				{
					Molecule m = null;
					if (molStr.startsWith("M"))
					{
						Mapping2 mapping2 = Repository.molecule.getMapping2(Integer.valueOf(molStr.substring(1)));
						if (mapping2 != null)
							m = mapping2.getMolecule();
					}
					else
						m = MoleculePeer.getMolecule(molStr);
					if (m == null)
					{
						System.out.println("Error, molecule "+molStr+" not found");
						continue;
					}
					ce.setMD5(SDFProcessor.getMD5SDF(m.getData(),false));
					ce.mp2 = m.mapping2.id;
				}
				else
				{
					// For external molecules, store the external ID in the MD5 field
					ce.mp2 = null;
					ce.setMD5(molStr);
				}

				String[] anames = descNames.toArray(new String[0]);
				float[] avalues = new float[descNames.size()];

				for (int j = 0; j < descNames.size(); j++)
				{
					try 
					{
						Float val = Float.valueOf(values.get(1 + j));
						avalues[j] = (val == null) ? Float.NaN : val;
					} catch(Exception e)
					{
						avalues[j] = Float.NaN;
						System.out.println("Warning: can not get Float value from " + (values.get(1 + j) == null? "empty cell" : values.get(1 + j)) + " for molecule " + 
								(externalIDsUpload? "M" + ce.mp2 : molStr) + " for entry #" + total + ", uploadind NaN instead");
					}
				}

				ce.setNamesAndValues(anames, avalues, false);
				entries.add(ce);
			}
			cache.saveCacheEntries(dConf, entries);
			finish = !parser.hasNext();
			Globals.restartAllTransactions(true);
		}
		CacheUploadResult result = new CacheUploadResult();
		result.count = total;
		result.externalIDsUpload = externalIDsUpload;
		return result;
	}

	public ExportThread getExportThread(final String configId, final String userLogin, ExportableSetConfiguration conf, String format)
	{
		ExportThread eThread = new ExportThread(format, conf)
		{
			@Override
			public void generateData() throws Exception
			{
				DescriptorsRepository repository = DescriptorsRepositoryFactory.getRepository();
				DescriptorConfigEntry dConfig = repository.getDescriptorConfigById(configId);
				List<CacheEntry> entries = repository.getDescriptors(dConfig);

				DataTable dtDescriptors = new DataTable(true);
				HashMap<String, Integer> colsMap = new HashMap<String, Integer>();
				eData.setDescriptors(dtDescriptors);

				for (CacheEntry entry : entries)
				{
					dtDescriptors.addRow();
					if (entry.error != null)
						dtDescriptors.getCurrentRow().setError(entry.error);
					else
					{
						String[] names = entry.getNames();
						float[] values = entry.getValues();

						for (int i = 0; i < names.length; i++)
						{
							Integer colNum = colsMap.get(names[i]);
							if (colNum == null)
							{
								dtDescriptors.addColumn(names[i]);
								colsMap.put(names[i], colNum = dtDescriptors.getColumnsSize() - 1);
							}
							dtDescriptors.setValue(colNum, values[i]);
						}
					}

					ExportableMolecule eMol = new ExportableMolecule();
					eData.addMolecule(eMol);
					if (entry.mp2 != null && Repository.molecule.getMapping2(entry.mp2) != null)
						eMol.setMolecule(Repository.molecule.getMapping2(entry.mp2).getMolecule());
					else
					{
						eMol.inhouseRecordId = entry.moleculeMD5;
					}
					eMol.descriptors = dtDescriptors.getCurrentRow();

				}

				setFileName(dConfig.type);
			}

		};
		return eThread;
	}

	Comparator<DescriptorConfigEntry> comparator = new Comparator<DescriptorConfigEntry>(){

		@Override
		public int compare(DescriptorConfigEntry arg0, DescriptorConfigEntry arg1)
		{
			return arg0.toString().compareTo(arg1.toString());
		}
	};

}
