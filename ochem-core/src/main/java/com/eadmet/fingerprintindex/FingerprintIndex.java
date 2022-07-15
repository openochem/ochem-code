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

package com.eadmet.fingerprintindex;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.Globals;
import qspr.dao.Various;

import com.eadmet.mmpa.domain.ThinPair;
import com.eadmet.utils.MemoryUtils;
import com.eadmet.utils.config.ConfigurableClass;
import com.eadmet.utils.config.ConfigurableProperty;

@ConfigurableClass(name = "fingerprintindex", comment = "Mapping2 structural fingerprint index for misc puproses")
@SuppressWarnings("unchecked")
public class FingerprintIndex 
{
	@ConfigurableProperty(name = "enabled", comment = "If enabled, the fingerprint index will be filled in database and in memory, and structural queries can be performed against it")
	public static boolean enabled = false;

	private static FingerprintIndex instance;

	private static final Logger logger = LogManager.getLogger(FingerprintIndex.class);


	long indexLoaded = 0;
	List<Integer> ids = new ArrayList<Integer>();
	List<byte[]> fps = new ArrayList<byte[]>();
	Map<Integer, byte[]> map = new HashMap<Integer, byte[]>();

	private FingerprintIndex()
	{

	}

	public List<Object[]> getMMPClusters(List<ThinPair> pairs, double similarity)
	{
		List<Integer> mappings = new ArrayList<Integer>();
		Set<Integer> mappingSet = new HashSet<Integer>();
		for (int i = 0; i < pairs.size(); i++) 
		{
			ThinPair pair = pairs.get(i);
			mappingSet.add(pair.mol1Id);
		}
		mappings.addAll(mappingSet);
		return getMappingClusters(mappings, similarity);
	}

	public List<Object[]> getMappingClusters(List<Integer> mappings, double threshold) //Naive quadratic-complexity tanimoto-based clustering, returns "centroids"
	{
		class Elem
		{
			Integer mp2;
			Set<Integer> neighbors = new HashSet<Integer>();
		}


		if (indexLoaded == 0)
			loadIndex();

		Map<Integer, Elem> elemMap = new HashMap<Integer, Elem>();

		for (Integer mp2 : mappings)
		{
			Elem e = new Elem();
			e.mp2 = mp2;
			elemMap.put(mp2, e);
		}

		int size = mappings.size();
		for (int i=0; i<size - 1; i++)
		{
			Integer mp2_1 = mappings.get(i);
			if (map.get(mp2_1) == null)
				continue;

			for (int j=i+1; j<size; j++)
			{
				Integer mp2_2 = mappings.get(j);
				if (map.get(mp2_2) == null)
					continue;

				double dissimilarity;
				try {
					dissimilarity = Various.molecule.getTanimoto(map.get(mp2_1), map.get(mp2_2));
				} catch (Exception e) {
					logger.error("Error while calculating dissimilarity.");
					e.printStackTrace();
					dissimilarity = 1;
				}
				if (dissimilarity < 1 - threshold)
				{
					elemMap.get(mp2_1).neighbors.add(mp2_2);
					elemMap.get(mp2_2).neighbors.add(mp2_1);
				}
			}
		}

		List<Elem> l = new ArrayList<Elem>();
		l.addAll(elemMap.values());

		List<Object[]> result = new ArrayList<Object[]>();
		Object[] unclustered = new Object[]{null, new ArrayList<Integer>()};
		while (l.size() > 0)
		{
			Collections.sort(l, new Comparator<Elem>(){
				@Override
				public int compare(Elem arg0, Elem arg1) 
				{
					Integer s0 = arg0.neighbors.size();
					Integer s1 = arg1.neighbors.size();
					return s1.compareTo(s0);
				}});

			Elem e = l.remove(0);
			List<Integer> maps = null; 
			if (e.neighbors.size() > 1)
				maps = new ArrayList<Integer>();
			else
				maps = ((List<Integer>)unclustered[1]);

			maps.addAll(e.neighbors);
			maps.add(e.mp2);

			if (e.neighbors.size() > 1)
				result.add(new Object[]{e.mp2, maps});

			int index = 0;
			while (index < l.size())
			{
				if (e.neighbors.contains(l.get(index).mp2))
					l.remove(index);
				else
				{
					l.get(index).neighbors.removeAll(e.neighbors);
					//					if (l.get(index).neighbors.size() == 0)
					//						l.remove(index);
					//					else
					index++;
				}
			}
		}
		result.add(unclustered);
		return result;
	}

	public List<Integer> queryIndex(String sdf, double similarity)
	{

		if (indexLoaded == 0)
			loadIndex();

		indexLoaded = System.nanoTime();

		long timer = System.nanoTime();
		List<Integer> matches = new ArrayList<Integer>();
		try 
		{
			

			for (int i=0; i<fps.size(); i++)
			{

				double dissimilarity = Various.molecule.getTanimoto(sdf, fps.get(i));

				if (dissimilarity < 1 - similarity)
					matches.add(ids.get(i));
			}
		} catch (Exception e) 
		{
			//e.printStackTrace();
		}
		logger.info("Query index with similairty " + similarity + ", found " + matches.size() + " mappings in " + (System.nanoTime() - timer) / 1E6 + "ms");
		return matches;
	}

	public void clearIndexIfOlder(int minutes)
	{
		if ((System.nanoTime() - indexLoaded) / 1E9 / 60 > minutes)
			clearIndex();
	}

	public synchronized void clearIndex()
	{
		if (indexLoaded > 0)
		{
			ids.clear();
			fps.clear();
			map.clear();
			indexLoaded = 0;
		}
	}

	public synchronized void loadIndex()
	{
		if (indexLoaded > 0 || !enabled)
			return;

		long clTimer = System.nanoTime();
		List<Object[]> pids =  Globals.session().createSQLQuery("select mapping2_id, fingerprint from Mapping2Fingerprints").list();

		for (Object[] objects : pids) 
		{
			Integer mp2 = (Integer)objects[0];
			byte[] fp = (byte[])objects[1];

			if (fp == null)
				continue;

			ids.add(mp2);
			fps.add(fp);
			map.put(mp2, fp);
		}
		indexLoaded = System.nanoTime();
		logger.info("Cache load ("+ids.size()+" entries) in " + (System.nanoTime() - clTimer)/1E9 + "s");
		logger.info(MemoryUtils.memorySummary());
	}

	public static FingerprintIndex getInstance()
	{
		if (instance != null)
			return instance;

		synchronized(FingerprintIndex.class)
		{
			if (instance != null)
				return instance;
			return instance = new FingerprintIndex();
		}
	}

	public void refreshIndex()
	{
		Set<Integer> idSet = new HashSet<Integer>();

		List<Integer> ids = Globals.session().createSQLQuery("select mapping2_id from Mapping2").list();
		idSet.addAll(ids);
		Globals.restartAllTransactions(false);

		List<Integer> pids =  Globals.session().createSQLQuery("select mapping2_id from Mapping2Fingerprints").list();
		idSet.removeAll(pids);
		Globals.restartAllTransactions(false);

		List<Integer> allIds = new ArrayList<Integer>();
		allIds.addAll(idSet);

		int size = allIds.size();
		long global = System.nanoTime();
		while (allIds.size() > 0)
		{
			long timer = System.nanoTime();
			List<Integer> subIds = allIds.subList(0, Math.min(allIds.size(), 1000));
			for (Integer mapping2id : subIds) 
			{
				try
				{
					List<byte[]> l =  Globals.session().createSQLQuery("select uncompress(molecule_data) from Molecule where mapping2_id = :mp2").setInteger("mp2", mapping2id).list();
					if (l.size() == 0)
						continue;
				
					byte[] data = Various.molecule.getFingerprint(new String(l.get(0)));
					Globals.session().createSQLQuery("insert ignore into Mapping2Fingerprints values (:mp2, :blob)").setInteger("mp2", mapping2id).setBinary("blob", data).executeUpdate();
				} catch (Exception e)
				{
					e.printStackTrace();
				}
			}
			subIds.clear();
			logger.info(size - allIds.size() + " out of " + size + " molecules processed; last batch: "+(System.nanoTime() - timer)/1E9+"s; total time: "+(System.nanoTime() - global)/1E9+"s");
			Globals.restartAllTransactions(true);
		}
	}
}
