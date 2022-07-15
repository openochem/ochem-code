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

package com.eadmet.mmpa;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import com.eadmet.mmpa.domain.MMPTransformation;
import com.eadmet.utils.Counter;
import com.eadmet.utils.Grouping;
/**
 * Experimental. Delete if the feature does not take off.
 */
public class MMPFragmentService
{

	public static List<long[]> getFragmentGraph(TransformationFilter filter, Integer cluster, Long fragmentId) throws Exception {
		return getFragmentGraph(new MMPQueryService().getTransformations(filter), cluster, fragmentId);
	}

	public void pruneGraph(List<long[]> edges) {
		Counter<Long> counter = new Counter<Long>();
		for (long[] ls : edges)
		{
			counter.add(ls[0]);
			counter.add(ls[1]);
		}
	}

	public static List<long[]> getFragmentGraph(long propertyId, Integer cluster, Long fragmentId) throws Exception
	{
		List<MMPTransformation> transformations;

		TransformationFilter filter = new TransformationFilter();
		filter.getEffects(0L, propertyId)[0] = true;
		transformations = new MMPQueryService().getTransformations(filter);
		//inverseWhomNecessary(transformations);
		//		No idea if it will work like this.

		//		filter = new TransformationFilter();
		//		filter.getEffects(0L, propertyId)[2] = true;
		//		List<MMPTransformation> transformationsOut = new MMPQueryService().getCriteriaFor(filter).list();
		//		inverse(transformationsOut);
		//		transformations.addAll(transformationsOut);

		return getFragmentGraph(transformations, cluster, fragmentId);
	}

	public static List<long[]> getFragmentGraph(List<MMPTransformation> transformations, Integer cluster, Long fragmentId) throws Exception {
		List<long[]> edges = new ArrayList<long[]>();
		Map<Long, Long> mins = new HashMap<Long, Long>();

		for (MMPTransformation mmpTransformation : transformations)
		{
			mins.put(mmpTransformation.getFrag1Id(), Math.min(mmpTransformation.getFrag1Id(), mmpTransformation.getFrag2Id()));
			mins.put(mmpTransformation.getFrag2Id(), Math.min(mmpTransformation.getFrag1Id(), mmpTransformation.getFrag2Id()));
		}

		for (int k = 0; k < 10; k++)
		{
			for (MMPTransformation mmpTransformation : transformations)
			{
				long val = Math.min(mins.get(mmpTransformation.getFrag1Id()), mins.get(mmpTransformation.getFrag2Id()));
				mins.put(mmpTransformation.getFrag1Id(), val);
				mins.put(mmpTransformation.getFrag2Id(), val);
			}
		}

		Set<Long> uqMinsSet = new HashSet<Long>();
		uqMinsSet.addAll(mins.values());
		List<Long> uqMins = new ArrayList<Long>();
		uqMins.addAll(uqMinsSet);

		reduce(transformations);

		Set<Long> fragFilter = null;
		if (fragmentId != null)
			fragFilter = identifyChain(transformations, fragmentId);

		for (MMPTransformation transformation : transformations)
			if (cluster == null || mins.get(transformation.getFrag1Id()).equals(uqMins.get(cluster)))
				if (fragFilter == null || (fragFilter.contains(transformation.getFrag1Id()) && fragFilter.contains(transformation.getFrag2Id())))
					edges.add(new long[]{transformation.getFrag1Id(), transformation.getFrag2Id()});

		return edges;	
	}

	private static Set<Long> identifyChain(List<MMPTransformation> transformations, long fragId) {
		Grouping<Long, Long> outgoing = new Grouping<Long, Long>();
		Grouping<Long, Long> incoming = new Grouping<Long, Long>();

		for (MMPTransformation tr : transformations)
		{
			outgoing.add(tr.getFrag1Id(), tr.getFrag2Id());
			incoming.add(tr.getFrag2Id(), tr.getFrag1Id());
		}

		reduce(transformations);

		Set<Long> chainFrags = new HashSet<Long>();

		identifyChain(fragId, chainFrags, outgoing, 0);
		identifyChain(fragId, chainFrags, incoming, 0);

		return chainFrags;
	}

	private static void reduce(List<MMPTransformation> transformations) {
		int cnt = 0;
		Grouping<Long, Long> outgoing = new Grouping<Long, Long>();
		Grouping<Long, Long> incoming = new Grouping<Long, Long>();

		for (MMPTransformation tr : transformations)
		{
			outgoing.add(tr.getFrag1Id(), tr.getFrag2Id());
			incoming.add(tr.getFrag2Id(), tr.getFrag1Id());
		}

		// Transitive reduction
		Iterator<MMPTransformation> i = transformations.iterator();
		while (i.hasNext())
		{
			MMPTransformation t = i.next();
			for (Long f : outgoing.get(t.getFrag1Id()))
				if (incoming.get(t.getFrag2Id()).contains(f))
				{
					cnt++;
					i.remove();
					break;
				}
		}

		logger.info("Removed " + cnt + " transfrmations by tranitive reduction");
	}

	private static void identifyChain(long fragId, Set<Long> fragList, Grouping<Long, Long> directions, int depth) {

		if (depth > 4)
			return;

		fragList.add(fragId);
		if (directions.get(fragId) != null)
			for (Long nextFrag : directions.get(fragId))
				identifyChain(nextFrag, fragList, directions, depth + 1);
	}

	private static final Logger logger = LogManager.getLogger(MMPFragmentService.class);
}

