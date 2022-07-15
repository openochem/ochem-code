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

import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;

import qspr.Globals;
import qspr.ThreadScope;
import qspr.entities.Basket;
import qspr.entities.PendingTask;
import qspr.entities.Property;
import qspr.entities.Session;
import qspr.modelling.applier.ModelApplier;
import qspr.util.RWriter;
import qspr.util.WrapperThread;

import com.eadmet.mmpa.domain.MMPSubset;
import com.eadmet.mmpa.domain.MMPTransformation;
import com.eadmet.mmpa.domain.ThinPair;
import com.eadmet.utils.MAILERConstants;

public class MMPSignificanceService
{
	public static List<MMPTransformation>[] cachedValue;
	public void main() throws IOException, ClassNotFoundException, Exception
	{
		Basket b = Basket.getBasket("Thesis_CYP3A4_Full", Session.getLastSession(MAILERConstants.TESTER));
		Property p = Property.getByName("CYP450 modulation");
		PendingTask[] pt = new PendingTask[]{PendingTask.getById(393943L)};

		List<MMPTransformation>[] sigData = getSignificanceData(b, p, pt);

		List<MMPTransformation> transformationsMeasured = sigData[0];
		List<MMPTransformation> transformationsCombined = sigData[1];

		RWriter writer = new RWriter(new FileOutputStream("/l/ndata/mmp-sig-study"));
		writeSignificanceChart(transformationsMeasured, p.isQualitative(), writer, "mmp.measured");
		writeSignificanceChart(transformationsCombined, p.isQualitative(), writer, "mmp.combined");
		writer.close();
	}

	@SuppressWarnings("unchecked")
	public static List<MMPTransformation>[] getSignificanceData(Basket basket, Property property, PendingTask[] pTasks) throws IOException, ClassNotFoundException, Exception {

		if (cachedValue != null)
			return cachedValue;
		MMPQueryService mmpQueryService = new MMPQueryService();

		MMPSubset measuredMMPs = mmpQueryService.getPairsByBasket(new MMPQuery().setBasket(basket.id).setSimilarity((short)50));
		measuredMMPs.molValues.put(property.id, MMPStatsService.getMolValues(basket.id, property.id, null));
		List<MMPTransformation> transformationsMeasured = MMPStatsService.getTransformationsStats(new MMPStatsRequest(measuredMMPs, 1.0).setMinPairs(0).setProperty(property));

		MMPSubset combinedMMPs = null;
		List<MMPTransformation> transformationsCombined = null;
		if (pTasks != null)
		{
			combinedMMPs = measuredMMPs.clone();

			for (PendingTask pTask : pTasks)
			{
				ModelApplier applier = new ModelApplier(pTask);
				applier.useCache = false;
				MMPSubset predictedMMPs = mmpQueryService.getPairsFromPrediction(new MMPQuery(applier, true).setSimilarity((short)50));
				mergeMMPs(predictedMMPs, combinedMMPs);
			}

			Globals.restartAllTransactions(true);
			transformationsCombined = MMPStatsService.getTransformationsStats(new MMPStatsRequest(combinedMMPs, 1.0).setMinPairs(0).setProperty(property));
		}

		return cachedValue = new List[]{transformationsMeasured, transformationsCombined};
	}

	public void writeSignificanceChart(List<MMPTransformation> transformations, boolean classification, RWriter writer, String var) {
		Collections.sort(transformations, new Comparator<MMPTransformation>()
		{

			@Override
			public int compare(MMPTransformation arg0, MMPTransformation arg1)
			{
				return arg0.id > arg1.id ? 1 : -1;
			}
		});

		for (MMPTransformation tr : transformations)
		{
			double sig = tr.statistics.pValue > 0 ? -Math.log(tr.statistics.pValue) / Math.log(10) : 16;
			double delta = 0;

			if (classification)
			{
				double np = tr.statistics.nNP * 1.0D / (tr.statistics.nNP + tr.statistics.nNN);
				double pn = tr.statistics.nPN * 1.0D / (tr.statistics.nPN + tr.statistics.nPP);

				if (Double.isNaN(np))
					delta = pn;
				else if (Double.isNaN(pn))
					delta = np;
				else
					delta = Math.max(np, pn);
			} else
				delta = Math.abs(tr.statistics.deltaMean);

			writer.addValue(var + ".sig", sig);
			writer.addValue(var + ".delta", delta);
			writer.addValue(var + ".count", Math.abs(tr.statistics.pairsCount));
			writer.addValue(var + ".tr", Math.abs(tr.id));
		}

		writer.write();
	}

	public static void mergeMMPs(MMPSubset source, MMPSubset destination) {
		Set<Long> newTransformationIds = source.getTransformationIds();
		for (Long tr : destination.getTransformationIds())
			if (newTransformationIds.contains(tr))
				for (ThinPair tp : source.getThinPairs(tr))
					destination.addPair(tr, tp.id, tp.mol1Id, tp.mol2Id);

		for (Long propId : destination.molValues.keySet())
			for (Integer mp2 : source.molValues.get(propId).keySet())
				destination.molValues.put(propId, mp2, source.molValues.get(propId).get(mp2));
	}

	public void printDiff(MMPSubset newm, MMPSubset oldm) {
		for (Long tr : oldm.getTransformationIds()) {
			System.out.println("" + tr + ": " + newm.getThinPairs(tr).size() + ", " + oldm.getThinPairs(tr).size() + ", " + (newm.getThinPairs(tr).size() - oldm.getThinPairs(tr).size()));
		}
	}


	public static void main(String[] args)
	{
		new WrapperThread()
		{

			@Override
			public void wrapped() throws Exception
			{
				ThreadScope.get().userSession = Session.getFirstSession(MAILERConstants.ADMIN);
				new MMPSignificanceService().main();

			}
		}.run();
	}
}
