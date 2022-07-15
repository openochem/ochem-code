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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.hibernate.type.DoubleType;
import org.hibernate.type.LongType;

import qspr.Globals;
import qspr.dao.Various;
import qspr.workflow.utils.QSPRConstants
;
import qspr.modelling.PointStatistics;
import qspr.modelling.SetStatistics;

public class PredictionSpaceKnnCalculator 
{
	public static List<PredictionNeighbor> getNeighbours(PointStatistics point, SetStatistics trainingSet, String type)
	{
		return getNeighbours(point.ensemblePredictions, trainingSet, type);
	}

	public static List<PredictionNeighbor> getNeighbours(float[] predictionVector, SetStatistics trainingSet, String type)
	{
		long timer = System.nanoTime();
		List<PredictionNeighbor> list = new ArrayList<PredictionNeighbor>();

		for (PointStatistics nps : trainingSet.points)
		{
			PredictionNeighbor de = PredictionNeighbor.getNeighbour(predictionVector, nps, type);

			if (de != null)
				list.add(de);
		}
		Collections.sort(list);
		System.out.println((System.nanoTime() - timer)/1000000+"ms to find kNN for point");
		return list;
	}

	@SuppressWarnings("unchecked")
	public static List<PredictionNeighbor> getStructuralNeighbours(SetStatistics trainingSet, String molecule, Long propertyId) throws IOException
	{
		String smiles = Various.molecule.convertToFormat(molecule, QSPRConstants.SMILES_FORMAT);
		List<PredictionNeighbor> list = new ArrayList<PredictionNeighbor>();
		long timer = System.nanoTime();
		List<Object[]> rows = propertyId == null ? Globals.session().createSQLQuery("select similarity('" + smiles + "', simscreen) s, exp_property_id epid from ExperimentalProperty " + 
				"join BasketEntry be using (exp_property_id) join Molecule using (molecule_id) join Mapping2 mp2 using (mapping2_id) " + 
				"join StructureQuery using (mapping2_id) where be.basket_id=:tsId order by s desc limit 11")
				.addScalar("s", DoubleType.INSTANCE)
				.addScalar("epid", LongType.INSTANCE)
				.setParameter("tsId", trainingSet.basketId)
				.list() 
				: Globals.session().createSQLQuery("select similarity('" + smiles + "', simscreen) s, exp_property_id epid from ExperimentalProperty epr " + 
						"join BasketEntry be using (exp_property_id) join Molecule using (molecule_id) join Mapping2 mp2 using (mapping2_id) " + 
						"join StructureQuery using (mapping2_id) where be.basket_id=:tsId and epr.property_id=:prId order by s desc limit 10")
						.addScalar("s", DoubleType.INSTANCE)
						.addScalar("epid", LongType.INSTANCE)
						.setParameter("tsId", trainingSet.basketId)
						.setParameter("prId", propertyId)
						.list();


				for (Object[] row : rows) {
					for (PointStatistics ps : trainingSet.points) {
						if ((Long) row[1] == ps.id)
						{
							list.add(new PredictionNeighbor(ps, ((Double) row[0])/100));
							continue;
						}
					}
					//throw new UserFriendlyException("Shooher! Coul not find point with id " + (Long) row[1]);
				}

				System.out.println((System.nanoTime() - timer)/1000000+"ms to find structural neighbors for point");
				return list;
	}
}