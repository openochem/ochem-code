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

package qspr.modelling;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.entities.Basket;
import qspr.entities.BasketEntry;
import qspr.entities.ModelMapping;
import qspr.interfaces.Descriptable;
import qspr.workflow.utils.QSPRConstants
;
import qspr.workflow.datatypes.DataTable;

@XmlRootElement(name = "statistics")
public class ModelStatistics implements Serializable, Descriptable {
	private static transient final Logger logger = LogManager.getLogger(ModelStatistics.class);

	public static final long serialVersionUID = 1;

	public boolean containsPredictions;

	@XmlElement(name = "set")
	public List<SetStatistics> sets = new ArrayList<SetStatistics>();

	public transient List<Long> deletedRecordsIds = new ArrayList<Long>();

	/**
	 * Number of the validation set of interest (starting from 1)
	 */
	public String validationSetId;

	public String toString() {
		String res = "Training set: " + sets.get(0).toString();
		if (sets.size() > 0)
			res += "\nValidation set: " + sets.get(1).toString();

		return res;
	}

	public static ModelStatistics get(ModelMapping mm) {
		ModelStatistics ms = (ModelStatistics) mm.statisticsRecalculated.getObject();
		if (!ms.containsPredictions)
			ms = (ModelStatistics) mm.statisticsOriginal.getObject();
		return ms;
	}

	public Map<String, Object> getParameters() {
		Map<String, Object> parameters = new HashMap<String, Object>();

		Map<String, Object> tsParameters = sets.get(0).getParameters();
		for (String param : tsParameters.keySet())
			parameters.put(param + " training set", tsParameters.get(param));

		if (sets.size() > 1) {
			Map<String, Object> validationSetParameters = sets.get(1).getParameters();
			for (String param : validationSetParameters.keySet())
				parameters.put(param + " test set", validationSetParameters.get(param));
		}

		return parameters;
	}

	public void setPredictions(DataTable dtPredictions, ModelMapping mm, List<Set<Long>> moleculesInSets) throws Exception {
		dtPredictions.reset();
		for (int i=0; i< sets.size(); i++)
		{
			SetStatistics ss = sets.get(i);
			dtPredictions.reset();
			ss.setPredictions(dtPredictions, mm, moleculesInSets == null || moleculesInSets.size() <= i? null : moleculesInSets.get(i));
		}
		containsPredictions = true;
	}

	public static ModelStatistics getEmptyStatistics(ModelMapping mm) {
		List<String> setIds = new ArrayList<String>();
		setIds.add(QSPRConstants.TRAINING);
		for (int i = 0; i < mm.model.getValidationSets().size(); i++)
			setIds.add(QSPRConstants.VALIDATION + i);
		setIds.add(QSPRConstants.EXCLUDED);

		ModelStatistics ms = new ModelStatistics();
		for (String setId : setIds) {
			Basket set = mm.model.getFilteredSet(setId);
			if (set != null) {
				SetStatistics ss = new SetStatistics(set, mm);
				ss.setId(setId);
				ms.sets.add(ss);
			}
		}

		ms.containsPredictions = false;

		return ms;
	}

	public ModelStatistics recalculateStatistics(ModelMapping modelMapping) {
		return recalculateStatistics(modelMapping, null);
	}

	public ModelStatistics recalculateStatistics(ModelMapping modelMapping, PointSelector selector) {
		// try{return recalculateStatistics(modelMapping,0.1,true);}catch(Exception e){}
		for (SetStatistics set : sets)
			set.recalculateStatistics(modelMapping, selector);

		return this;
	}

	public void setBootstrapReplicas(int replicas) {
		for (SetStatistics set : sets)
			set.bootstrapReplicas = replicas;
	}

	@SuppressWarnings("unchecked")
	public void actualizeStatistics(ModelMapping modelMapping) {

		long time = Calendar.getInstance().getTimeInMillis();

		// Remove irrelevant validation sets
		if (!"all".equals(validationSetId)) {
			int num = validationSetId == null ? 0 : Integer.valueOf(validationSetId) - 1;
			if (modelMapping.model.getValidationSets().size() > num) {
				modelMapping.model.selectedValidationSet = modelMapping.model.getValidationSets().get(num);
				Long selectedBasketId = modelMapping.model.selectedValidationSet.id;
				Iterator<SetStatistics> iSets = sets.iterator();
				while (iSets.hasNext()) {
					SetStatistics ss = iSets.next();
					if (ss.setId.startsWith(QSPRConstants.VALIDATION))
						if (ss.basketId == null || ss.basketId.equals(selectedBasketId))
							ss.setId = QSPRConstants.VALIDATION;
						else
							iSets.remove();
				}
			}
		}

		if (deletedRecordsIds == null)
			deletedRecordsIds = new ArrayList<Long>();
		else
			deletedRecordsIds.clear();
		SetStatistics ssTraining = getSetStatistics(QSPRConstants.TRAINING);
		Set<Long> tsIds = new HashSet<Long>( Globals.session().createQuery("select record.id from BasketEntry be left join be.ep record where be.basket=:basket")
				.setParameter("basket", modelMapping.model.trainingSet).list() );

		if (ssTraining.points.size() != tsIds.size()) {
			// Search for deleted points.. This part can probably be optimized
			for (PointStatistics ps : ssTraining.points)
				if (!tsIds.contains(ps.id)) {
					// logger.info("Records " + ps.id + " has been removed from the basket");
					deletedRecordsIds.add(ps.id);
					ps.deleted = true;
				}
		}

		Set<Integer> currentlyExcludedRecords = modelMapping.model.microattachment.getObject().excludedBasketEntries;

		List<Integer> globalExcludedIds = Globals.session().createCriteria(BasketEntry.class).add(Restrictions.eq("basket", modelMapping.model.trainingSet))
				.add(Restrictions.eq("exclude", Boolean.TRUE)).setProjection(Projections.id()).list();

		currentlyExcludedRecords.addAll(globalExcludedIds);
		Set<Long> currentlyExcludedIds = null;

		if (currentlyExcludedRecords.size() > 0) {
			currentlyExcludedIds = new HashSet<Long>( Globals.session()
					.createQuery("select expProperty.id from BasketEntry be join be.ep expProperty where be.id in (:entriesIds)")
					.setParameterList("entriesIds", currentlyExcludedRecords).list() );
		}

		SetStatistics ssModifiedOnTheFly = new SetStatistics();
		ssModifiedOnTheFly.setId = QSPRConstants.MODIFIED_POINTS_ID;
		ssModifiedOnTheFly.distancesToModel = ssTraining.distancesToModel;

		// Check for post-excluded points
		if (currentlyExcludedIds != null && !currentlyExcludedIds.isEmpty()) {
			Iterator<PointStatistics> iterPoints = ssTraining.points.iterator();
			while (iterPoints.hasNext()) {
				PointStatistics ps = iterPoints.next();
				if (currentlyExcludedIds.contains(ps.id)) {
					ssModifiedOnTheFly.points.add(ps);
					iterPoints.remove();
				}
			}
		}

		// Check for post-included points
		SetStatistics ssExcluded = getSetStatistics(QSPRConstants.EXCLUDED);
		if (ssExcluded != null && !ssExcluded.points.isEmpty()) {
			Iterator<PointStatistics> iterPoints = ssExcluded.points.iterator();
			while (iterPoints.hasNext()) {
				PointStatistics ps = iterPoints.next();
				if (currentlyExcludedIds == null || !currentlyExcludedIds.contains(ps.id)) {
					iterPoints.remove();
					ssModifiedOnTheFly.points.add(ps);
				}
			}

			if (ssExcluded.points.isEmpty())
				sets.remove(ssExcluded);
		}

		if (!ssModifiedOnTheFly.points.isEmpty())
			sets.add(ssModifiedOnTheFly);

		logger.info("Model statistics actualized in " + (Calendar.getInstance().getTimeInMillis() - time) + "ms.");
	}

	@XmlElement(name = "deletedRecordsCount")
	protected Integer getDeletedRecordsCount() {
		if (deletedRecordsIds == null)
			return null;
		return deletedRecordsIds.size();
	}

	@XmlElement(name = "hasDuplicates")
	protected Integer getHasDuplicates() {
		long time = Calendar.getInstance().getTimeInMillis();
		String sets[] = {QSPRConstants.TRAINING,QSPRConstants.MODIFIED_POINTS_ID};
		Set<Long> mol = new HashSet<Long>();
		for(String set:sets) {
			logger.info("Processing: " + set);
			SetStatistics ss = getSetStatistics(set);
			if(ss==null || ss.points == null) {
				logger.info("ModelStat set is absent " + set);
				continue;
			}
			for (PointStatistics ps : ss.points) {
				if(mol.contains(ps.moleculeId))return 1;
				mol.add(ps.moleculeId);
			}
		}
		logger.info("hasDuplicates took " + (Calendar.getInstance().getTimeInMillis() - time) + "ms.");
		return 0;
	}

	public PointStatistics getPointById(Long epId) {
		for (SetStatistics ss : sets)
			for (int i = 0; i < ss.points.size(); i++)
				if (epId.equals(ss.points.get(i).id)) {
					PointStatistics ps = ss.points.get(i);
					ps.parent = ss;
					ps.numInSet = i;
					return ps;
				}

		return null;
	}

	public void deleteValidationSet(Long id) {
		Iterator<SetStatistics> iSets = sets.iterator();
		while (iSets.hasNext()) {
			SetStatistics ss = iSets.next();
			if (ss.setId.startsWith(QSPRConstants.VALIDATION) && (id == null || ss.basketId.equals(id)))
				iSets.remove();
		}
	}

	public SetStatistics getSetStatistics(String id) {
		for (SetStatistics ss : sets)
			if (id.equals(ss.setId))
				return ss;

		return null;
	}

	public SetStatistics getSetStatisticsByBasket(Long basketId) {
		for (SetStatistics ss : sets)
			if (basketId.equals(ss.basketId))
				return ss;

		return null;
	}

	@XmlTransient
	public int getRowsSize() {
		int size = 0;
		for (SetStatistics ss : sets)
			size += ss.points.size();
		return size;
	}
}
