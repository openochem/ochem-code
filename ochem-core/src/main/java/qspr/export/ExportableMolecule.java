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

package qspr.export;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

import qspr.Globals;
import qspr.dao.MetalBondParserSdf;
import qspr.dao.Repository;
import qspr.dao.Various;
import qspr.entities.Article;
import qspr.entities.ConditionSet;
import qspr.entities.ExperimentalProperty;
import qspr.entities.ModelAttachment;
import qspr.entities.ModelMapping;
import qspr.entities.Molecule;
import qspr.entities.MoleculeName;
import qspr.entities.Property;
import qspr.entities.PropertyValue;
import qspr.entities.Unit;
import qspr.modelling.PointStatistics;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.utils.QSPRConstants;

/**
 * Encapsulates a molecule exportable from OCHEM with all the available columns
 * @author midnighter
 *
 */
public class ExportableMolecule
{
	public String smiles;
	public String casRN;
	public Long recordId;
	public String molNum;
	public Integer moleculeId;
	public String name;
	public String name2;
	public String sdf;
	public Article article;
	public String introducer, modifier;
	public String error;
	public int sheet;
	public double molWeight;
	public String inhouseRecordId;
	public String inchiKey;
	public Integer artPageNum;
	public String artTableNum;
	public String artNum;
	public String comments;
	public boolean ownRecord; // This record belongs to the downloading user. He can export it without restrictions
	public boolean ousideOfAD = false;

	public void clearSensitive() {
		moleculeId = null;
		casRN = name = name2 = sdf = inhouseRecordId =inchiKey =  comments = null;
	}

	/**
	 * This is a freely available exportable item. It has no limitations for export
	 */
	public boolean freelyAvailable = false;

	Map<Long, Prediction> predictedValues = new HashMap<Long, Prediction>();

	Map<Property, ExportableValue> expValues = new HashMap<Property, ExportableValue>();
	//	Unit expValueUnit;
	Map<Property, ExportableValue> expValuesConverted = new HashMap<Property, ExportableValue>();

	Map<Property, PropertyValue> conditions = new HashMap<Property, PropertyValue>();
	public AbstractDataRow descriptors;

	protected ExportableSet parent;

	private ThreadLocal<Pattern> casPattern = new ThreadLocal<Pattern>();

	public ExportableMolecule setExperimentalProperty(ExperimentalProperty ep)
	{
		recordId = ep.id;
		if (ep.molecule != null)
			setMolecule(ep.molecule);
		if (ep.artMolId != null)
			molNum = ep.artMolId;
		article = ep.article;
		introducer = ep.introducer != null ? ep.introducer.login : QSPRConstants.ANONYMOUS;
		modifier = ep.owner != null ? ep.owner.login : QSPRConstants.ANONYMOUS;
		inhouseRecordId = ep.externalId;
		comments = ep.other;
		ownRecord = (Globals.userSession().user != null) && Globals.userSession().user.equals(ep.introducer);
		freelyAvailable = ep.rights != null && ep.rights == Globals.RIGHTS_FREELY_AVAILABLE;

		if (ep.moleculenames != null)
			setMoleculeNames(ep.moleculenames);
		//		if (ep.property.isNumeric())
		setExpValue(ep.property, new ExportableValue(ep));
		//		else
		//			setExpValue(ep.property, ep.getStringValue());
		//		expValueUnit = ep.unit;
		artPageNum = ep.artPageNum;
		artTableNum = ep.artTableNum;

		if (ep.conditions != null)
			addConditions(ep.conditions);

		return this;
	}

	/**
	 * Awkward filling from Point
	 * @param point
	 * @param property
	 * @param unit
	 * @param modelAttachment
	 * @return
	 */
	public ExportableMolecule setExperimentalProperty(PointStatistics point, Property property, Unit unit, ModelAttachment modelAttachment) {
		recordId = -1l;
		setMolecule(Repository.molecule.getMolecule(point.moleculeId));
		freelyAvailable = true;
		ExportableValue val = new ExportableValue();
		if (!property.isNumeric())
			val.stringValue = modelAttachment.getOptionFromPrediction(point.real, property).name;
		else {
			val.value = point.real;
			val.unit = unit;
		}
		setExpValue(property, val);
		return this;
	}

	public void setMoleculeNames(List<MoleculeName> names)
	{
		if (parent.selectedColumns.contains(ExportableColumn.NAMES) || parent.selectedColumns.contains(ExportableColumn.CASRN))
		{	
			if (casPattern.get() == null)
				casPattern.set(Pattern.compile("^\\d{2,6}-\\d{2}-\\d$"));
			for (MoleculeName mName : names)
			{
				if (casPattern.get().matcher(mName.name).matches())
				{
					if ( ! parent.selectedColumns.contains(ExportableColumn.CASRN))
						parent.selectedColumns.add(ExportableColumn.CASRN);
					casRN = mName.name;
				}
				else
					if (name == null)
						name = mName.name;
					else
						name2 = mName.name;
			}
		}
	}

	private void setErrorMolecule(){
		moleculeId = 0;
		sdf = QSPRConstants.EMPTY_MOL; // just H atom, to avoid error with the order
		molWeight = 1;
		inchiKey = "";
		smiles = QSPRConstants.ERROR + " Structure cannot be converted to smiles";
	}

	public void setMolecule(Molecule m)
	{
		if(m == null){
			setErrorMolecule();
			return;
		}
		if (!Globals.session().contains(m) && m.id > 0)
			m = Repository.molecule.getMolecule(m.id);
		moleculeId = m.mapping2.id;
		sdf = m.getData();
		if (m.molWeight != null)
			molWeight = m.molWeight;
		inchiKey = m.mapping2.inchi2;

		try
		{
			smiles = 
				Various.molecule.convertToSmilesOrSmart(MetalBondParserSdf.substituteMetalBondwithSingleBond(sdf));
		}catch(Exception e)
		{
			setErrorMolecule();
			e.printStackTrace();
		}
		Globals.session().evict(m);
	}

	public void addConditions(ConditionSet cs)
	{
		for (PropertyValue pv : cs.values)
		{
			if (!parent.conditions.contains(pv.property))
				parent.conditions.add(pv.property);
			conditions.put(pv.property, pv);
		}
	}

	public void setPrediction(ModelMapping mm, Double numericValue)
	{
		getPrediction(mm).predictedValue = numericValue;
	}

	public void setPrediction(ModelMapping mm, String value, Double numericValue)
	{
		getPrediction(mm).predictedValue = value;
		getPrediction(mm).numericPredictedValue = numericValue;
	}

	public void setExpValue(Property p, ExportableValue value)
	{
		if (!parent.properties.contains(p))
			parent.properties.add(p);
		expValues.put(p, value.clone());
		expValuesConverted.put(p, value.clone());
	}

	public void setAccuracy(ModelMapping mm, Double value)
	{
		getPrediction(mm).accuracy = value;
	}

	public void setDM(ModelMapping mm, String dmName, double value)
	{
		List<String> dmNames = parent.dmNames.get(mm);
		if (dmNames == null)
			parent.dmNames.put(mm, dmNames = new ArrayList<String>());
		if (!dmNames.contains(dmName))
			dmNames.add(dmName);
		getPrediction(mm).dms.put(dmName, value);
	}

	public Prediction getPrediction(ModelMapping mm)
	{
		Prediction pred = predictedValues.get(mm.id);
		if (pred == null)
		{
			parent.addModel(mm.model);
			predictedValues.put(mm.id, pred = new Prediction());
		}
		return pred;
	}

	static class Prediction
	{
		/**
		 * Can be a number for regression or a string for classification models
		 */
		Object predictedValue;

		/**
		 * Relevant for classification models
		 */
		Double numericPredictedValue;

		Double accuracy;
		Map<String, Double> dms = new HashMap<String, Double>();
	}

	public class ExportableValue
	{
		String predicate;
		String shortPredicate;

		Double value;
		Double secondValue;
		Unit unit;
		String stringValue;

		public ExportableValue()
		{

		}
		public ExportableValue(ExperimentalProperty ep)
		{
			if (!ep.property.isNumeric())
				stringValue = ep.getStringValue();
			else
			{
				predicate = ep.predicate.name;
				shortPredicate = ep.predicate.shortName;
				value = ep.value;
				secondValue = ep.secondValue;
				unit = ep.unit;
			}
		}

		public String getPredicatedValue()
		{
			if (stringValue != null)
				return stringValue;

			if (predicate == null || predicate.equals("="))
				return "" + value;
			else if (shortPredicate.equals("+-") || shortPredicate.equals("-"))
				return value + " " + predicate+ " " + secondValue;
			else
				return predicate + value;
		}

		public ExportableValue clone()
		{
			ExportableValue ev = new ExportableValue();
			ev.predicate = predicate;
			ev.shortPredicate = shortPredicate;
			ev.value = value;
			ev.secondValue = secondValue;
			ev.unit = unit;
			ev.stringValue = stringValue;
			return ev;
		}


	}


}


