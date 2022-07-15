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

import java.io.IOException;
import java.io.Serializable;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import qspr.Globals;
import qspr.MarshallingOption;
import qspr.dao.Repository;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Model;
import qspr.entities.ModelMapping;
import qspr.util.unitconversion.UnitConversion;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.NumericalValueStandardizer;

// This class represents a "model dot". It is stored in serialiable form as a part of ModelStatistics.SetStatistics
// In future, rather than to store this info as a Serializable, it would be useful to make it a database entity to allow easy filtering
// Midighter/modified by IVT

@XmlRootElement(name = "point")
public class PointStatistics implements Serializable {
	private static final long serialVersionUID = 2L;

	public double real; // its an average value in case of an interval
	public double real2; // the right boundary of an interval
	public double predicted;

	/**
	 * Prediction values returned by an ensebmle of models
	 */
	public float[] ensemblePredictions;

	public String predicate; // we need it for calculating statistics

	public long id; // this is experimental property id in ExperimentalProperty

	//TODO change to mapping2 id
	@XmlTransient // currently this is mapping1
	public Long moleculeId;

	// For dynamic highlighting of the points on the graph, filtered by some criteria
	@XmlAttribute
	public Boolean selected;

	// @XmlElement(name = "dm")
	@XmlTransient
	public List<Double> distancesToModel = new ArrayList<Double>();

	@XmlElement
	protected ADData ad; // For XML only. Calculated on the fly when necessary, not stored

	@XmlAttribute
	public String error;

	@XmlAttribute
	public Boolean deleted;

	@XmlAttribute
	public Boolean virtual;

	public transient SetStatistics parent;

	public transient int numInSet;
	public transient Double molWeight;

	@XmlElement(name = "moleculeId")
	public Long getMoleculeId() {
		if (moleculeId != null)
			return moleculeId;
		else {
			return moleculeId = (Long) Globals.session()
					.createQuery("select mp1.id from ExperimentalProperty as ep join ep.molecule as mol join mol.mapping1 as mp1 where ep.id=:id")
					.setParameter("id", id).uniqueResult();
		}
	}

	@XmlElement(name = "articleId")
	public Long getArticleId() {
		if (!Globals.getMarshallingOption(MarshallingOption.ARTICLE_IN_MODEL_DOT))
			return null;

		return (Long) Globals.session().createQuery("select art.id from ExperimentalProperty as ep join ep.article as art where ep.id=:id")
				.setParameter("id", id).uniqueResult();
	}

	public void setADData(Double dmValue, Double error) {
		ad = new ADData();
		ad.setValue(dmValue);
		ad.setError(error);
	}

	public static class DMComparator implements Comparator<PointStatistics> {
		public int compare(PointStatistics p1, PointStatistics p2) {
			return p1.ad.compare(p2.ad);
		}
	}

	public String getMoleculeHashWithValue(){
		return moleculeId + " " + predicted;
	}

	public PointStatistics() {
	}

	public PointStatistics(double predicted) {
		this.predicted = NumericalValueStandardizer.getSignificantDigitsDouble(predicted,NumericalValueStandardizer.SIGNIFICANT_DIGITS);
	}

	public double getMolWeight() throws IOException {
		if (molWeight == null)
			molWeight = Repository.record.getWeight(id); // by recordid
		return molWeight;
	}

	public PointStatistics(PointStatistics ps) {
		predicted = 0;
		id = ps.id;
		moleculeId = ps.moleculeId; //TODO Simplify weight assignment
		try {
			molWeight = ps.getMolWeight();
		}catch(IOException e) {
			molWeight = 0.0;
		}
		real = ps.real;
		real2 = ps.real2;
		predicate = ps.predicate;
		virtual = ps.virtual;
	}

	public PointStatistics(ExperimentalProperty ep, ModelMapping modelMapping) {
		predicted = 0;
		id = ep.id;
		moleculeId = ep.molecule.mapping1.id;
		molWeight = ep.molecule.molWeight;

		if (modelMapping.property.isQualitative())
		{
			Long realVal = modelMapping.model.getMappedOption(ep.option.id);
			if (realVal != null)
				real = realVal;
			else
				// This can happen if a validation set contains an option not present in the training set. We have to discuss how to treat such cases 
				throw new UserFriendlyException("A set contains an unsupported property option: \"" + ep.option.name + "\" for property \"" + ep.property.getName()+ "\"");
		}
		else
			try {
				real = ep.getConvertedAverageValue(modelMapping.unit);
				if (ep.predicate.isInterval())
					real2 = ep.getConvertedValue(ep.secondValue, modelMapping.unit);

				predicate = ep.predicate.shortName;
				if (ep.predicate.isOrdering() && ep.unit.isOppositeTo(modelMapping.unit))
					// invert the predicate
					predicate = ep.predicate.getOpposite();
				else if (ep.predicate.isInterval() && real2 < real) {
					// swap interval values
					double tmp = real;
					real = real2;
					real2 = tmp;
				}
			} catch (Exception e) {
			}
	}

	public double getRealValue(boolean considerPredicates) {
		if (considerPredicates) {
			// If we are within accepted range, consider that real = predicted
			if (predicate != null && !predicate.equals("="))
				if (predicate.startsWith("<") && (predicted < real))
					return predicted;
				else if (predicate.startsWith(">") && (predicted > real))
					return predicted;
				else if (predicate.equals("-") && (predicted > real && predicted < real2)) //
					return predicted;
		}

		return real;
	}

	public void transferDataToEP(ExperimentalProperty ep, Model model, ModelMapping modelMapping, SetStatistics ss)
	{
		DecimalFormat formatter = (DecimalFormat)NumberFormat.getInstance();
		formatter.applyPattern("###,##0.00");
		ep.real = formatter.format(real);
		ep.predicted = formatter.format(predicted);
		ep.modelUnit = modelMapping.unit;
		if (ep.property.isNumeric())
		{
			String currentValueInModelUnits = formatter.format(UnitConversion.convert(ep.value, ep.unit, modelMapping.unit, ep.molecule.molWeight));
			ep.valueInModelUnits = formatter.format(real);

			if (!currentValueInModelUnits.equals(ep.valueInModelUnits))
				ep.valueInModelUnits = currentValueInModelUnits + "_" + ep.valueInModelUnits;
		}
		if (modelMapping.property.isQualitative())
			ep.predictedInOriginalUnits = model.attachment.getObject().getOptionFromPrediction(predicted, modelMapping.property).name;
		else
			ep.predictedInOriginalUnits = formatter.format(UnitConversion.convert(predicted, modelMapping.unit, ep.unit, ep.molecule.molWeight));
		if (distancesToModel.size() > 0)
			ep.setDMValue(ss.distancesToModel.get(0), formatter.format(distancesToModel.get(0)));

		ep.predictionVector = new ArrayList<Float>();
		if (ensemblePredictions != null)
			for (float prediction : ensemblePredictions)
				ep.predictionVector.add(prediction);
	}
}

class ADData {
	@XmlAttribute
	public String error;

	@XmlAttribute
	public String dmValue;

	private float dm;

	void setValue(double val){
		this.dm = (float) val;
		dmValue = NumericalValueStandardizer.getSignificantDigits(val);
	}

	public int compare(ADData ad) {
		return dm > ad.dm ? 1 : dm == ad.dm ? 0 : 1;
	}

	void setError(double val){
		error = NumericalValueStandardizer.getSignificantDigits(val);
	}


}