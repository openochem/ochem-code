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

package qspr.modelling.applier;

import java.util.List;

import org.hibernate.Query;

import qspr.Globals;
import qspr.dao.Repository;
import qspr.entities.Basket;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Molecule;
import qspr.entities.Property;
import qspr.entities.Unit;
import qspr.util.unitconversion.UnitConversion;

public class PropertyPrediction
{
	private double value;
	private double accuracy;
	private double dm;
	private Boolean insideAD;
	private String unit;
	private String property;
	private String moleculeId;
	private String error;
	private Double realValue;
	private String predictedValueString;
	private String realValueString;
	
	public String getPredictedValueString() {
		return predictedValueString;
	}

	public void setPredictedValueString(String predictedValueString) {
		this.predictedValueString = predictedValueString;
	}

	public String getRealValueString() {
		return realValueString;
	}

	public void setRealValueString(String realValueString) {
		this.realValueString = realValueString;
	}
	
	public double getDm()
	{
		return dm;
	}

	public void setDm(double dm)
	{
		this.dm = dm;
	}

	public String getError() {
		return error;
	}

	public void setError(String error) {
		this.error = error;
	}

	public Double getRealValue() {
		return realValue;
	}

	public void setRealValue(Double realValue) {
		this.realValue = realValue;
		this.realValueString = "" + realValue;
	}

	public String getMoleculeId() {
		return moleculeId;
	}

	public void setMoleculeId(String moleculeId) {
		this.moleculeId = moleculeId;
	}

	public double getValue()
	{
		return value;
	}
	
	public void setValue(double value)
	{
		this.value = value;
		this.predictedValueString = "" + value;
	}
	
	public double getAccuracy()
	{
		return accuracy;
	}
	
	public void setAccuracy(double accuracy)
	{
		this.accuracy = accuracy;
	}
	
	public String getUnit()
	{
		return unit;
	}
	
	public void setUnit(String unit)
	{
		this.unit = unit;
	}

	public void setProperty(String property)
	{
		this.property = property;
	}

	public String getProperty()
	{
		return property;
	}	
	
	public Boolean getInsideAD()
	{
		return insideAD;
	}

	public void setInsideAD(Boolean insideAD)
	{
		this.insideAD = insideAD;
	}
	
	public String toString()
	{
		return "Property:"+property+", Value:"+value+", Accuracy:"+accuracy+", Unit:"+unit;
	}

	public void addRealValues(Property property, Molecule molecule, Unit targetUnit, Basket trainingSet) {
		String qq =  "select exp_property_id from ExperimentalProperty  natural join BasketEntry  natural join Molecule where basket_id=" +  trainingSet.id +
				" and property_id=" + property.id + "  and mapping2_id=" + molecule.mapping2.id;
		
		System.out.println(qq);
		
	    Query q = Globals.session().createSQLQuery(qq);
	    @SuppressWarnings("unchecked")
		List<Integer> entities = q.list();
	    if(entities.size() == 0) return;

	    long id = entities.get(0);
	    
	    ExperimentalProperty ep = Repository.record.getRecord(id);
	    
	    if(property.isQualitative())
	    	this.realValueString = Repository.option.getPropertyOptionById(ep.option.id).name;
	    else
	    	this.realValueString = "" + UnitConversion.convert(ep.value, ep.unit, targetUnit, molecule.molWeight);
	}
}
