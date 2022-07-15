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

package qspr.entities;

import java.util.Comparator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.FetchType;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.Transient;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import com.eadmet.exceptions.UserFriendlyException;


@Entity
@XmlRootElement(name = "property-value")
public class PropertyValue 
{	

	@Id
	@GeneratedValue
	@Column(name = "propertyvalue_id")
	@XmlAttribute
	public Long id;

	@XmlTransient
	@ManyToOne(fetch=FetchType.LAZY)
	@JoinColumn(name = "con_set_id")
	public ConditionSet conditionSet;

	@ManyToOne
	@JoinColumn(name = "property_id")	
	public Property property;	

	@ManyToOne
	@JoinColumn(name = "poption_id")
	public PropertyOption option;

	@Column
	@XmlElement
	public Double value;

	@Column(name = "second_value")
	@XmlElement
	public Double secondValue; 

	//	@XmlAttribute
	//	@Column(name = "predicate")
	//	public String oldPredicate; // =, <, >, +-, -

	@ManyToOne
	@JoinColumn(name = "predicate_id")
	public Predicate predicate;

	@ManyToOne
	@JoinColumn(name = "unit_id")
	@XmlElement
	public Unit unit;

	@Transient
	@XmlTransient
	public Long old_id;

	@XmlAttribute
	@Transient
	public boolean multi; // means different values for one quantitative property / condition value in batch editing

	@Column(name = "text_value")
	@XmlElement
	public String textualValue;

	public PropertyValue()
	{

	}

	public PropertyValue(PropertyOption option)
	{
		this.property = option.property;
		this.option = option;
		this.unit = this.property.defaultUnit;
	}

	public PropertyValue(Property property, String textualValue)
	{
		this.property = property;
		this.textualValue = textualValue;
	}

	public PropertyValue(Property cond, Double value)
	{
		this(cond, value, null);
	}

	public PropertyValue(Property condition, Double condition_value, Unit condition_unit)
	{
		property = condition;
		value = condition_value;
		if (condition_unit == null)
			unit = property.defaultUnit;
		else 
			unit = condition_unit;
	}

	public PropertyValue(Property condition, Predicate pred, Double condition_value, Double second_value, Unit condition_unit)
	{
		property = condition;
		value = condition_value;
		secondValue = second_value;
		predicate = pred;

		if (condition_unit == null)
			unit = property.defaultUnit;
		else 
			unit = condition_unit;
	}

	public static PropertyValue clone(PropertyValue original)
	{
		PropertyValue epv = new PropertyValue();
		epv.value = original.value;
		epv.secondValue = original.secondValue;
		epv.predicate = original.predicate;
		epv.option = original.option;
		epv.unit = original.unit;
		epv.textualValue = original.textualValue;
		epv.multi = original.multi;
		epv.property = original.property;
		return epv;
	}

	public void setValueWithPredicate(String newValue) throws Exception
	{

		try 
		{
			Object[] vals = PropertyValue.parseTextIntoPredicateAndValues(newValue);
			predicate = (Predicate)vals[0];
			value = (Double)vals[1];
			secondValue = (Double)vals[2];
		}
		catch (Exception e)
		{
			if (newValue == null || newValue.trim().equals(""))
				throw new UserFriendlyException("No value for property \""+property.getName()+"\" provided");
			else
				throw new UserFriendlyException("Invalid value for property \""+property.getName()+"\": "+newValue);	
		}
	}

	public static Object[] parseTextIntoPredicateAndValues(String textValue) throws Exception
	{
		if (textValue == null)
			throw new Exception("Error parsing empty value");

		textValue = textValue.replaceAll("\\s+", "");

		if (textValue.equals(""))
			throw new Exception("Error parsing empty value");
		//TODO: Allowed predicates are hard-coded in this pattern... we can form the pattern on the fly theoretically
		String predicatePattern = "=|>|<|>>|<<|>=|<=|~|~="; 
		String decimalValuePattern = "([+-]?\\d*)(\\.\\d*)?([eE][+-]?\\d+)?";
		Pattern p = Pattern.compile("("+predicatePattern+")?("+decimalValuePattern+")((-|\\+-)("+decimalValuePattern+"))?");
		Matcher m = p.matcher(textValue);
		if (!m.matches())
			throw new Exception("Unrecognized format for property value");
		String predicate = m.group(1);
		String value1 = m.group(2);
		String subpredicate = m.group(7);
		String value2 = m.group(8);

		if (value2 != null)
			if (predicate != null)
				if (!predicate.equals("="))
					throw new Exception("Predicates and range/accuracy expressions are not allowed together");

		if (value2 != null)
		{
			Object[] result = {Predicate.get(subpredicate), Double.valueOf(value1), Double.valueOf(value2)};
			return result;
		} else
		{
			if (predicate == null)
				predicate = "=";

			Object[] result = {Predicate.get(predicate), Double.valueOf(value1), null};
			return result;
		}

	}

	@XmlAttribute(name = "printableValueFull")
	public String getStringValueFull()
	{
		if (property == null)
			return null;

		if (property.isQualitative())
			return " = " + option.name;

		if (property.isTextual())
			return " = "+ textualValue;

		String res;

		if (predicate == null)
			res = " = "+value;
		else
			if (predicate.shortName.equals("+-") || predicate.shortName.equals("-"))
				res = " = "+value+" "+predicate.name+" "+secondValue;
			else
				res = predicate.name+" "+value;

		res += " "+unit.getName();

		return res;
	}

	@XmlAttribute(name = "printableValue") //No equal as predicate, no unit
	public String getStringValue()
	{
		if (property == null)
			return null;

		if (property.isQualitative())
			return option.name;

		if (property.isTextual())
			return textualValue;

		String res;

		if (predicate == null)
			res = ""+value;
		else
			if (predicate.shortName.equals("+-") || predicate.shortName.equals("-"))
				res = value+" "+predicate.name+" "+secondValue;
			else
				res = predicate.name+" "+value;

		return res;
	}

	public String toString()
	{
		return property.getName() + ":" + getStringValueFull();
	}

	public boolean isValid()
	{
		if (null != property)
		{
			if (property.isQualitative())
				return null != option;
			else if (property.isTextual())
				return textualValue != null;
			else
				return null != unit && null != value;
		}
		return false;
	}


	@Override
	public int hashCode()
	{
		return property.id.hashCode();
	}

	public boolean equals(PropertyValue pv)
	{
		if  (!pv.property.id.equals(property.id))
			return false;

		if (property.isQualitative())
			return pv.option.id.equals(option.id);
		else if (property.isTextual())
			return textualValue.equalsIgnoreCase(pv.textualValue);
		else
			return pv.value.equals(value) && pv.unit.equals(unit) && ((pv.predicate == null) ? (predicate == null) : (pv.predicate.equals(predicate))) && ((pv.secondValue == null) ? (secondValue == null) : (pv.secondValue.equals(secondValue)));
	}

	public static Comparator<PropertyValue> condNameComp = new Comparator<PropertyValue>() {
		public int compare(PropertyValue pv1, PropertyValue pv2) {
			return pv1.property.getName().compareTo(pv2.property.getName());
		}
	};
}