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

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import javax.persistence.CascadeType;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.OneToMany;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.util.ChangesTracker;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.utils.OCHEMUtils;

// TODO : Update variables' names, from condition to property. Refactor it.

@Entity
@XmlRootElement(name = "conditions")
public class ConditionSet implements ChangesTracker
{
	@Id
	@Column(name="con_set_id")
	@GeneratedValue
	@XmlAttribute
	public Long id;

	@OneToMany(mappedBy="conditionSet", cascade={CascadeType.ALL})
	@XmlElement(name = "property-value")
	public List<PropertyValue> values = new ArrayList<PropertyValue>();

	@Column(name = "con_set_md5")
	@XmlTransient
	public String md5;


	@SuppressWarnings("rawtypes")
	public ConditionSet get()
	{
		if(values.size() == 0) return null;

		String _md5 = getMD5();

		List l = Globals.session().createCriteria(ConditionSet.class).add(Restrictions.eq("md5", _md5)).list();
		if (l.size() > 0)
			return (ConditionSet) l.get(0);
		else
		{
			this.md5 = _md5;
			Globals.session().saveOrUpdate(this);
			return this;
		}
	}

	public String getMD5() {
		String hash ="";

		Collections.sort(values, PropertyValue.condNameComp); // sort by name of conditions in the condition set

		for (PropertyValue cv : values)
		{
			cv.conditionSet = this;
			hash += "." + cv.property.id +".";
			if (cv.property.isQualitative())
				hash += cv.option.id;
			else if (cv.property.isNumeric())
			{
				hash += cv.unit.id + "." + cv.value;
				if (cv.predicate != null && !cv.predicate.shortName.equals("="))
					hash += "." + cv.predicate.id;
				if (cv.secondValue != null)
					hash += "." + cv.secondValue;
			}
			else
				hash += cv.textualValue;
		}
		hash = hash.replaceFirst("\\.", "");

		return OCHEMUtils.getMD5(hash);	
	}


	public void addValues(ConditionSet cs)
	{
		for (PropertyValue pv : cs.values) {
			PropertyValue clone = PropertyValue.clone(pv);
			clone.conditionSet = this;
			values.add(clone);
		}
	}

	public String getChanges(ChangesTracker arg) 
	{
		ConditionSet newConditionSet = (ConditionSet)arg;

		List<Property> oldConditions = new ArrayList<Property>();
		List<Property> newConditions = new ArrayList<Property>();

		String log = "";

		for (PropertyValue cv : this.values)
			oldConditions.add(cv.property);

		for (PropertyValue cv : newConditionSet.values)
			newConditions.add(cv.property);			

		oldConditions.removeAll(newConditions);
		if (oldConditions.size() > 0)
			log += "Deleted conditions: "+oldConditions+"\n";
		///
		oldConditions.clear();
		newConditions.clear();

		for (PropertyValue cv : this.values)
			oldConditions.add(cv.property);	

		for (PropertyValue cv : newConditionSet.values)
			newConditions.add(cv.property);

		newConditions.removeAll(oldConditions);

		if (newConditions.size() > 0)
			log += "Added conditions: "+newConditions+"\n";

		oldConditions.clear();
		newConditions.clear();

		for (PropertyValue cv : this.values)
			oldConditions.add(cv.property);	

		for (PropertyValue cv : newConditionSet.values)
			newConditions.add(cv.property);

		oldConditions.retainAll(newConditions);

		for (Property cond : oldConditions) 
		{
			PropertyValue a = null, b = null;
			for (int i=0; i<this.values.size(); i++)
				if (this.values.get(i).property.equals(cond))
					a = this.values.get(i);

			for (int i=0; i<newConditionSet.values.size(); i++)
				if (newConditionSet.values.get(i).property.equals(cond))
					b = newConditionSet.values.get(i);

			if (!a.equals(b))
				log += cond.getName()+" changed from "+a.getStringValue()+" to "+b.getStringValue()+" \n";
		}

		if (!log.equals(""))
			return log;
		else
			return null;
	}

	public PropertyValue getValue(String conditionName)
	{
		for (PropertyValue pv : values)
			if ((Property.shortName(pv.property.getName()).equals(Property.shortName(conditionName))))
				return pv;

		return null;

	}

	public PropertyValue getValue(Property condition)
	{
		for (PropertyValue pv : values)
			if (pv.property.equals(condition))
				return pv;

		return null;
	}

	public PropertyValue getPropertyValue(PropertyValue conditionPropertyValue)
	{
		for (PropertyValue pv : values)
			if (pv.property.id.equals(conditionPropertyValue.property.id))
				return pv;

		return null;
	}

	public PropertyValue getPropertyValue(Long prop_id)
	{
		for (PropertyValue pv : values)
			if (pv.property.id.equals(prop_id))
				return pv;

		return null;
	}

	public ConditionSet mergeAddOrUpdateWith(ConditionSet update)
	{
		// create container set (not to modify existing conditions)
		ConditionSet cs = new ConditionSet();

		// get all existing conditions
		for (PropertyValue opv: this.values)
			cs.values.add(PropertyValue.clone(opv));

		for (PropertyValue upv : update.values)
		{
			// update existing condition
			PropertyValue epv = cs.getValue(upv.property);
			if (epv != null)
				cs.values.remove(epv); // remove existing one if it exists

			cs.values.add(PropertyValue.clone(upv)); // add the update condition 
		}

		return cs.get();
	}

	@Override
	public String toString() {
		String s ="";
		for(PropertyValue v: values)
			s += v + ";\n";
		return s;
	}


	public void addMixtures(String value) throws Exception {
		PropertyValue pv = new PropertyValue();
		pv.property = Property.getByName(QSPRConstants.MIXTURE_CONDITION);
		if(pv.property == null)throw new Exception("Condition " + QSPRConstants.MIXTURE_CONDITION+" is absent");
		pv.textualValue = value;
		pv.predicate = Predicate.get("=");
		values.add(pv);
	}
}
