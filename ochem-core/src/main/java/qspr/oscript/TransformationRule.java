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

package qspr.oscript;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.antlr.runtime.ANTLRInputStream;
import org.antlr.runtime.CommonTokenStream;
import org.antlr.runtime.RecognitionException;

import qspr.Globals;
import qspr.ThreadScope;
import qspr.dao.Repository;
import qspr.entities.ConditionSet;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Property;
import qspr.entities.PropertyOption;
import qspr.entities.PropertyValue;

import com.eadmet.exceptions.UserFriendlyException;

public class TransformationRule 
{
	//	private static transient final Logger logger = Logger.getLogger(TransformationRule.class);

	List<AtomicRule> rules = new ArrayList<AtomicRule>();
	public int countMatches = 0;

	public void addRule(AtomicRule rule)
	{
		rules.add(rule);
	}

	public ExperimentalProperty apply(ExperimentalProperty ep)
	{
		ExperimentalProperty transformed = ep.fullClone();
		transformed.moleculenames = null;
		transformed.rights = Globals.RIGHTS_NONE;
		transformed.connectedProperty = ep;
		transformed.owner = transformed.introducer = ThreadScope.get().userSession.user;

		boolean wasModified = false;

		for (AtomicRule rule : rules)
			if (rule.transform(transformed))
			{
				wasModified = true;
				countMatches++;
			}

		if (!wasModified)
			return null;

		return transformed;
	}

	public static TransformationRule compile(String script) throws IOException, RecognitionException
	{
		final ANTLRInputStream input = 
				new ANTLRInputStream(new ByteArrayInputStream(script.getBytes()));

		OScriptParser parser = new OScriptParser(
				new CommonTokenStream(new OScriptLexer(input)));

		TransformationRule rule = parser.start();
		if (parser.failed() || parser.getNumberOfSyntaxErrors() > 0)
			throw new UserFriendlyException("Failed to compile the script. Please, check the syntax");

		return rule;
	}

	public RuleValidationResult validate()
	{
		RuleValidationResult validation = new RuleValidationResult();

		for (AtomicRule rule : rules)
			rule.validate(validation);

		return validation;
	}
}


class AtomicRule
{
	//private static transient final Logger logger = Logger.getLogger(AtomicRule.class);
	public Specifier left;
	public Specifier right;
	public String rawRule;

	public boolean transform(ExperimentalProperty ep)
	{
		if (!left.matches(ep))
			return false;
		else 
		{
			//logger.info("Matched rule " + rawRule);
			ep.property = right.property;

			if(ep.option != null) { //N.B.! Very long, can be speeded up: we have to re-map options to the new property or add them 
				ep.option = Repository.option.getPropertyOptionByName(ep.option.name, ep.property.id, true, true);
			}

			if (right.condition != null)
			{
				ConditionSet cSet = new ConditionSet();
				if (ep.conditions != null)
					cSet.addValues(ep.conditions);
				if (right.option != null)
					cSet.values.add(new PropertyValue(right.option));
				else
					cSet.values.add(new PropertyValue(right.condition, right.number, right.condition.defaultUnit));
				ep.conditions = cSet.get();
			}
		}

		return true;
	}

	public void validate(RuleValidationResult validation)
	{
		left.validate(validation);
		right.validate(validation);
		
		if(!left.property.unitCategory.equals(right.property.unitCategory)) { // new property ! Need to be saved
			//right.property = Repository.property.getPropertyById(right.property.id);
			right.property.unitCategory = left.property.unitCategory;
			right.property.defaultUnit = left.property.defaultUnit;
			right.property.type = left.property.type;
			right.property.owner = right.property.introducer = ThreadScope.get().userSession.user;
			right.property.conditionsUsed = left.property.conditionsUsed;
			right.property.obligatoryConditions.addAll(left.property.obligatoryConditions);
			Globals.session().saveOrUpdate(right.property);
		}
		
		if (left.property != null && right.property != null && !left.property.unitCategory.equals(right.property.unitCategory))
			validation.critical("Incompatible properties: " + left.propertyName + " and " + right.propertyName);
	}
}

class Specifier
{
	public String propertyName;
	public String conditionName;
	public String value;
	public String predicate;
	public String unitName;

	public Property property;
	public Property condition;
	public PropertyOption option;
	public Double number;

	public Specifier(String property)
	{
		this.propertyName = property.replaceAll("\"", "");;
	}

	public Specifier(String property, String condition, String predicate, String value, String unitName)
	{
		this.propertyName = property.replaceAll("\"", "");
		this.conditionName = condition;
		this.predicate = predicate;
		this.value = value;
		this.unitName = unitName;
	}

	public boolean matches(ExperimentalProperty ep)
	{
		if (!propertyName.equals("*") && !ep.property.getName().equals(propertyName)) // property matched?
			return false;

		if (condition == null)
			return true;

		if (!matches(ep.conditions.getValue(condition)))
			return false;

		return true;
	}

	public boolean matches(PropertyValue pv)
	{
		if (pv == null)
			return false;
		if (pv.property.isNumeric())
			return pv.value.equals(new Double(value));
		else
			return pv.option.name.equals(value);
	}

	public void validate(RuleValidationResult validation)
	{
		property = Repository.property.getProperty(propertyName, true);
		
//		if (!propertyName.equals("*") && (property = Repository.property.getProperty(propertyName, true)) == null)
//			validation.unknownProperties.add(propertyName);
		
		if (conditionName != null)
		{
			condition = Repository.property.getProperty(conditionName, false);
			if (condition == null)
				validation.unknownConditions.add(conditionName);
			else
			{
				if (condition.isQualitative())
				{
					option = condition.getOptionByName(value);
					if (option == null)
						validation.critical("\""+value+"\" is not a valid option for condition \"PROP$"+condition.id+"\"");
					if (!predicate.equals("="))
						validation.critical("PROP$"+condition.id  + " is a qualitative condition! Are you OK? Only 'equals' are supported!");
				}
				else
				{
					try
					{
						number = Double.valueOf(value);
					}
					catch (Exception e)
					{
						validation.critical("" + value + " is not a valid number.");
					}
				}
			}
		}
	}
}


