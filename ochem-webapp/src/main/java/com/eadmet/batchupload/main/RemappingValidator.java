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

package com.eadmet.batchupload.main;

import qspr.entities.Article;
import qspr.entities.Property;
import qspr.entities.PropertyOption;
import qspr.entities.Unit;
import qspr.util.UploadContext;

import com.eadmet.batchupload.entityschema.AbstractPropertyRemapping;
import com.eadmet.batchupload.entityschema.ConditionRemapping;
import com.eadmet.batchupload.entityschema.EntitiesRemapping;
import com.eadmet.batchupload.entityschema.PropertyRemapping;
import com.eadmet.batchupload.entityschema.RemappedValue;

public class RemappingValidator 
{
	public static void validate (EntitiesRemapping remapping, UploadContext context)
	{
		for (RemappedValue article : remapping.articles) 
			validateArticle(article, context);

		for (PropertyRemapping property : remapping.properties)
		{
			validateAbstractProperty(property, context);
			for (ConditionRemapping condition : property.conditions)
				validateAbstractProperty(condition, context);
		}
	}

	private static void validateArticle(RemappedValue article, UploadContext context)
	{
		article.messages.clear();
		if (article.name.equals("unpublished"))
			return;
		try {
			Article a = Article.getArticle(article.name, context);
			if (a == null)
				article.messages.add(BatchUploadMessage.newError("Article "+article.name+" was not found in the database"));
		} catch (Exception e)
		{
			article.messages.add(BatchUploadMessage.newError(e));
		}
	}

	private static void validateAbstractProperty(AbstractPropertyRemapping property, UploadContext context)
	{
		property.messages.clear();
		try {
			Property p = Property.getByName(property.name);
			if (p == null)
				property.messages.add(BatchUploadMessage.newError("Property "+property.name+" was not found in the database"));
			else
			{
				if (p.isNumeric())
				{
					if (property.options.size() > 0)
						property.messages.add(BatchUploadMessage.newError("Property "+property.name+" is numeric, but you try to upload qualitative values, e.g. " + property.options.get(0).name));
					for (RemappedValue unit : property.units) 
						validateUnit(unit, p, context);
				} else
					if (p.isQualitative())
					{
						if (property.units.size() > 0)
							property.messages.add(BatchUploadMessage.newError("Property "+property.name+" is qualitative,  but you try to upload numeric values, e.g. "
									+ property.units.get(0).minValue + (property.units.get(0).minValue != property.units.get(0).maxValue ? " " + property.units.get(0).maxValue : "" )));
						for (RemappedValue option : property.options) 
							validateOption(option, p, context);
					}
			}
		} catch (Exception e)
		{
			property.messages.add(BatchUploadMessage.newError(e));
		}

	}

	private static void validateUnit(RemappedValue unit, Property p, UploadContext context)
	{
		unit.messages.clear();
		if (unit.name.equals("default"))
			return;
		try {
			Unit u = Unit.getByNameAndCategory(unit.name, p.unitCategory.name, false);
			if (u == null)
				unit.messages.add(BatchUploadMessage.newError("Unit "+unit.name+" was not found for category "+p.unitCategory.name));
		} catch (Exception e)
		{
			unit.messages.add(BatchUploadMessage.newError(e));
		}
	}

	private static void validateOption(RemappedValue option, Property p, UploadContext context)
	{
		option.messages.clear();
		try {
			PropertyOption o = p.getOptionByName(option.name);
			if (o == null)
				option.messages.add(BatchUploadMessage.newError("Option "+option.name+" was not found for property "+p.getName()));
		} catch (Exception e)
		{
			option.messages.add(BatchUploadMessage.newError(e));
		}		
	}
}
