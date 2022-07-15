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

package com.eadmet.useractions;

import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.annotation.XmlRootElement;

import qspr.entities.Basket;
import qspr.entities.Model;
import qspr.entities.ModelMapping;
import qspr.entities.Property;
import qspr.modelling.configurations.CDSConfiguration;

/**
 * A documentation data for a model creation action
 * @author midnighter
 *
 */
@XmlRootElement
public class ModelCreationAction extends AbstractUserAction
{
	public List<Property> properties;
	public String method;
	public String descriptors;
	public long tsID;
	public long tsSize;
	public Long modelID;

	public int modelsCount;

	public ModelAction action = ModelAction.START; 


	@Override
	public String getLogLine()
	{
		if (modelsCount == 0)
			if (method != null && method.startsWith("Descriptors"))
				return "has started calculation of descriptors " + descriptors + " for " + tsSize + " compounds";
			else
			{
				if (action == ModelAction.START)
					return "has started model creation for " + properties.toString() + " using method " + method + " and descriptors " + descriptors + " using a training set of " + tsSize + " compounds";
				else if (action == ModelAction.DELETE)
					return " has deleted a model for " + properties;
				else if (action == ModelAction.RECALCULATE)
					return " has started recalculation of a model for " + properties;
				else if (action == ModelAction.SAVE)
					return " has saved a private model for " + properties;
				else
					return "";

			}
		else
			return String.format("has started creation of %d models using a training set with %d compounds", modelsCount, tsSize);
	}

	public ModelCreationAction()
	{

	}

	public ModelCreationAction(Basket ts, int modelsCount)
	{
		this.modelsCount = modelsCount;
		properties = new ArrayList<Property>();

		for (Property property : ts.getProperty())
		{
			Property p = new Property();
			p.id = property.id;
			p.setName(property.getName());
			properties.add(p);
		}

		tsID = ts.id;
		tsSize = ts.getRowsSize();
	}

	public ModelCreationAction(Model model, ModelAction action)
	{
		modelID = model.id;
		this.action = action;
		properties = new ArrayList<Property>();
		method = model.template.name; 

		try {
			if (model.attachment.getObject().configuration instanceof CDSConfiguration)
				descriptors = ((CDSConfiguration) model.attachment.getObject().configuration).descriptors.types.toString();
		}catch(RuntimeException e) {}

		for (ModelMapping mm : model.modelMappings)
		{
			Property p = new Property();
			p.id = mm.property.id;
			p.setName(mm.property.getName());

			properties.add(p);
		}

		tsID = model.trainingSet.id;
		tsSize = model.trainingSet.getRowsSize();
	}

	public static enum ModelAction{
		START, DELETE, SAVE, RECALCULATE
	};

}
