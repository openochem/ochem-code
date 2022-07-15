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

import org.apache.commons.lang.StringUtils;

import qspr.entities.Model;
import qspr.modelling.applier.ModelApplier;
import qspr.modelling.applier.ModelApplierTaskProcessor;

@XmlRootElement
public class PredictionAction extends AbstractUserAction
{
	public List<Model> models;
	public long compoundsSize;
	
	@Override
	public String getLogLine()
	{
		List<String> names = new ArrayList<String>();
		for (Model m : models)
			names.add(m.name);
		return String.format("is predicting %d compounds using model %s", compoundsSize, StringUtils.join(names.toArray(), ", "));
	}
	
	public PredictionAction(ModelApplier applier)
	{
		models = new ArrayList<Model>();
		for (ModelApplierTaskProcessor mt : applier.modelTasks)
		{
			Model m = new Model();
			m.name = mt.model.name;
			m.id = mt.model.id;
			models.add(m);
		}
		
		compoundsSize = applier.compoundsProvider.getCompoundsNum();
	}
	
	public PredictionAction()
	{
		
	}

}
