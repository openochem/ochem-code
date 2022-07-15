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

package com.eadmet.mmpa;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import qspr.dao.Repository;
import qspr.entities.Model;

public class MMPOptmiserSetup implements Serializable 
{
	private static final long serialVersionUID = 1L;
	public TransformationFilter filter;
	public List<Long> modelIds;
	public List<Long> transformationIds;
	public Set<Long> inverseTransformationIds;
	
	public Short similarity;
	public Integer rounds;
	
	public List<Model> getModels() {
		List<Model> models = new ArrayList<Model>();
		for (Long id : modelIds)
			models.add(Repository.model.getByPublicId(id));
		
		return models;
	}
}