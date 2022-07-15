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

package qspr.dao;

import java.util.List;

import qspr.entities.Basket;
import qspr.entities.Model;

/**
 * Data access object for Model entity
 * @author midnighter
 *
 */
public interface ModelDAO 
{
	/**
	 * Get a model by its internal ID
	 * @param id
	 * @return Model entity
	 */
	public Model getById(long id);
	
	/**
	 * Get a model by its public ID
	 * @param publicId
	 * @return Model entity
	 */
	public Model getByPublicId(long publicId);
	
	public List<Model> getFeaturedModels();
	
	public List<Model> getPublishedApprovedModels();
	
	public List<Model> getPublishedModels();
	
	public List<Model> getAllModels();

	/**
	 * Get a list of models associated with a specific training set basket 
	 */
	public List<Model> getModelsByTrainingSet(Basket trainingSet);
}
