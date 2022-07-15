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

import org.hibernate.criterion.Order;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.entities.Basket;
import qspr.entities.Model;


@SuppressWarnings("unchecked")
public class ModelDAOImpl implements ModelDAO {

	@Override
	public Model getById(long id) {
		return (Model) Globals.session().get(Model.class, id);
	}

	@Override
	public Model getByPublicId(long publicId) {
		List<Model> models = Globals.session().createCriteria(Model.class).add(Restrictions.eq("publicId", publicId)).list();
		if (models.size() > 0)
			return models.get(0);
		return null;
	}

	@Override
	public List<Model> getFeaturedModels() {
		return Globals.session().createCriteria(Model.class)
				.add(Restrictions.eq("published", true))
				.add(Restrictions.isNotNull("featuredName"))
				.addOrder(Order.asc("publicId"))
				.list();
	}

	@Override
	public List<Model> getModelsByTrainingSet(Basket trainingSet)
	{
		return Globals.session().createCriteria(Model.class)
				.add(Restrictions.eq("trainingSet", trainingSet))
				.list();
	}

	@Override
	public List<Model> getPublishedApprovedModels() {
		return Globals.session().createCriteria(Model.class)
				.add(Restrictions.eq("published", true))
				.add(Restrictions.eq("approved", true))
				.addOrder(Order.asc("publicId"))
				.list();	
	}

	@Override
	public List<Model> getPublishedModels() {
		return Globals.session().createCriteria(Model.class)
				.add(Restrictions.eq("published", true))
				.addOrder(Order.asc("publicId"))
				.list();	
	}

	@Override
	public List<Model> getAllModels() {
		return Globals.session().createCriteria(Model.class)
				.list();	
	}
}
