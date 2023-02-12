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

package qspr.business;


import java.util.List;

import org.hibernate.Criteria;
import org.hibernate.criterion.Disjunction;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.entities.Basket;
import qspr.entities.Model;
import qspr.entities.User;
import qspr.util.CriteriaWrapper;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.exceptions.UserFriendlyException;

public class ModelOperation
{
	public Criteria getFilterCriteria(ModelFilter filter)
	{
		Criteria c = Globals.session().createCriteria(Model.class);
		CriteriaWrapper cw = new CriteriaWrapper(c);
		if (filter.proQuery != null)
		{
			Disjunction d = Restrictions.disjunction();
			d.add(Restrictions.like("name", "%" + filter.proQuery + "%"));

			cw.createAlias("modelMappings", "mm");
			cw.createAlias("mm.property", "property");
			d.add(Restrictions.like("property.name", "%" + filter.proQuery + "%"));

			c.add(d);
		}

		if (filter.publicId != null)
			c.add(Restrictions.eq("publicId", filter.publicId));

		if (filter.toBeDeleted)
			c.add(Restrictions.isNotNull("deleteAfterDays")).add(Restrictions.eqOrIsNull("published", false));

		c.add(Restrictions.isNull("taskId"));
		Model.addAuthRestrictions(c, filter.groupModels, filter.awaitingApproval);

		return c;
	}

	public static void approveModel(Model model, boolean publishedAndCited, int qualityGrade) {
		if (model.approved)
			throw new UserFriendlyException("The model is already approved!");

		User publisher = User.getById(QSPRConstants.PUBLISHER_ID);

		model.trainingSet.user=publisher;
		Globals.session().saveOrUpdate(model.trainingSet); 

		for(Basket b:model.getValidationSets()) {
			b.user = publisher;
			Globals.session().saveOrUpdate(b); 
		}
		model.approved = true;
	}

	@SuppressWarnings("unchecked")
	public int processMarkedModels(boolean deleteAll) {
		ModelFilter filter = new ModelFilter();
		filter.toBeDeleted = true;
		List<Model> models = getFilterCriteria(filter).list();

		for (Model model : models)
			if (deleteAll)
				model.delete();
			else
				model.markAccessed();

		return models.size();
	}
}
