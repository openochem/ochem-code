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

import org.hibernate.Criteria;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.entities.ModelTemplate;

public class ModelTemplateDAOImpl implements ModelTemplateDAO {

	@SuppressWarnings("unchecked")
	@Override
	public ModelTemplate getByName(String name) {
		Criteria c = Globals.session().createCriteria(ModelTemplate.class);
		c.add(Restrictions.eq("name", name));
		List<ModelTemplate> templates = (List<ModelTemplate>)c.list();
		if (templates.size()>0)
			return templates.get(0);
		else
			return null;
	}

}
