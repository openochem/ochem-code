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

package qspr.frontend;

import java.util.List;

import javax.xml.bind.annotation.XmlRootElement;

import org.hibernate.criterion.Order;
import org.hibernate.criterion.Projections;

import qspr.Globals;
import qspr.entities.Article;
import qspr.entities.Property;
import qspr.entities.SubstructureAlert;

@XmlRootElement(name = "available-alerts")
public class AvailableAlertsData
{
	public List<Article> articles;
	public List<Property> endpoints;
	public int selectionSize;

	@SuppressWarnings("unchecked")
	public AvailableAlertsData()
	{
		articles = Globals.session()
				.createCriteria(SubstructureAlert.class)
				.createAlias("article", "a")
				.setProjection(Projections.groupProperty("article")).addOrder(Order.asc("a.publicationDate")).list();
		endpoints = Globals.session()
				.createCriteria(SubstructureAlert.class)
				.setProjection(Projections.groupProperty("property")).list();
		selectionSize = Globals.userSession().selectedAlerts.size();
	}
}
