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

import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.entities.SubstructureAlert;
import qspr.fragmententities.Fragment;
import qspr.workflow.utils.QSPRConstants;
import qspr.workflow.datatypes.AbstractDataRow;

/**
 * A class to represent a under-represented or over-represented descriptor in SetCompare UI
 * @author midnighter
 *
 */
@XmlRootElement(name = "setcompare-result")
public class SetCompareResult
{
	public int id;
	public String depictionPath;
	public String href;
	public String name;
	public String displayName;
	public int inSet1;
	public int inSet2;
	public double pValue1;
	public double pValue2;

	public transient Long alertId;

	public SetCompareResult setId(int id)
	{
		this.id = id;
		return this;
	}

	/**
	 * Create from the resulting DataRow
	 * @param row
	 */
	@SuppressWarnings("unchecked")
	public SetCompareResult(AbstractDataRow row)
	{
		displayName = name = (String) row.getValue(0);
		inSet1 = new Integer("" + row.getValue(1));
		inSet2 = new Integer("" + row.getValue(2));
		pValue1 = new Double("" + row.getValue(3));
		pValue2 = new Double("" + row.getValue(4));

		if (name.startsWith("Alert"))
		{
			String smart = name.split("_", 2)[1];
			List<SubstructureAlert> alerts = 
					Globals.session().createSQLQuery("select * from SubstructureAlert where smart like :smart")
					.addEntity(SubstructureAlert.class)
					.setParameter("smart", smart)
					.list();
			if (!alerts.isEmpty() && alerts.get(0).image != null)
			{
				alertId = alerts.get(0).id;
				href = "alerts/show.do?id=" + alerts.get(0).id;
				displayName = alerts.get(0).name;
				depictionPath = "alerts/image.do?id="+alerts.get(0).id+"&render-mode=popup&hasImage=true";
			}
		}
		else if (name.length() == 27 && name.matches("[A-Z]{14}-[A-Z]{10}-N"))
		{
			String inchi1 = name.substring(0, 14);
			// Look up a fragment by inchi key
			Long id = (Long) Globals.alternateSession().createCriteria(Fragment.class)
					.add(Restrictions.like("inchi1", inchi1))
					.setProjection(Projections.id())
					.uniqueResult();

			if (id != null)
				depictionPath = "depiction.jsp?frag_id=" + id;

		}
		else if (name.startsWith(QSPRConstants.SMILES_FORMAT + ":"))
		{
			int p = QSPRConstants.SMILES_FORMAT.length();
			int p2 = name.lastIndexOf(':');
			if (p != p2)
				depictionPath = "depiction.jsp?mol=" + name.substring(p + 1, p2);
			else
				depictionPath = "depiction.jsp?mol=" + name.substring(p + 1);
		}

	}

	/**
	 * An empty constructor to make JAXB happy
	 */
	public SetCompareResult()
	{

	}
}
