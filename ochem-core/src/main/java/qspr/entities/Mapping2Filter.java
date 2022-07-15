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

package qspr.entities;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.FetchType;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlTransient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.Globals;

/**
 * A database entity to store temporary filters for molecules
// Its stored in an in-memory table for speed and for automatic expiery
 * 
 * Feb 16, 2012
 * @author midnighter
 *
 */
@Entity
public class Mapping2Filter 
{
	@Id
	@Column(name = "mf_id")
	@javax.persistence.OrderBy
	@XmlAttribute
	public Integer id;
	
	@Column(name = "filter_id")
	@XmlTransient
	public Integer filterId;
	
	@ManyToOne (fetch=FetchType.LAZY)
	@JoinColumn(name="mapping2_id", referencedColumnName="mapping2_id")
	public Mapping2 mapping2;
	
	/**
	 * Generate a new random filter identifier
	 * @return
	 */
	public static int generateFilterID()
	{
		return Long.valueOf(Math.round(Math.random() * 1000000)).intValue();
	}
	
	public Mapping2Filter()
	{
		
	}
	
	/**
	 * Add the given compound to the filter.
	 * For speed, you can use direct queries:
	 * 
	 * Example:
	 * insert into Mapping2Filter(filter_id, mapping2_id) select ... from Mapping2 where ...
	 * 
	 * @param filterID
	 * @param compoundID
	 */
	public static void addCompoundToFilter(int filterID, int compoundID)
	{
		Globals.session().createSQLQuery("insert into Mapping2Filter(filter_id, mapping2_id) values(:f, :m)").setInteger("f", filterID).setInteger("m", compoundID).executeUpdate();
	}
	
	public static void clear(int filterID)
	{
		logger.info("Clearing filter " + filterID);
		Globals.session().createSQLQuery("delete from Mapping2Filter where filter_id = :f").setInteger("f", filterID).executeUpdate();	
	}
	
	public static void addBasketCompoundsToFilter(int filterID, long basketId)
	{
		Globals.session().createSQLQuery("insert into Mapping2Filter(filter_id, mapping2_id) select distinct :f, m.mapping2_id from Basket b join BasketEntry be using (basket_id) join ExperimentalProperty ep using (exp_property_id) natural join Molecule m where b.basket_id = :b").setInteger("f", filterID).setLong("b", basketId).executeUpdate();
	}
	
	public static void addCompoundToFilter(int filterID, Collection<Integer> compoundIDs) {
		List<Integer> ids = new ArrayList<Integer>();
		ids.addAll(compoundIDs);
		addCompoundToFilter(filterID, ids);
	}
	
	public static void addCompoundToFilter(int filterID, List<Integer> compoundIDs)
	{
		int batchSize = 1000;
		
		// Insert the first batch
		StringBuffer sb = new StringBuffer();
		int max = Math.min(compoundIDs.size(), batchSize);
		for (int i = 0; i < max; i++)
		{
			sb.append("("+filterID + "," + compoundIDs.get(i) + ")");
			if (i < max - 1)
				sb.append(",");
		}
		
		logger.info(String.format("Adding %d compounds to filter %d", max, filterID));
		Globals.session().createSQLQuery("insert into Mapping2Filter(filter_id, mapping2_id) values " + sb.toString()).executeUpdate();
		
		// Do the next batch
		if (compoundIDs.size() > batchSize)
			addCompoundToFilter(filterID, compoundIDs.subList(batchSize, compoundIDs.size()));
	}
	
	private static final Logger logger = LogManager.getLogger(Mapping2Filter.class);
}
