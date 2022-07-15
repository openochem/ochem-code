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

package com.eadmet.mmpa.domain;

import java.io.Serializable;
import java.util.List;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.Id;

import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;
import org.hibernate.exception.ConstraintViolationException;

import qspr.Globals;

@Entity
public class MMPIndex implements Serializable
{
	private static final long serialVersionUID = 1L;

	@Id
	@Column(name = "scaffold_id")
	public String scaffoldId;
	
	@Id
	@Column(name = "fragment_id")
	public long fragmentId;
	
	@Id
	@Column(name = "mapping2_id")
	public int mapping2Id;
		
	public static MMPIndex saveIndex(String scaffoldId, long fragmentId, int mapping2Id) {
		MMPIndex entry = new MMPIndex();
		entry.scaffoldId = scaffoldId;
		entry.fragmentId = fragmentId;
		entry.mapping2Id = mapping2Id;
		
		try
		{
			if (entryExists(scaffoldId, fragmentId, mapping2Id))
				return null;
			Globals.session().save(entry);
			return entry;
		} catch (ConstraintViolationException e) {
			return null;
		}
	}
	
	private static boolean entryExists(String scaffoldId, long fragmentId, int mapping2Id) {
		return ((Long)Globals.session().createCriteria(MMPIndex.class)
			.add(Restrictions.eq("scaffoldId", scaffoldId))
			.add(Restrictions.eq("fragmentId", fragmentId))
			.add(Restrictions.eq("mapping2Id", mapping2Id))
			.setProjection(Projections.rowCount())
			.uniqueResult()) > 0;
	}
	
	@SuppressWarnings("unchecked")
	public static List<MMPIndex> getEntries(String scaffold) {
		return Globals.session().createCriteria(MMPIndex.class).add(Restrictions.eq("scaffoldId", scaffold)).list();
	}

	public static int newScaffolds;	
}
