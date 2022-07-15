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

package com.eadmet.business;

import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.Globals;
import qspr.entities.Tag;

public class TagsService
{
	@SuppressWarnings("unchecked")
	public static void addTagByMP2(Tag tag, List<Integer> mp2Ids)
	{
		// Check if we already have some molecules in the tag
		List<Integer> existentMp1Ids = Globals.session().createSQLQuery("select mapping1_id from MoleculeTag where tag_id=" + tag.id).list();

		if (existentMp1Ids.isEmpty())
			Globals.session().createSQLQuery("insert into MoleculeTag(tag_id, mapping1_id) select "+tag.id+", mapping1_id from Mapping2 where mapping2_id in (:ids) group by mapping1_id")
			.setParameterList("ids", mp2Ids).executeUpdate();
		else
		{
			logger.info("Tag " + tag.name + " already has " + existentMp1Ids.size() + " entries");
			Globals.session().createSQLQuery("insert into MoleculeTag(tag_id, mapping1_id) select "+tag.id+", mapping1_id from Mapping2 where mapping2_id in (:ids) and not mapping1_id in (:existent_ids) group by mapping1_id")
			.setParameterList("ids", mp2Ids)
			.setParameterList("existent_ids", existentMp1Ids)
			.executeUpdate();
		}
	}

	private static final Logger logger = LogManager.getLogger(TagsService.class);
}
