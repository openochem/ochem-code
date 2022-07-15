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

package qspr.util;

import org.hibernate.event.spi.PostInsertEvent;
import org.hibernate.event.spi.PostInsertEventListener;
import org.hibernate.persister.entity.EntityPersister;

import qspr.Globals;
import qspr.annotations.Loggable;
import qspr.entities.Action;
import qspr.entities.Action.ActionType;


public class InsertListener implements PostInsertEventListener 
{
	private static final long serialVersionUID = 1L;

	@SuppressWarnings({ "rawtypes", "unchecked" })
	public void onPostInsert(PostInsertEvent event) 
	{
		try
		{
			Object entity = event.getEntity();
			Class eClass = entity.getClass();
			//logger.info("onPreInsert: "+eClass.getName());
			Loggable classLoggable = (Loggable) eClass.getAnnotation(Loggable.class);
			if (classLoggable != null)
			{
				Action action = new Action(entity);
				action.comment = "Item has been created";
				action.type = ActionType.CREATE;
				Globals.session().save(action);
				action.publish(entity);
			}
		} catch (Exception e)
		{
			e.printStackTrace();
		}
	}

	@Override
	public boolean requiresPostCommitHanding(EntityPersister arg0) 
	{
		return false;
	}

}
