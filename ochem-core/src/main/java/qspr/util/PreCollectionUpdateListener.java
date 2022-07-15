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

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.engine.spi.SessionFactoryImplementor;
import org.hibernate.event.spi.PreCollectionUpdateEvent;
import org.hibernate.event.spi.PreCollectionUpdateEventListener;
import org.hibernate.persister.collection.CollectionPersister;

import qspr.Globals;
import qspr.ThreadScope;
import qspr.annotations.Loggable;
import qspr.entities.Action;


public class PreCollectionUpdateListener implements PreCollectionUpdateEventListener 
{
	private static final long serialVersionUID = 1L;
	private static transient final Logger logger = LogManager.getLogger(PreCollectionUpdateListener.class);

	@SuppressWarnings({ "unchecked", "rawtypes" })
	public void onPreUpdateCollection(PreCollectionUpdateEvent event) 
	{
		try {

			if (Boolean.TRUE.equals(ThreadScope.get().disableTrackChanges))
				return;

			logger.debug("[Listener] Collection "+event.getCollection().getRole()+" is being updated.");
			Object entity = event.getCollection().getOwner();
			Class eClass = entity.getClass();
			Loggable classLoggable = (Loggable) eClass.getAnnotation(Loggable.class);

			if (classLoggable == null)
				return;

			String[] elements = event.getCollection().getRole().split("\\.");
			String propertyName = elements[elements.length - 1];

			Action action = new Action(entity);
			String name = action.getPropertyNameInLog(propertyName);
			if (name == null)
				return;

			CollectionPersister persister = ((SessionFactoryImplementor)event.getSession().getSessionFactory()).getCollectionPersister(event.getCollection().getRole());

			Collection collectionBefore = getCollection(event.getCollection().getStoredSnapshot());
			Collection collectionAfter = getCollection(event.getCollection().getSnapshot(persister));
			List tmp = new ArrayList(collectionBefore);
			tmp.removeAll(collectionAfter);
			if (tmp.size() > 0)
				action.comment += "Deleted "+name+": "+tmp+"\n";
			tmp = new ArrayList(collectionAfter);
			tmp.removeAll(collectionBefore);
			logger.info(tmp);
			if (tmp.size() > 0)
				action.comment += "Added "+name+": "+tmp;

			if (!action.comment.equals(""))
			{
				Globals.session().save(action);
				action.publish(event.getCollection().getOwner());
			}

		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	}

	@SuppressWarnings("rawtypes")
	private Collection getCollection(Object object)
	{
		if (object instanceof HashMap)
			return ((HashMap) object).values();
		else
			return (Collection) object;
	}

}
