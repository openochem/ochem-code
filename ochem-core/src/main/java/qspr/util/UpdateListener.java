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
import java.util.Iterator;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.engine.spi.CollectionEntry;
import org.hibernate.event.spi.PreUpdateEvent;
import org.hibernate.event.spi.PreUpdateEventListener;

import qspr.Globals;
import qspr.ThreadScope;
import qspr.annotations.Loggable;
import qspr.entities.Action;
import cern.colt.Timer;


public class UpdateListener implements PreUpdateEventListener {

	private static final long serialVersionUID = 1L;
	private static transient final Logger logger = LogManager.getLogger(UpdateListener.class);

	@SuppressWarnings({ "rawtypes", "unchecked" })
	public boolean onPreUpdate(PreUpdateEvent event) 
	{
		try
		{
			if (Boolean.TRUE.equals(ThreadScope.get().disableTrackChanges))
				return false;

			Object entity = event.getEntity();
			Class eClass = entity.getClass();


			Loggable classLoggable = (Loggable) eClass.getAnnotation(Loggable.class);
			if (classLoggable != null)
			{
				logger.info("[Listener] Entity "+eClass.getName()+" being updated");
				Action action = new Action(entity);
				boolean modified = false;

				Object[] oldValues = event.getOldState();
				Object[] newValues = event.getState();
				String[] properties = event.getPersister().getPropertyNames();

				// Workaround on bug HHH-2763
				// get CollectionEntries before traversing the relationship 
				Timer list = new Timer(); list.start();
				List entriesBeforeLoad = new ArrayList();
				entriesBeforeLoad.addAll(event.getSource().getPersistenceContext().getCollectionEntries().values());
				list .stop();

				// Iterate through all fields of the updated object
				String log = "";
				for (int i = 0; i < properties.length; i++)
				{
					String name;
					if ((name = action.getPropertyNameInLog(properties[i])) != null)
					{
						if (!(oldValues[i] instanceof Collection))
							if (newValues[i] != null)
							{
								if (newValues[i] instanceof ChangesTracker && oldValues[i] != null)
								{
									String changes = ((ChangesTracker)oldValues[i]).getChanges((ChangesTracker)newValues[i]);
									if (changes != null)
									{
										modified = true;
										log += name+":"+changes+"\n";
									}
								} else
									if (!newValues[i].equals(oldValues[i]))
									{
										modified = true;
										log += name+": "+oldValues[i]+" to "+newValues[i]+"\n";
									}
							}
							else
							{
								if (oldValues[i] != null)
									log += "Field '"+name+"'"+" has been cleared\n";
							}
					}
				}

				// Workaround on hibernate bug HHH-2763
				// get CollectionEntries after traversing the relationship
				List entriesAfterLoad = new ArrayList();
				entriesAfterLoad.addAll(event.getSource().getPersistenceContext().getCollectionEntries().values());

				// remove only if there was a new entity loaded
				if (entriesBeforeLoad.size() != entriesAfterLoad.size())
				{
					entriesAfterLoad.removeAll(entriesBeforeLoad);
					Iterator iter = entriesAfterLoad.iterator();
					while (iter.hasNext())
					{
						CollectionEntry ce = (CollectionEntry) iter.next();
						ce.setProcessed(true);
					}
				}

				action.comment = log;
				if (modified)
				{
					Globals.session().save(action);
					action.publish(event); //Additional logic for logging
				}

			}
		} catch (Exception e)
		{
			e.printStackTrace();
		}
		return false;
	}

}
