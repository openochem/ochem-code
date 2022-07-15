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

import org.hibernate.event.spi.PreUpdateEvent;
import org.hibernate.event.spi.PreUpdateEventListener;

import qspr.ThreadScope;
import qspr.entities.Basket;

import com.eadmet.exceptions.UserFriendlyException;


// Should be (mostly) a safety measure... normally we should block creation of data in the UI for guest users
public class AnonymousUpdateListener implements PreUpdateEventListener  
{
	private static final long serialVersionUID = 1L;
	public boolean onPreUpdate(PreUpdateEvent event) 
	{
		if (Boolean.TRUE.equals(ThreadScope.get().disableTrackChanges))
        	return false;
		
        Object entity = event.getEntity();
    	if (entity instanceof UserContributedEntity)
   		if (((UserContributedEntity)entity).getIntroducer() == null && !(entity instanceof Basket))
   			throw new UserFriendlyException("Guest users are not allowed to create data in OCHEM. Please, register to upload your data.");
    	
		return false;
    }
}
