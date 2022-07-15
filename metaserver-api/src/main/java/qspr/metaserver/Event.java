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

package qspr.metaserver;

import java.util.ArrayList;
import java.util.List;

// Simple event target
// Midnighter

public class Event<T>
{
	public Object parent;
	private List<EventListener<T>> listeners = new ArrayList<EventListener<T>>();
	
	public Event(Object parent)
	{
		this.parent = parent;
	}
	
	public void fire()
	{
		for (EventListener<T> listener : listeners) {
			listener.onEvent(this, null);
		}
	}
	
	public void fire(T arg)
	{
		for (EventListener<T> listener : listeners) {
			listener.onEvent(this, arg);
		}
	}
	
	public void addListener(EventListener<T> listener)
	{
		if (!listeners.contains(listener))
			listeners.add(listener);
	}
}
