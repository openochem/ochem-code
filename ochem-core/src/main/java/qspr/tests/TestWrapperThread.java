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

package qspr.tests;

import qspr.ThreadScope;
import qspr.entities.Session;
import qspr.util.WrapperThread;

@Deprecated
public abstract class TestWrapperThread extends WrapperThread 
{
	abstract protected void testWrapped();
	
	@Override
	public void wrapped()
	{
		ThreadScope.get().threadClassLoader = TestWrapperThread.class.getClassLoader();
		ThreadScope.get().disableTrackChanges = true;
		ThreadScope.get().userSession = Session.getFirstSession("test");
		
		testWrapped();
		
	}
}
