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

package qspr.business;

import java.util.HashSet;
import java.util.Set;


/**
 * A filter for batch upload preview browser
 */
public class BatchUpoadBrowserFilter
{
	public StatusFilter status;
	public Set<Integer> lines = new HashSet<Integer>();
	
	public BatchUpoadBrowserFilter(String status, String lines)
	{
		try {
			if (status != null)
				this.status = StatusFilter.valueOf(status);
			else
				this.status = StatusFilter.all;
		} catch (Exception e)
		{
			this.status = StatusFilter.all;
		}
		if (lines == null)
			return;
		String[] pcs = lines.split(",");
		for (String p : pcs)
		{
			try {
				this.lines.add(Integer.valueOf(p.trim()));
			} catch (Exception e)
			{
				e.printStackTrace();
			}
		}
	}
	
	public enum StatusFilter {all,valid,invalid,error,fatal_error,warning,duplicate_internal,duplicate_external};
}
