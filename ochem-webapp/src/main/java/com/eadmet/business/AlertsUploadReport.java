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

import java.util.ArrayList;
import java.util.List;

/**
 * A simple report of the alert upload
 * @author midnighter
 */
public class AlertsUploadReport
{
	public int dublicates = 0;
	public int errors = 0;
	public int successes = 0;
	public List<String> warnings = new ArrayList<String>();
	
	public String toString()
	{
		String msg = String
				.format("%s new alerts added\n%s errors skipped\n%s dublicates skipped\n\n",
						successes, errors, dublicates);
		for (String string : warnings)
			msg += string + "\n";
		
		
		return msg;
	}
}
