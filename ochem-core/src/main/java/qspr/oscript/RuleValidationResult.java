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

package qspr.oscript;

import java.util.HashSet;
import java.util.Set;

import qspr.util.SmartText;

public class RuleValidationResult
{
	Set<String> unknownProperties = new HashSet<String>();
	Set<String> unknownConditions = new HashSet<String>();
	Set<String> criticalErrors = new HashSet<String>();
	
	public boolean isValid()
	{
		return unknownConditions.isEmpty() && unknownProperties.isEmpty() && criticalErrors.isEmpty();
	}
	
	public String toString()
	{
		String st = "";
		if (!unknownProperties.isEmpty())
			st += "Unrecognized property names: " + unknownProperties + "\n";
		if (!unknownConditions.isEmpty())
			st += "Unrecognized condition names: " + unknownConditions + "\n";
		
		if (!criticalErrors.isEmpty())
			st += "Critical errors: " + criticalErrors + "\n";
		
		return st;
	}
	
	public void critical(String message)
	{
		criticalErrors.add(SmartText.transform(message));
	}
}