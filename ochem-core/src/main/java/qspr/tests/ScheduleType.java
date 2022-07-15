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

/**
 * Scheduling options for the periodic tests
 * @author midnighter
 *
 */
public enum ScheduleType {
	MINUTE,
	HOURLY,
	BIHOURLY,
	TWELVE_HOURS,
	SIX_HOURS,
	DAILY,
	
	AT_FIVE_AM,
	
	STARTUP;
	
	public int getMilliseconds()
	{
		switch (this)
		{
		case MINUTE:
			return 1000 * 60;
		case HOURLY:
			return 1000 * 60 * 60;
		case TWELVE_HOURS:
			return 1000 * 60 * 60;
		case BIHOURLY:
			return 1000 * 60 * 60 * 2;
		case DAILY:
			return 1000 * 60 * 60 * 24;
		case SIX_HOURS:
			return 1000 * 60 * 60 * 6;
		case AT_FIVE_AM:
			return 1000 * 60 * 60 * 5;
		case STARTUP:
			return -1;
		}
		return -1;
	}
	
}
