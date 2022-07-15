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

package com.eadmet.exceptions;

/**
 * An "expected" exception, reported to the user in a friendly manner rather than via a stacktrace
 * @author midnighter
 */
public class CriticalException extends UserFriendlyException 
{

	public static final String CRITICAL_ERROR = "ERROR: Critical ";

	private static final long serialVersionUID = 1L;

	public CriticalException(String st)
	{
		super(st != null && !st.toUpperCase().contains(CRITICAL_ERROR)?CRITICAL_ERROR+st:st);
	}

	public CriticalException(Exception e)
	{
		super(e.getMessage());
	}

	public CriticalException()
	{

	}
}
