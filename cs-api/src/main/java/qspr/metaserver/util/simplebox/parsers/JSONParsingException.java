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

package qspr.metaserver.util.simplebox.parsers;

public class JSONParsingException extends Exception {

	private static final long serialVersionUID = -1393364105190162572L;

	public JSONParsingException() {
		super();
	}

	public JSONParsingException(String message, Throwable cause) {
		super(message, cause);
	}

	public JSONParsingException(String message) {
		super(message);
	}

	public JSONParsingException(Throwable cause) {
		super(cause);
	}

}
