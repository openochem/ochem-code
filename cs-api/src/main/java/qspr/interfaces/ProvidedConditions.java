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

package qspr.interfaces;

import java.io.Serializable;
import java.util.List;

import qspr.metaserver.util.ShortCondition;

public interface ProvidedConditions extends Serializable{

	public static final String SOLVENT = "solvent";
	public static final double MISSED_VALUE = -1.;

	public List<ShortCondition> getConditions();

	boolean hasConditions();

	static public String isReplaceableCondition(String type) {
		if(type.contains(":")) {
			String parts[] = type.split(":");
			if(parts.length>1)
				type = parts[1];
		}
		return type.toLowerCase().startsWith(SOLVENT)? SOLVENT : null ;
	}

}
