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

package qspr.dao;

import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement(name = "standardizer")
public enum ChemInfEngine {
	CDK,
	CHEMAXON,
	NONE
	;

	public static ChemInfEngine getEnum(String s){
		if(CDK.name().equalsIgnoreCase(s)){
			return CDK;
		}else if(CHEMAXON.name().equalsIgnoreCase(s))
			return CHEMAXON;
		else if(NONE.name().equalsIgnoreCase(s))
			return NONE;
		throw new IllegalArgumentException("No Enum specified for this string " + s);
	}

};