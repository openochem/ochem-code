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

package qspr.frontend;

import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlElements;
import javax.xml.bind.annotation.XmlRootElement;

import qspr.entities.Alert;
import qspr.entities.XlsColumnHeader;

@XmlRootElement(name = "list")
@SuppressWarnings("rawtypes")
public class MarshalableList
{
	@XmlAttribute
	public String name;
	
	public MarshalableList()
	{
		list = new ArrayList();
	}
	
	public MarshalableList(List arg)
	{
		list = arg;
	}
	
	public MarshalableList setName(String name)
	{
		this.name = name;
		return this;
	}

	@XmlElements
	({
		@XmlElement(name = "integer", type=Integer.class),
		@XmlElement(name = "string", type=String.class),
		@XmlElement(name = "message", type=Alert.class),
		@XmlElement(name = "xlscolumnheader", type=XlsColumnHeader.class),
	    @XmlElement
	})
	
	public List list;
}
