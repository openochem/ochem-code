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

package qspr.workflow.structure;

import java.io.Serializable;

import javax.xml.bind.annotation.XmlAttribute;

public class Connection implements Serializable
{
	private static final long serialVersionUID = 1L;

	@XmlAttribute
	public String from = "in";
	
	@XmlAttribute
	public String to = "out";
	
	@XmlAttribute
	public Integer fromPort = 0;
	
	@XmlAttribute
	public Integer toPort = 0;
	
	public Connection()
	{
		
	}
	
	public Connection(String from, String to)
	{
		this.from = from;
		this.to = to;
	}
}
