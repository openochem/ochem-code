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

package com.eadmet.mmpa.domain;

import java.util.Collections;
import java.util.List;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlRootElement;

import org.apache.commons.lang.StringUtils;

@XmlRootElement
public class MergedAnnotationSet extends MMPAnnotationSet
{
	List<Long> ids;
	
	public MergedAnnotationSet()
	{
		
	}
	
	public MergedAnnotationSet(List<Long> ids)
	{
		Collections.sort(ids);
		this.ids = ids;
		this.name = "Combined set";
	}
	
	@Override
	public List<Long> getIdList()
	{
		return ids;
	}
	
	@XmlAttribute
	public String getId()
	{
		return StringUtils.join(ids, ".");
	}
}
