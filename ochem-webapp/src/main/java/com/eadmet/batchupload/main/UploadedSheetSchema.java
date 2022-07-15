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

package com.eadmet.batchupload.main;

import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import com.eadmet.batchupload.main.BatchUploadMessage.BatchUploadMessageType;
import com.eadmet.batchupload.main.UploadedColumnSchema.UploadedColumnType;

@XmlRootElement
public class UploadedSheetSchema 
{
	@XmlAttribute
	public Boolean ignore = true;

	@XmlAttribute
	public Boolean addDummyProperty = false;

	@XmlAttribute
	public Boolean dummyPropertyAdded = false;
	
	@XmlAttribute
	public Boolean missingProperty = false;
	
	@XmlAttribute
	public String name;
	
	@XmlAttribute
	private String originalName;
	
	@XmlAttribute
	private Integer number;
	
	@XmlElement
	public List<BatchUploadMessage> messages = new ArrayList<BatchUploadMessage>();
	
	@XmlTransient
	private List<UploadedColumnSchema> columns = new ArrayList<UploadedColumnSchema>();
	
	public String toString()
	{
		return number+":"+name+":"+originalName+":"+ignore+":"+addDummyProperty+":"+columns;
	}
	
	public Integer getNumber()
	{
		return number;
	}
	
	@XmlAttribute
	public BatchUploadMessageType getSheetStatus()
	{
		BatchUploadMessageType type = BatchUploadMessageType.notice;
		for (BatchUploadMessage m : messages)
			if (type == null || m.type.ordinal() < type.ordinal())
				type = m.type;
		return type;
	}
	
	@XmlElement
	public List<UploadedColumnSchema> getColumns()
	{
		List<UploadedColumnSchema> columns = new ArrayList<UploadedColumnSchema>();
		columns.addAll(this.columns);
		return columns;
	}
	
	public void remapColumnToKnown(Integer columnIndex, UploadedColumnType type)
	{
		remapColumnToKnown(getColumns().get(columnIndex), type);
	}
	
	public void remapColumnToKnown(UploadedColumnSchema column, UploadedColumnType type)
	{
		column.name = type.name();
		column.ignore = false;
	}
	
	public void remapColumnToProperty(Integer columnIndex, String propertyName)
	{
		remapColumnToProperty(getColumns().get(columnIndex), propertyName);
	}
	
	public void remapColumnToProperty(UploadedColumnSchema column, String propertyName)
	{
		column.name = propertyName;
		column.ignore = false;
	}
	
	protected void setColumns(List<UploadedColumnSchema> columns)
	{
		this.columns = columns;
	}
	
	public UploadedColumnSchema addColumn(UploadedColumnSchema c)
	{
		return addColumn(columns.size(), c);
	}
	
	public UploadedColumnSchema addColumn(Integer position, UploadedColumnSchema c)
	{
		columns.add(position, c);
		for (int i=0; i<columns.size(); i++)
			columns.get(i).remappedNumber = i;
		return c;
	}
	
	public UploadedColumnSchema addColumn(Integer position, String name, String defaultValue)
	{
		UploadedColumnSchema c = new UploadedColumnSchema();
		c.type = UploadedColumnType.undefined;
		c.name = name;
		c.defaultValue = defaultValue;
		return addColumn(position, c);
	}
	
	public UploadedColumnSchema addColumn(UploadedColumnType name, String defaultValue)
	{
		return addColumn(columns.size(), name.toString(), defaultValue);
	}
	
	public UploadedColumnSchema addColumn(String name, String defaultValue)
	{
		return addColumn(columns.size(), name, defaultValue);
	}
	
	public void removeColumn(int position)
	{
		columns.remove(position);
	}
	
//	public UploadedColumn addColumnAfter(UploadedColumn column, String name, String defaultValue)
//	{
//		int index = -1;
//		for (int i=0; i<columns.size(); i++)
//			if (columns.get(i).equals(column))
//			{
//				index = i;
//				break;
//			}
//		return addColumn(index, type, name, defaultValue);
//	}
//	
//	public UploadedColumn addColumnAfter(UploadedColumn column, String name, List<String> rowValues)
//	{
//		UploadedColumn c = addColumnAfter(column, type, name, (String) null);
//		c.rowValues = rowValues;
//		return c;
//	}
//
//	
//	public UploadedColumn addColumn(Integer position, String name, List<String> rowValues)
//	{
//		UploadedColumn c = addColumn(position, type, name, (String) null);
//		c.rowValues = rowValues;
//		return c;		
//	}
//	

//	
//	public UploadedColumn addColumn(String name, List<String> rowValues)
//	{
//		return addColumn(columns.size()-1, type, name, rowValues);	
//	}
//
//	
//	public List<UploadedColumn> getColumns(UploadedColumnType type)
//	{
//		List<UploadedColumn> result = new ArrayList<UploadedColumn>();
//		for (UploadedColumn column : columns)
//			if (column.type == type)
//				result.add(column);
//		return result;
//	}

	public UploadedSheetSchema()
	{
	}
	
	public UploadedSheetSchema(String sheetName, Integer sheetNumber, List<UploadedColumnSchema> columns)
	{
		this.name = this.originalName = sheetName;
		this.number = sheetNumber;
		this.columns = columns;
	}
}
