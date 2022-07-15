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

/*
 * Produced by RecordStubBuilder from UploadedRow (which is in turn produced by BatchUploadParser)
 * Should act as a pseudo-record, which holds all the values in an uploaded row in appropriate fields
 * 
 * The idea is to lose no information from this sheet on this step (yet)
 * Therefore "ColumnValue" instead of just "String" - as we want to keep 
 * the e.g column number information to show to the user in case of errors
 * 
 * Same goes for ValueList<ColumnValue> - there may be potentially e.g. two 
 * UNIT columns or two NAME columns for the same VALUE column.
 * The first case is OK, for the for the second one we want to show a warning.
 */
public class RecordStub 
{
	public PropertyStub property = new PropertyStub();

	ValueList<ColumnValue> molecules = new ValueList<ColumnValue>();
	public ValueList<ColumnValue> article = new ValueList<ColumnValue>();
	public ValueList<ColumnValue> basket = new ValueList<ColumnValue>();
	public ValueList<ColumnValue> name = new ValueList<ColumnValue>();
	public ValueList<ColumnValue> casrn = new ValueList<ColumnValue>();
	public ValueList<ColumnValue> comments = new ValueList<ColumnValue>();

	public ValueList<ColumnValue> hidden = new ValueList<ColumnValue>();
	public ValueList<ColumnValue> line = new ValueList<ColumnValue>();
	public ValueList<ColumnValue> page = new ValueList<ColumnValue>();
	public ValueList<ColumnValue> table = new ValueList<ColumnValue>();
	public ValueList<ColumnValue> n = new ValueList<ColumnValue>();
	public ValueList<ColumnValue> externalid = new ValueList<ColumnValue>();
	public ValueList<ColumnValue> evidence = new ValueList<ColumnValue>();


	public RecordStub()
	{

	}

	public RecordStub(RecordStub p) // Create a new and share all the "common" fields
	// Used in multiple record per row situations
	{
		super();
		article = p.article;
		basket = p.basket;
		name = p.name;
		casrn = p.casrn;
		comments = p.comments;
		molecules = p.molecules;
	}

	public static class NameValueStub
	{
		public ColumnValue name;
		public ValueList<ValueStub> value = new ValueList<ValueStub>();
	}

	public static class ValueStub
	{
		public ColumnValue value;
		public ValueList<ColumnValue> predicate = new ValueList<ColumnValue>();
		public ValueList<ColumnValue> accuracy = new ValueList<ColumnValue>();
		public ValueList<ColumnValue> interval = new ValueList<ColumnValue>();
		public ValueList<ColumnValue> unit = new ValueList<ColumnValue>();
	}

	public static class PropertyStub extends NameValueStub
	{
		public ValueList<NameValueStub> conditions = new ValueList<NameValueStub>();
	}

	public static class ColumnValue
	{
		public UploadedColumnSchema column;
		public String value;
		public ColumnValue(UploadedColumnSchema column, String value)
		{
			this.column = column;
			this.value = value;
		}

		public ColumnValue()
		{

		}

		public void setColumnValue(UploadedColumnSchema column, String value)
		{
			this.column = column;
			this.value = value;
		}

		public String toString()
		{
			return value+" ("+column+")";
		}
	}

	public static class ValueList<T> extends ArrayList<T>
	{
		private static final long serialVersionUID = 1L;

		public T getLast()
		{
			if (size() == 0)
				return null;
			else
				return get(size()-1);
		}

		public T getFirst()
		{
			if (size() == 0)
				return null;
			else
				return get(0);
		}
	}
}
