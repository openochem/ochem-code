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

package qspr.workflow.datatypes;

import static java.lang.String.format;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Unmarshaller;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlElementWrapper;
import javax.xml.bind.annotation.XmlElements;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import qspr.dao.ChemInfEngine;
import qspr.dao.Various;
import qspr.metaserver.cs.WorkflowNodeServer;
import qspr.metaserver.protocol.DataSize;
import qspr.util.ClassCompressor;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.NumericalValueStandardizer;

@XmlRootElement(name = "datatable")
public class DataTable implements Serializable, Cloneable, Iterable<AbstractDataRow>
{
	private static final long serialVersionUID = 2L;
	public static final String DEFAULT_LABEL = "default";

	@XmlAttribute
	protected boolean compactRowFormat = false;

	@XmlAttribute
	public boolean compressStrings = false;

	@XmlTransient
	public Serializable attachment;

	/**
	 * Required to indicate which columns are conditions
	 */
	@XmlTransient
	public Map<String, Map<String, Serializable>> columnAttachments = new HashMap<String, Map<String, Serializable>>(); //Column Name -> Map Of Attachments

	@XmlTransient
	public List<Integer> storedUncompressedToCompressedRows;

	@XmlTransient
	public List<Map<String, Serializable>> storedAttachements;

	@XmlTransient
	public int currentRow = 0;

	@XmlAttribute
	public String id;

	@XmlElementWrapper(name = "columns")
	@XmlElement(name = "column")
	public List<String> columns = new ArrayList<String>();

	//	@XmlTransient
	//	This was supposed to be a map to help avoid the columns.indexOf call, which is slow and in moderately big tables takes noticable time
	//	But everywhere the columns list is used directly and uncontrollably, so keeping the map consistent is impossible without changing lotsa stuff
	//  TODO novserj: do it some day 
	//	protected Map<String,Integer> columnMap = new HashMap<String, Integer>();

	@XmlElementWrapper(name = "rows")
	@XmlElements({ @XmlElement(name = "row", type = DataRow.class), @XmlElement(name = "crow", type = CompactDataRow.class),
		@XmlElement(name = "srow", type = SparseDataRow.class) })
	private List<AbstractDataRow> data = new ArrayList<AbstractDataRow>();

	public static DataTable fromXml(JAXBContext jContext, String filename) throws Exception
	{
		Unmarshaller ml = jContext.createUnmarshaller();
		return (DataTable) ml.unmarshal(new File(filename));
	}

	public void setAttachment(Serializable attachment)
	{
		this.attachment = attachment;
	}

	public DataTable(Serializable val)
	{
		try
		{
			addColumn(DEFAULT_LABEL);
			addRow();
			setValue(val);
		} catch (Exception e)
		{

		}
	}

	public DataTable()
	{

	}

	public DataTable compress(){

		DataTable newTab  = getEmptyCopy();
		newTab.storedUncompressedToCompressedRows = new ArrayList<Integer>();
		newTab.storedAttachements = new ArrayList<Map<String, Serializable>>(); 
		Map<Integer,Integer> hashes = new HashMap<Integer,Integer>();

		for(int i = 0; i<data.size();i++){
			int hash = data.get(i).getHashValues();
			if(!hashes.containsKey(hash)) {
				hashes.put(hash, newTab.data.size());
				newTab.data.add(data.get(i));
			}
			newTab.storedUncompressedToCompressedRows.add(hashes.get(hash));
			newTab.storedAttachements.add(data.get(i).attachments); // all stored, usually unique per row
		}

		if(newTab.data.size() < data.size()) {
			System.out.println("Compressed: " + data.size() + " " + newTab.data.size());
			return newTab;
		}

		return this;

	}

	public DataTable deCompress(){
		if(storedUncompressedToCompressedRows != null) {
			List<AbstractDataRow> newData = new ArrayList<AbstractDataRow>();
			for(int i = 0; i<storedUncompressedToCompressedRows.size();i++) {
				AbstractDataRow r = data.get(storedUncompressedToCompressedRows.get(i)).getDeeperCopy();
				r.attachments = storedAttachements.get(i);
				newData.add(r);
			}
			data = newData;
			storedUncompressedToCompressedRows = null;
			storedAttachements = null;
		}
		return this;
	}

	public DataTable(boolean compactRowFormat)
	{
		this.compactRowFormat = compactRowFormat;
	}

	/**
	 * Convenience constructor from numeric array.
	 * 
	 * Usage example: new DataTable(new float[][]{{a,b,...},{c,d,...},...}, new String[]{"a","b",...}).
	 *  
	 * @param data Numeric table as content.
	 * @param columnNames Names of columns.
	 */
	public DataTable(final double[][] data, final String[] columnNames)
	{
		for (final String name : columnNames)
			addColumn(name);
		for (final double[] row : data)
		{
			addRow();
			for (int col = 0; col < row.length; ++col)
				setValue(col, row[col]);
		}
	}

	/**
	 * Convenience constructor from a datatable from an arrays.
	 * 
	 * Usage example: new DataTable(new String[]{"abc","xzy",...}, "name").
	 * 
	 * @param data String array as content.
	 * @param columnName Name of the only column.
	 */
	public DataTable(final Serializable[] data, final String columnName)
	{
		addColumn(columnName);
		for (final Serializable entry : data)
		{
			addRow();
			setValue(0, entry);
		}
	}

	public boolean isCompactRowFormat()
	{
		return compactRowFormat;
	}

	/**
	 * Return number of rows in the Table
	 * @return
	 */

	public int getRowsSize()
	{
		return data == null ? 0 : data.size();
	}

	/**
	 * Return number of rows without errors in the Table
	 * @return
	 */

	public int getRowsNoErrorsSize()
	{
		if(data == null) return 0;
		int n=0;
		for(int i=0;i<data.size();i++)
			if(!data.get(i).isError())n++;
		return n;
	}

	/**
	 * @return returns number of columns in the Table
	 */
	public int getColumnsSize()
	{
		return columns.size();
	}

	@XmlTransient
	private static AbstractDataRow stubRow;

	public AbstractDataRow getStubRow()
	{
		if (stubRow == null)
			stubRow = createRow(0).setStatus(WorkflowNodeServer.STUBROW);
		return stubRow;
	}

	private AbstractDataRow createRow(int size)
	{
		if (compactRowFormat)
			return new SparseDataRow();
		else
		{
			DataRow dr = new DataRow(size);
			dr.compressStrings = this.compressStrings;
			return dr;
		}
	}

	/**
	 * Creates row which has type supported by DataTable but does not store it in the table
	 * @return
	 */

	public AbstractDataRow createRow()
	{
		return createRow(columns.size());
	}

	public AbstractDataRow addStubRow()
	{
		return addRow(getStubRow());
	}

	public AbstractDataRow addRow()
	{
		return addRow(createRow());
	}

	public void deleteRow(int row) throws IOException
	{
		if (data.size() <= row)
			throw new IOException("Row index out of bounds.");
		data.remove((int) row);
		if (currentRow >= data.size())
			currentRow = data.size() - 1;
	}

	public void insertRow(int index, AbstractDataRow row)
	{
		data.add(index, row);
	}

	public void compact()
	{
		for (AbstractDataRow adr : data)
			adr.compact();
	}

	public void addColumn(String name)
	{
		if (getColumnIndex(name) != -1)
			throw new UserFriendlyException("Duplicated descriptors are not allowed."
					+ " This means that you have same descriptor twice: (" + name + ") at " + getColumnIndex(name) + " and at " + columns.size());

		columns.add(name);
	}

	public void addColumns(int numOfColumns)
	{
		for (int i = 0; i < numOfColumns; i++)
			addColumn(Integer.valueOf(i).toString());
	}

	public void setValue(int row, String col, Serializable value)
	{
		int i = columns.indexOf(col);
		if (i == -1)
		{
			// If the column is not there - create it automatically
			addColumn(col);
			i = columns.size() - 1;
		}
		setValue(row, i, value);
	}

	public void setValue(int row, int col, Serializable value)
	{
		if (col >= columns.size())
			throw new RuntimeException("DataTable error: cannot set column " + col + " with total columns num " + columns.size());
		AbstractDataRow drow = data.get(row);
		drow.setValue(col, value);
	}

	public void setValue(String col, Serializable value)
	{
		setValue(currentRow, col, value);
	}

	public void setValue(int col, Serializable value)
	{
		setValue(currentRow, col, value);
	}

	public void setValue(Serializable value)
	{
		setValue(0, value);
	}

	public Serializable getValue(int row, String col)
	{
		return getValue(row, getColumnIndex(col));
	}

	public int getColumnIndex(String col)
	{
		return columns.indexOf(col);
	}

	public Serializable getValue(int row, int col)
	{
		if (col >= columns.size())
			throw new RuntimeException("Datatable column number out of bounds. Index: " + col + ", columns" + columns.size());

		AbstractDataRow drow = data.get(row);

		return drow.getValue(col);
	}

	public AbstractDataRow getCurrentRow()
	{
		return data.get(currentRow);
	}

	public void setCurrentRow(int n)
	{
		currentRow = n;
	}

	public void setSize(int n) throws IOException
	{
		while (data.size() > n)
			deleteRow(n);
		currentRow = n - 1;
	}

	public AbstractDataRow getRow(int index)
	{
		return data.get(index);
	}

	/**
	 * Adding rows with the same columns (ASSUMED -- dangerous operation)!
	 * The type of columns is not checked and it is assumed to be the same as in the DataTable
	 * @param row
	 * @return
	 */

	public AbstractDataRow addRow(AbstractDataRow row)
	{
		if (compactRowFormat != row instanceof CompactDataRow && row.status != WorkflowNodeServer.STUBROW) // we do not care about the STUBROW
			throw new UserFriendlyException("Trying to combine compact and non-compact data rows");
		if (row.size() > columns.size())
			throw new UserFriendlyException("Trying to add row with " + row.size() + " columns in table with " + columns.size() + " columns");
		data.add(row);
		currentRow = data.size() - 1;
		return row;
	}

	public DataTable mergeColumnsWith(DataTable source) throws IOException
	{
		return mergeColumnsWith(source, "");
	}

	/**
	 * Merge two types of descriptors for the same rows
	 * It is assumed that descriptor names are different ones
	 * @param source
	 * @param suffix
	 * @return
	 * @throws IOException 
	 */

	public DataTable mergeColumnsWith(DataTable source, String suffix) throws IOException
	{
		if (source.getRowsSize() != getRowsSize())
			throw new RuntimeException("Merging datatables with unequal numbers of rows");

		int initialSize = getColumnsSize(); // number of columns in the initial table before the merge

		for (String column : source.columns)
		{
			if(columns.contains(column + suffix)) throw new IOException("The DataTable already contains this column: " + column + suffix);
			addColumn(column + suffix);
			columnAttachments.put(column + suffix, source.columnAttachments.get(column));
		}

		//columns.addAll(source.columns);
		reset();
		source.reset();
		while (nextRow())
		{
			source.nextRow();
			getCurrentRow().addColumns(source.getCurrentRow(), initialSize);
			if(getCurrentRow().attachments == null)getCurrentRow().attachments = source.getCurrentRow().attachments;
			else
				if(source.getCurrentRow().attachments != null)getCurrentRow().attachments.putAll(source.getCurrentRow().attachments);			
		}

		return this;
	}


	/**
	 *  Adding new Table rows at the end of this table
	 *  The newTable may have another number of columns
	 * @param newTable
	 * @return
	 * @throws IOException
	 */

	public DataTable addRowsFrom(DataTable newTable) throws IOException
	{	
		if (!columns.equals(newTable.columns))
			return addRowsWithDifferentColumns(newTable);

		for (newTable.reset();newTable.nextRow();){
			data.add(newTable.getCurrentRow());
			currentRow = data.size() - 1;
		}
		return this;
	}

	private DataTable addRowsWithDifferentColumns(DataTable newTable) throws IOException
	{
		Map<String,Integer> cachedColumns = new HashMap<String,Integer>();

		for(int i=0;i<columns.size();i++)
			cachedColumns.put(columns.get(i),i);

		for (String col:newTable.columns)
			if(!cachedColumns.containsKey(col)){
				cachedColumns.put(col,columns.size());
				columns.add(col);
			}

		for (newTable.reset();newTable.nextRow();){
			AbstractDataRow row = newTable.getCurrentRow();

			if(row.isError()){
				addRow(row);
				continue;
			}

			AbstractDataRow newrow = newTable.compactRowFormat? addRowFastCompact((CompactDataRow)row, newTable.columns, cachedColumns) :
				addRowFast(row, newTable.columns, cachedColumns);

			newrow.attachments = row.attachments;
			newrow.status = row.status;
			newrow.detailedStatus = row.detailedStatus;
		}

		return this;
	}


	public void removeErrors()
	{
		if(getRowsNoErrorsSize() == getRowsSize()) return;

		List<AbstractDataRow> dataNew = new ArrayList<AbstractDataRow>();

		for(int i=0;i<data.size();i++)
			if(!data.get(i).isError())
				dataNew.add(data.get(i));
		data = dataNew;
	}


	/**
	 * Provides fast addition of new rows with potentially different columns 
	 * cachedColumns should contain column name - id for both new and old columns
	 */
	private AbstractDataRow addRowFast(AbstractDataRow newRow, List<String> newColumns, Map<String, Integer> cachedColumns)
	{
		AbstractDataRow row = createRow();
		for (int i = 0; i < newColumns.size(); i++)
			row.setValue(cachedColumns.get(newColumns.get(i)), newRow.getValue(i));
		addRow(row);
		return row;
	}

	private AbstractDataRow addRowFastCompact(CompactDataRow newRow, List<String> newColumns, Map<String, Integer> cachedColumns)
	{
		Map<Integer,Float> vals= new TreeMap<Integer,Float>();

		AbstractDataRow row = createRow();
		float oldValues[] = newRow.toArray();
		for(int i=0;i<oldValues.length;i++)
			if(oldValues[i] != 0)
				vals.put(cachedColumns.get(newColumns.get(i)), oldValues[i]);

		Set<Entry<Integer, Float>> set = vals.entrySet();

		for(Iterator<Entry<Integer, Float>> i = set.iterator(); i.hasNext();) {
			Map.Entry <Integer, Float> me = (Map.Entry<Integer, Float>)i.next();
			row.setValue(me.getKey(), me.getValue());
		}

		addRow(row);
		return row;
	}

	public void reset()
	{
		currentRow = -1;
	}

	public boolean nextRow()
	{
		if (currentRow == data.size() - 1)
			return false;
		else
		{
			currentRow++;
			return true;
		}
	}

	public void forceNextRow()
	{
		if (currentRow >= data.size() - 1)
			throw new UserFriendlyException("There are not sufficient rows in the datatable!");
		currentRow++;
	}

	public Serializable getValue(String colname)
	{
		return getValue(currentRow, colname);
	}

	public Serializable getValue(int column)
	{
		return getValue(currentRow, column);
	}

	@XmlTransient
	public Serializable getValue()
	{
		return getValue(0);
	}

	@XmlTransient
	public String getStatus()
	{
		return data.get(currentRow).status;
	}

	public DataTable clone() throws CloneNotSupportedException
	{
		DataTable dtResult = (DataTable) super.clone();
		dtResult.currentRow = 0;
		return dtResult;
	}

	public DataTable sort(Comparator<AbstractDataRow> c)
	{
		Collections.sort(data, c);
		return this;
	}

	public static DataTable fromArray(Serializable[] values)
	{
		DataTable dt = new DataTable();
		dt.addRow();
		Integer num = 0;
		for (Serializable value : values)
		{
			String column = (++num).toString();
			dt.addColumn(column);
			dt.setValue(column, value);
		}
		return dt;
	}

	public List<Serializable> columnToListUnique(int column)
	{
		List<Serializable> result = new ArrayList<Serializable>();
		reset();
		while (nextRow())
		{
			Serializable value = getValue(column);
			if (!result.contains(value))
				result.add(value);
		}
		return result;
	}

	public void normalize()
	{
		reset();
		while (nextRow())
			getCurrentRow().setWidth(columns.size());
	}

	/**
	 *  Round values to have better compression
	 * @param dtDescriptors
	 * @throws IOException
	 */

	public void roundValues() throws IOException
	{
		if(!isCompactRowFormat())throw new IOException("Not isCompactRowFormat");	
		for(int row = 0; row < data.size() ; row++)
		{
			CompactDataRow drow = (CompactDataRow)data.get(row);
			float vals[] = drow.toArray();
			for (int i = 0; i < vals.length; i++)
				if(vals[i]!=0)
					drow.setValue(i, NumericalValueStandardizer.getSignificantDigits(vals[i]));
		}
	}

	public void fillNullsWith(Serializable filling) throws Exception
	{
		normalize();
		reset();
		while (nextRow())
		{
			for (int i = 0; i < columns.size(); i++)
			{
				Object value = getValue(i);
				if (value == null || value instanceof EmptyCell)
					setValue(i, filling);
			}
		}
	}

	private int statusCount(String status)
	{
		int count = 0;
		reset();
		while (nextRow())
			if (status.equals(getCurrentRow().status))
				count++;
		return count;
	}

	public int errorCount()
	{
		return statusCount(QSPRConstants.ERROR_STATUS);
	}

	public String toString()
	{
		return columns.toString() + "\n" + data.toString();
	}

	public String toStringColumns()
	{
		return columns.toString();
	}

	public void printDebugOneMol(PrintStream stream)
	{
		for(int i=0;i<columns.size();i++)
			stream.println(columns.get(i)+" "+getValue(0, i));
		stream.flush();
	}

	public void print(PrintStream stream)
	{
		reset();
		stream.println(this);
		stream.flush();
	}

	public DataTable setId(String id)
	{
		this.id = id;
		return this;
	}

	public boolean hasMoreRows()
	{
		return currentRow < data.size() - 1;
	}

	public DataTable getCopy()
	{
		return replicateDataTable(0, getRowsSize(), true, true, true);
	}

	/**
	 * Provides a completely new table 
	 * @return
	 */

	public DataTable getDeepCopy() {
		return (DataTable) ClassCompressor.byteToObject(ClassCompressor.objectToByte(this));
	}

	public DataTable getSlice(int fromIndex, int toIndex)
	{
		return replicateDataTable(fromIndex, toIndex, true, true, true);
	}

	// Create an empty copy of this DataTable
	public DataTable getEmptyCopy()
	{
		return replicateDataTable(0, 0, false, false, false);
	}

	/**
	 *  Provides full set of options to replicate DataTable 
	 *  Checks consistency of the options
	 * @param fromIndex
	 * @param toIndex
	 * @param addRows
	 * @param copyRowsStatus
	 * @param copyRowsData
	 * @return
	 */
	public DataTable replicateDataTable(int fromIndex, int toIndex, boolean addRows, boolean copyRowsStatus, boolean copyRowsData)
	{
		DataTable dt = new DataTable(compactRowFormat);
		dt.id = id;
		dt.columns.addAll(columns);
		dt.attachment = attachment;
		dt.columnAttachments = columnAttachments;

		if (copyRowsStatus && !addRows)
			throw new UserFriendlyException("CopyRowsStatus should be set only together with addRows");

		if (copyRowsData && !copyRowsStatus)
			throw new UserFriendlyException("copyRowsData should be set only together with copyRowsStatus");

		if (!addRows)
			return dt;

		if (fromIndex < 0)
			throw new UserFriendlyException("fromIndex = " + fromIndex + " < 0");
		if (fromIndex > toIndex)
			throw new UserFriendlyException("fromIndex = " + fromIndex + " > toIndex = " + toIndex);
		if (toIndex > getRowsSize())
			throw new UserFriendlyException("toIndex = " + toIndex + " > getRowsSize() = " + getRowsSize());

		for (int i = fromIndex; i < toIndex; i++)
		{
			if (copyRowsData)
				dt.addRow(getRow(i));
			else
			{
				dt.addRow();
				if (copyRowsStatus)
				{
					dt.getRow(i - fromIndex).status = getRow(i).status;
					dt.getRow(i - fromIndex).detailedStatus = getRow(i).detailedStatus;
				}
			}
		}
		return dt;
	}

	public DataTable setColumnAttachment(String column, String id, Serializable object)
	{
		Map<String, Serializable> map = columnAttachments.get(column);
		if (map == null)
		{
			map = new HashMap<String, Serializable>();
			columnAttachments.put(column, map);
		}
		map.put(id, object);
		return this;
	}

	public Serializable getColumnAttachment(String column, String id)
	{
		if (columnAttachments.get(column) == null)
			return null;
		return columnAttachments.get(column).get(id);
	}

	/** Returns submatrix from row "from" (inclusive) to row "to" (exclusive). */
	public final float[][] asFloatArray(final int from, final int to) throws Exception
	{
		final int n = getRowsSize();
		if (n == 0)
			throw new Exception("Expected DataTable to be a matrix, but it has no rows.");
		if (from < 0 || to > n || from > to)
			throw new Exception("Invalid number of rows specified.");
		final int m = columns.size();
		if (m == 0)
			throw new Exception("Expected DataTable to be a matrix, but it has no columns.");

		final float[][] res = new float[to - from][m];
		for (int i = from; i < to; ++i)
			for (int j = 0; j < m; ++j)
				res[i - from][j] = ((Double) getValue(i, j)).floatValue();

		return res;
	}

	public float[][] asFloatArray(final int numrows) throws Exception
	{
		return asFloatArray(0, numrows);
	}

	public float[][] asFloatArray() throws Exception
	{
		return asFloatArray(0, getRowsSize());
	}

	public String getSDF() throws IOException{
		Serializable s = getValue(currentRow, 0);
		if(!(s instanceof String)) throw new IOException("not String");
		return (String)s;
	}

	public String getSDF(int row) throws IOException{
		Serializable s = getValue(row, 0);
		if(!(s instanceof String)) throw new IOException("not String");
		return (String)s;
	}

	/** Returns structures in column and a given range of rows of the data table as one string. */
	private String asSdfString(final int from, final int to) throws Exception
	{
		// Check parameters.
		final int n = getRowsSize();
		if (n == 0)
			throw new Exception("Expected DataTable to contain structures, but it has no rows.");
		if (from < 0 || to > n || from > to)
			throw new Exception("Invalid number of rows specified.");

		// This DataTable can come from different sources, all of which may give a slightly different table layot, column names, etc.
		int colInd = 0; // Default is first column
		if (columns.size() != 1)
		{
			// Look for a column named SDF.
			colInd = -1;
			for (int i = 0; i < columns.size(); ++i)
			{
				final String colName = columns.get(i);
				if (colName.equals(QSPRConstants.SDF_COLUMN))
				{
					colInd = i;
					break;
				}
			}

			// Look for a column containing strings.
			colInd = -1;
			for (int i = 0; i < columns.size(); ++i)
			{
				if (getValue(i, 0).getClass().equals(String.class))
				{
					colInd = i;
					break;
				}
			}

			if (colInd == -1)
				throw new Exception("More than one column and could not identify SDF column.");
		}
		if (!(getValue(0, 0) instanceof String))
			throw new Exception("Column does not contain strings.");

		// Build string.
		final StringBuilder sb = new StringBuilder();
		for (int i = from; i < to; ++i)
		{
			sb.append((String) getValue(i, colInd));
			sb.append("\n\n$$$$\n");
		}

		return sb.toString();
	}

	/** Returns structures starting from row "from" (inclusive) until row "to" (exclusive) as an array of Strings. */
	public final String[] asStructures(final int from, final int to) throws Exception
	{
		if (from > to)
			throw new RuntimeException(format("Invalid range specification (%d,%d)", from, to));
		if (from == to)
			return new String[] {};
		final String[] result = asSdfString(from, to).split("\\$\\$\\$\\$\n");
		for (int i = 0; i < to - from; ++i)
			result[i] = result[i].concat("$$$$\n");
		result[result.length - 1] = result[result.length - 1].substring(0, result[result.length - 1].length() - 1);
		return result;
	}

	/** Returns structures of the data table as an array of Strings. */
	public final String[] asStructures(final int numrows) throws Exception
	{
		return asStructures(0, numrows);
	}

	/** Returns structures of the data table as an array of Strings. */
	public final String[] asStructures() throws Exception
	{
		return asStructures(0, getRowsSize());
	}

	public static DataTable getDummy(int numRows)
	{
		DataTable dt = new DataTable();
		dt.id = "dummyTableToFillPort";
		for (int i = 0; i < numRows; i++)
			dt.addRow();
		return dt;
	}

	// TODO remove as soon as possible!

	public int getClassificationID(int molecule)
	{
		return (int) Math.rint((Double) getValue(molecule, "CLASS"));
	}

	/**
	 * Delete columns which are not required
	 * @param deleteColumns
	 */

	public void deleteByList(Set<String> deleteColumns)
	{
		List<String> keepColumns = new ArrayList<String>(columns.size());
		for (Integer i = 0; i < columns.size(); i++)
			if (!deleteColumns.contains(columns.get(i)))
				keepColumns.add(columns.get(i));

		keepByList(keepColumns);
	}

	public void addColumns(List <String> columns){
		for(String col:columns)
			addColumn(col);
	}

	/**
	 * Rearranges columns according to the provided list of descriptors
	 * In case if some descriptors are not available, add 0 values for all of them
	 */

	public void keepByList(List<String> descriptors)
	{
		if (descriptors == null)
			return;

		// it is possible that nothing needs to be changed
		boolean theSame = descriptors.size() == columns.size();
		for (int i = 0; i < descriptors.size() && theSame; i++)
			theSame = descriptors.get(i).equals(columns.get(i));

		if (theSame)
			return;

		List<String> oldColumns = columns;
		List<AbstractDataRow> oldData = data;

		columns = new ArrayList<String>(descriptors);
		data = new ArrayList<AbstractDataRow>();

		int newToOld[] = new int[columns.size()];
		for (int i = 0; i < columns.size(); i++)
			newToOld[i]=oldColumns.indexOf(columns.get(i));

		for (int r = 0; r < oldData.size(); r++)
		{
			CompactDataRow oldRow = (CompactDataRow)oldData.get(r);
			float oldValues[]=oldRow.toArray();

			AbstractDataRow newRow = addRow();
			newRow.attachments = oldRow.attachments;

			if(oldRow.isError())newRow.setError(oldRow.detailedStatus);
			else
				for (int i = 0; i < columns.size(); i++)
				{
					int index = newToOld[i];
					if(index >= oldValues.length || index == -1)continue;
					newRow.setValue(i, oldValues[index]);
				}

			oldData.set(r, null); // in order do not keep two tables in memory simultaneously
		}
	}

	@Override
	public Iterator<AbstractDataRow> iterator() 
	{
		return new DataTableIterator();
	}

	class DataTableIterator implements Iterator<AbstractDataRow>
	{
		public DataTableIterator()
		{
			reset();
		}

		@Override
		public boolean hasNext() 
		{
			return hasMoreRows();
		}

		@Override
		public AbstractDataRow next() 
		{
			nextRow();
			return getCurrentRow();
		}

		@Override
		public void remove() {
			throw new RuntimeException("Unimplemented");

		}

	}

	/**
	 * This can be rather long operation
	 * Avoid using this function often
	 * @return
	 */

	public DataSize getValuesCount() {
		DataSize size = new DataSize();
		for(int row=0;row<getRowsSize();row++)
			for(int col=0;col<getColumnsSize();col++){
				size.all++;
				Serializable val = getValue(row, col);
				if(val != null && val instanceof Double && (Double)val != 0.)
					size.nonZero++;
			}

		return size;
	}

	/**
	 *  Substitute row with a new one
	 * @param index
	 * @param row
	 */

	public void setRow(int index, AbstractDataRow row) {
		data.set(index, row);
	}

	/**
	 *  Substitute row with a new one
	 * @param index
	 * @param row
	 */

	public void setColumnName(int column, String name) {
		columns.set(column, name);	
	}

	public boolean containsColumn(String column) {
		return columns.contains(column);
	}

	public String getColumn(int column) {
		return columns.get(column);
	}

	public int indexOfColumn(String name) {
		return columns.indexOf(name);
	}

	final public List<String> getColumns() {
		return columns;
	}

	public void allowNonStored() {
		for (reset();nextRow();)
			if(getCurrentRow().isError())
			{
				if(getCurrentRow().detailedStatus.contains(QSPRConstants.NOSTORED) || getCurrentRow().detailedStatus.contains(QSPRConstants.EMPTY_MOLECULE_ERROR)){
					getCurrentRow().setError(null);
					getCurrentRow().setStatus(null);
				}else
					System.out.println("ignoring row with error: " + getCurrentRow().detailedStatus);
			}
	}

	public void markNoDescriptorsAsErrors() {
		for (reset();nextRow();){
			AbstractDataRow row = getCurrentRow();
			if(row.size() == 0 && !getCurrentRow().isError())
				getCurrentRow().setError(QSPRConstants.NOSTORED + " or were calculated");
		}
	}

	public String compareColumns(DataTable compareColumns) {
		if(compareColumns.columns.size() != columns.size()) {
			String s = "";
			for(String column : compareColumns.columns) if(!columns.contains(column))
				s += "new: " + column+ ";";
			for(String column : columns) if(!compareColumns.columns.contains(column))
				s += "old: " + column+ ";";

			return " different columns sizes " + s;
		}

		for(int i =0; i< columns.size();i++)
			if(!columns.get(i).equals(compareColumns.columns.get(i)))
				return columns.get(i) + " != " + compareColumns.columns.get(i);
		return null;
	}

	public void filterNaN(int row) {
		int good = 0;
		CompactDataRow drow = (CompactDataRow)data.get(row);
		if(drow.isError())return;
		float vals[] = drow.toArray();
		if(vals.length == 0)return; // good, but all values are 0
		for (int i = 0; i < vals.length; i++)
			if(Double.isFinite(vals[i])) good++;
		if(good == 0) {
			drow.setError("All descriptors were NaN for row: " + row);
			return;
		}

		vals = drow.toArray();
		String erorr="";
		for (int i = 0; i < vals.length; i++)
			if(!Double.isFinite(vals[i])) 
				erorr += columns.get(i)+",";
		if(erorr.length()>1)
			drow.setError("The value of descriptor(s) \"" + erorr + "\" is not a number (NaN)");
	}

	public void filterNaN() {
		for(int row = 0; row < data.size() ; row++)
			filterNaN(row);
	}

	public DataTable deleteFailed(boolean[] skip) {
		DataTable table = getEmptyCopy();
		for(int i=0;i < getRowsSize()  && (skip == null || i<skip.length) ;i++)
			if(skip == null || !skip[i])
				table.addRow(getRow(i));
		return table;
	}

	public int getMinValueRowId(int column, int startFrom) {
		if(getColumnsSize()<column)return 0;
		double val = 1e10; int iter=0;
		for(int i = startFrom;i<getRowsSize()-1;i++) {
			Double newval = (Double)getValue(i, column);
			if(newval == null)continue;
			if(newval<val*0.99){
				val = newval;
				iter = i;
			}
		}		
		return iter;
	}

	public DataTable decompress(int[] compressing) {
		DataTable res =  getEmptyCopy();
		for(int i=0;i<compressing.length;i++)
			res.addRow((AbstractDataRow)ClassCompressor.byteToObject(ClassCompressor.objectToByte(data.get(compressing[i]))));
		return res;
	}

	public int getMaxValueRowId(int column, int startFrom) {
		if(getColumnsSize()<column)return 0;
		double val = -1e10; int iter=0;
		for(int i = startFrom;i<getRowsSize()-1;i++) {
			Double newval = (Double)getValue(i, column);
			if(newval == null)continue;
			if(val>newval)continue;
			val = newval;
			iter = i;
		}		
		return iter;
	}

	public void keepResults() {
		List<String> keep = new ArrayList<String>();
		for(String s: getColumns())
			if(s.startsWith(QSPRConstants.PREDICTION_RESULT_COLUMN))
				keep.add(s);
		keepByList(keep);
	}

	/**
	 * Prediction contains non supported descriptors, e.g. ADRIANA
	 * They are ignored to avoid model supports
	 * @return
	 */

	public boolean containsMissedDescriptors() {
		reset();
		while (nextRow())
			if(getCurrentRow().isError() && QSPRConstants .MISSED_DESCRIPTOR.equals(getCurrentRow().detailedStatus)) return true;

		return false;
	}

	public ChemInfEngine getChemInfEngine() { // should have different names (set/get) are interpreted by Java bean classes incorrectly
		String engine =  columnAttachments.get(QSPRConstants.SDF_COLUMN) == null? null: 
			(String) columnAttachments.get(QSPRConstants.SDF_COLUMN).get(QSPRConstants.CHEMENGINE);
		engine = engine == null ? "" + ChemInfEngine.CHEMAXON : engine;
		return ChemInfEngine.valueOf(engine);
	}

	public void setInfEngine(ChemInfEngine engine) throws Exception{
		if(columns == null || !columns.contains(QSPRConstants.SDF_COLUMN))
			throw new Exception("datatable does not contain " + QSPRConstants.SDF_COLUMN);
		Map<String, Serializable> map = new HashMap<String,Serializable>();
		map.put(QSPRConstants.CHEMENGINE, "" + Various.molecule.engine);
		columnAttachments.put(QSPRConstants.SDF_COLUMN, map);	
	}

}
