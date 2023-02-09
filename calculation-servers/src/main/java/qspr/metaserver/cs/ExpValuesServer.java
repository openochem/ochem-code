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

package qspr.metaserver.cs;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.List;

import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.configurations.DescriptorsExpValuesConfiguration;
import qspr.metaserver.protocol.Task;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;
import qspr.workflow.utils.QSPRConstants;

public class ExpValuesServer extends DescriptorsAbstractServer
{
	public String dbUrl = "jdbc:mariadb://iprior/ochem_clean";
	public String dbUsername = "root";
	public String dbPassword = "TnOtB1982";

	public ExpValuesServer()
	{
		supportedTaskType = "ExpValues";
		setInputFlowGroup(0);
		setOutputFlowGroup(0);

		try
		{
			Class.forName(QSPRConstants.DEFAULTDATABASEDRIVER).newInstance();
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	}

	public void setParam(String name, String value)
	{
		super.setParam(name, value);
		name = name.toUpperCase();
		if ("DBURL".equals(name))
			dbUrl = value;
		else if ("DBUSERNAME".equals(name))
			dbUsername = value;
		else if ("DBPASSWORD".equals(name))
			dbPassword = value;
	}

	Connection conn;
	Long connectionOpenedTime;
	private Connection getConnection() throws SQLException, InterruptedException
	{
		for (int i = 0; i < 10; i++)
			try
		{
				connectionOpenedTime = Calendar.getInstance().getTimeInMillis();
				if (conn == null || conn.isClosed())
				{
					out.println("Connecting to " + dbUrl);
					conn = DriverManager.getConnection(dbUrl, dbUsername, dbPassword);
				}
				conn.prepareStatement("select mapping2_id from Mapping2 limit 1").execute(); //Test query to make sure the connection is not broken
				return conn;
		}
		catch (Exception e)
		{
			conn = null;
			if (i == 9)
				throw new RuntimeException("Could not create a connection after 10 attempts", e);
			Thread.sleep(1000);
		}
		throw new RuntimeException("Could not create a connection after 10 attempts");
	}

	@Override
	public WorkflowNodeData calculateDescriptors(WorkflowNodeData task,
			DescriptorsAbstractConfiguration configuration) throws Exception
	{
		// Fetch exp. values by molecule IDs
		DescriptorsExpValuesConfiguration conf = (DescriptorsExpValuesConfiguration) configuration;
		DataTable dtMolecules = task.ports.get(0);
		String propertyNames = "";
		for (String propName : conf.properties) {
			propertyNames += "\"" + propName + "\", ";
		}
		propertyNames = propertyNames.substring(0, propertyNames.length() - 2);

		setStatus("Quering the database");
		out.println("Property names: " + propertyNames);
		out.println("Basket ID: " + conf.basketId);
		PreparedStatement statement = getConnection().prepareStatement("select mapping2_id, Property.name, Property.type, canonical_value, PropertyOption.name from ExperimentalProperty natural left join Molecule inner join Property using (property_id) left join PropertyOption using (poption_id) inner join BasketEntry using (exp_property_id) where Property.property_id in (select property_id from Property pr_inner where pr_inner.name in ("+propertyNames+")) and basket_id=?");
		statement.setLong(1, conf.basketId);
		ResultSet rs = statement.executeQuery();

		List<Integer> molIds = new ArrayList<Integer>();
		int[] molRefs = new int[dtMolecules.getRowsSize()];

		dtMolecules.reset();
		while (dtMolecules.nextRow())
		{
			Integer mp2 = (Integer)dtMolecules.getCurrentRow().getAttachment(QSPRConstants.MOLECULE_ID_STEREOCHEM);
			int ref = molIds.indexOf(mp2);
			if (ref != -1)
				molRefs[dtMolecules.currentRow] = ref;
			else
			{
				molIds.add(mp2);
				molRefs[dtMolecules.currentRow] = molIds.size() - 1;
			}
		}

		float[][] values = new float[molIds.size()][conf.properties.size()]; // [molecule][property] = property value
		for (float[] fs : values) {
			Arrays.fill(fs, Float.MIN_VALUE);
		}

		setStatus("Processing the query results");
		while (rs.next())
		{
			Integer molNum = molIds.indexOf(rs.getInt(1));
			Integer propNum = conf.properties.indexOf(rs.getString(2));
			if (molNum != -1 && propNum != -1)
			{
				float existingValue = values[molNum][propNum];
				int type = rs.getInt(3);
				float dbValue;
				if (type == 0)
					dbValue = rs.getFloat(4);
				else
					dbValue = numericValueOf(rs.getString(5));

				if (existingValue == Float.MAX_VALUE) // This value is an error. Ignore new values
					continue;

				if (existingValue != Float.MIN_VALUE) // We already have a value, consider several scenarios
				{
					if ("FAIL".equals(conf.multipleExpValues))
					{
						values[molNum][propNum] = Float.MAX_VALUE;
					}
					else if ("USE_MAXIMUM".equals(conf.multipleExpValues))
					{
						if (existingValue <  dbValue)
							values[molNum][propNum] = dbValue;
					}
					else if ("USE_MINIMUM".equals(conf.multipleExpValues))
					{
						if (existingValue > dbValue)
							values[molNum][propNum] = dbValue;
					}
				}
				else
					values[molNum][propNum] = dbValue; // We have a value
			}
		}


		DataTable dtResults = new DataTable(true);
		for (String propertyName : conf.properties)
			dtResults.addColumn(propertyName);

		for (int i = 0; i < dtMolecules.getRowsSize(); i++) {
			dtResults.addRow();
			for (int pNum = 0; pNum < conf.properties.size(); pNum++) {
				if (values[molRefs[i]][pNum] == Float.MIN_VALUE)
					dtResults.getCurrentRow().setError("No value for property " + conf.properties.get(pNum));
				else if (values[molRefs[i]][pNum] == Float.MAX_VALUE)
					dtResults.getCurrentRow().setError("Dublicate value for property " + conf.properties.get(pNum));
				else
					dtResults.setValue(pNum, 1D * values[molRefs[i]][pNum]);
			}
		}

		return new WorkflowNodeData(dtResults);
	}

	// A "patch" function to retrieve numeric value based on qualitative option (so far, "Yes" and "No")
	// This is not a universal solution, but is enough for Ahmed's needs
	float numericValueOf(String option)
	{
		if ("Yes".equalsIgnoreCase(option))
			return 1;
		else if ("No".equalsIgnoreCase(option))
			return 0;
		else
		{
			out.println("WARNING: Unknown property option " + option);
			return Float.MAX_VALUE; // error
		}

	}

	public static void main(String[] args) throws Exception
	{
		DescriptorsExpValuesConfiguration evConf = new DescriptorsExpValuesConfiguration();
		evConf.properties.add("Solidus (P450):Solidus_P450");
		evConf.basketId = 1000000678L;

		DataTable dtMols = new DataTable();
		dtMols.addColumn("MOL");
		dtMols.addRow();
		dtMols.getCurrentRow().addAttachment(QSPRConstants.MOLECULE_ID_STEREOCHEM, 1000000044L);

		Task task = new Task("ExpValues", evConf, new WorkflowNodeData(dtMols));
		new ExpValuesServer().calculate(task);
	}
}

