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

package qspr.modelling.configurators;

import java.io.BufferedInputStream;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import javax.servlet.http.HttpServletRequest;

import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;

import qspr.entities.Model;
import qspr.frontend.MarshalableList;
import qspr.frontend.WebModel;
import qspr.metaserver.configurations.CompressedObject;
import qspr.metaserver.configurations.MLRAConfiguration;
import qspr.modelling.configurations.CDSConfiguration;

public class MLRAConfigurator extends DescriptorsConfigurator 
{

	@Override
	public void setModel(Model model)
	{
		super.setModel(model);
		((CDSConfiguration)model.attachment.getObject().configuration).modelConfiguration = new MLRAConfiguration();
		this.stepAfterDescriptors = "mlra";
	}

	public WebModel mlra()
	{
		return new WebModel(getDefaultTemplate()).setTemplate("modeller/configurators/mlra");
	}

	public void mlraSubmit(HttpServletRequest request)
	{
		MLRAConfiguration conf = (MLRAConfiguration) 
				((CDSConfiguration)model.attachment.getObject().configuration).modelConfiguration;
		conf.alpha = Float.parseFloat(request.getParameter("alpha"));
		conf.nvariables = Integer.parseInt(request.getParameter("nvariables"));
		currentPage = startStep;

		if (request.getParameter("leverage") != null)
			conf.addDM("Leverage");

		if (request.getParameter("limitrange") != null)
			conf.limitRange = true;

	}

	@SuppressWarnings("unchecked")
	@Override
	public void processModelUpload(File f) throws IOException 
	{
		MarshalableList preview = new MarshalableList();

		descriptorOverrides.clear();

		Map<String, Double> modelMap = new HashMap<String, Double>();

		double bias = 0.0; // default value

		byte[] xlsData = new byte[(int) f.length()];

		BufferedInputStream bis = new BufferedInputStream(new FileInputStream(f));
		bis.read(xlsData);
		bis.close();

		Workbook wb = new HSSFWorkbook(new ByteArrayInputStream(xlsData));  
		Sheet sheet = wb.getSheet("Linear regression data");

		if (sheet == null)
			sheet = wb.getSheetAt(0);

		Row header = sheet.getRow(0);

		if (!header.getCell(0).toString().equalsIgnoreCase("Descriptor"))
			throw new IOException("Incorrect format: in Excel first cell in the first row should be Descriptor");

		if (!header.getCell(1).toString().equalsIgnoreCase("Coefficient"))
			throw new IOException("Incorrect format: in Excel second cell in the first row should be Coefficient");

		Row biasRow = sheet.getRow(1);

		if (!biasRow.getCell(0).toString().equalsIgnoreCase("Bias"))
			throw new IOException("Incorrect format: in Excel second row should be Bias value");

		int i = 0;
		for (Row row : sheet)
		{
			i++;
			if (i == 1)
				continue;

			String descriptorName = row.getCell(0).toString();

			if (descriptorName.length() == 0)
				break; // there are no more coefficients -- stop!

			String stringValue = row.getCell(1).toString();
			Double value = stringValue.equals("") ? 0 : Double.valueOf(stringValue);

			if (value == 0)
				continue; // the values with 0 will not be added

			modelMap.put(descriptorName, value);

			if (!descriptorName.equalsIgnoreCase("Bias")) // Valid and not the bias
			{
				descriptorOverrides.add(descriptorName);

				MarshalableList mlrow = new MarshalableList();
				mlrow.list.add(descriptorName);
				mlrow.list.add(""+value);
				preview.list.add(mlrow);
			} else
				bias = value;
		}	

		// Adding Bias
		MarshalableList row = new MarshalableList();
		row.list.add("Bias");
		row.list.add(""+bias);
		preview.list.add(0, row); // insert as first row
		MarshalableList prev = new MarshalableList();
		prev.list.add(preview);
		model.preview = prev;
		((MLRAConfiguration)((CDSConfiguration)model.attachment.getObject().configuration).modelConfiguration).uploadedModelData = new CompressedObject<Object>(modelMap);
	}
}
