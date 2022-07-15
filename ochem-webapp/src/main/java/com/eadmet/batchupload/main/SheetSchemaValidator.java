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

public class SheetSchemaValidator 
{
	UploadedSheetSchema sheet;
	public void validate(UploadedSheetSchema sheet)
	{
		this.sheet = sheet;
		sheet.messages.clear();
		List<String> fakeValues = new ArrayList<String>();
		for (UploadedColumnSchema column : sheet.getColumns())
			fakeValues.add("stub_"+column.name);
		UploadedRow fakeRow = new UploadedRow(1, fakeValues, sheet);
		RecordStubBuilder rsb = new RecordStubBuilder();
		List<RecordStub> stubs = rsb.getRowStubs(fakeRow);
		
		if (stubs.size() > 1)
		{
			List<String> propertyName = new ArrayList<String>();
			for (RecordStub stub : stubs)
				if (stub.property.name != null)
					propertyName.add(stub.property.name.value);
			sheet.messages.add(BatchUploadMessage.newNotice("Several records per row will be created for this sheet: "+propertyName));
		}
		
		for (RecordStub stub : stubs)
			validateStub(stub);
	}
	
	private void validateStub(RecordStub stub)
	{
		sheet.messages.clear();
		sheet.missingProperty = false;
		if (stub.property.name == null)
		{
			sheet.messages.add(BatchUploadMessage.newError("The sheet does not have a recognized property"));
			sheet.missingProperty = true;
		}
		else
		if (stub.property.value.size() == 0 || stub.property.value.getFirst().value == null)
			sheet.messages.add(BatchUploadMessage.newError("The sheet does not have a recognized value column for property, column #"+(stub.property.name.column.realNumber + 1)));

		
		if (stub.molecules.size() == 0 && stub.name.size() == 0 && stub.externalid.size() == 0)
			sheet.messages.add(BatchUploadMessage.newError("The MOLECULE or NAME or EXTERNALID columns are missing"));
		
		if (stub.molecules.size() + stub.name.size() > 1)
			sheet.messages.add(BatchUploadMessage.newNotice("Several MOLECULE or NAME columns are present."));
		
		if (stub.article.size() == 0 || stub.article.getFirst().column == null)
			sheet.messages.add(BatchUploadMessage.newNotice("The ARTICLE column is missing, the stub unpublished article will be assigned by default"));
	}
}
