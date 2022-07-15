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

package com.eadmet.useractions;

import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.annotation.XmlRootElement;

import com.eadmet.batchupload.entityschema.EntitiesRemapping;
import com.eadmet.batchupload.entityschema.PropertyRemapping;
import com.eadmet.batchupload.main.UploadPreview;

/**
 * A documentation data for a bacth upload action
 * @author midnighter, novserj
 *
 */
@XmlRootElement
public abstract class BatchUpload30Action extends AbstractUserAction
{
	public BatchUpload30Action()
	{
		
	}
	
	public static BatchUpload30Action previewStart(EntitiesRemapping schema)
	{
		BatchUpload30StartPreviewAction a = new BatchUpload30StartPreviewAction(schema);
		return a;
	}
	
	public static BatchUpload30Action uploadStart(EntitiesRemapping schema, UploadPreview preview)
	{
		BatchUpload30StartUploadAction a = new BatchUpload30StartUploadAction(schema, preview);
		return a;
	}
	
	public static BatchUpload30Action uploadFinish(EntitiesRemapping schema, UploadPreview preview)
	{
		BatchUpload30FinishUploadAction a = new BatchUpload30FinishUploadAction(schema, preview);
		return a;
	}
}

@XmlRootElement
class BatchUpload30StartPreviewAction extends BatchUpload30Action
{
	public List<String> properties;
	public long count;
	
	public boolean preview;
	
	@Override
	public String getLogLine()
	{
		return String.format("is starting upload preview with %d records for properties %s", count, properties.toString());
	}
	
	public BatchUpload30StartPreviewAction()
	{
		
	}

	public BatchUpload30StartPreviewAction(EntitiesRemapping schema)
	{
		properties = new ArrayList<String>();
		for (PropertyRemapping property : schema.properties)
			properties.add(property.name);
		this.count = schema.records;
	}
}

@XmlRootElement
class BatchUpload30StartUploadAction extends BatchUpload30Action
{
	public List<String> properties;
	public long count;
	public String summary;
	
	@Override
	public String getLogLine()
	{
		return String.format("is starting upload with %d records for properties %s, summary = %s", count, properties.toString(), summary);
	}
	
	public BatchUpload30StartUploadAction()
	{
		
	}

	public BatchUpload30StartUploadAction(EntitiesRemapping schema, UploadPreview preview)
	{
		properties = new ArrayList<String>();
		for (PropertyRemapping property : schema.properties)
			properties.add(property.name);
		this.count = schema.records;
		this.summary = preview.getSummary().toString();
	}
}

@XmlRootElement
class BatchUpload30FinishUploadAction extends BatchUpload30Action
{
	public List<String> properties;
	public long count;
	public String summary;
	
	@Override
	public String getLogLine()
	{
		return String.format("has finished upload with %d records for properties %s, summary = %s", count, properties.toString(), summary);
	}
	
	public BatchUpload30FinishUploadAction()
	{
		
	}

	public BatchUpload30FinishUploadAction(EntitiesRemapping schema, UploadPreview preview)
	{
		properties = new ArrayList<String>();
		for (PropertyRemapping property : schema.properties)
			properties.add(property.name);
		this.count = schema.records;
		this.summary = preview.getSummary().toString();
	}
}
	
