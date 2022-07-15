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

import java.io.IOException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.Globals;
import qspr.util.WrapperThread;

import com.eadmet.batchupload.entityschema.EntitiesRemapping;

public class MultithreadBatchUploadProcessor extends BatchUploadProcessor
{
	WrapperThread currentThread;
	IOException exc;
	private static final Logger logger = LogManager.getLogger(MultithreadBatchUploadProcessor.class);
	public void waitToFinish() throws InterruptedException
	{
		Thread.sleep(100);
		if (isThreadRunning())
			currentThread.join();
	}
	
	public boolean isThreadRunning()
	{
		return (currentThread != null && currentThread.isRunning());
	}
	
	@Override
	boolean isRemappingReady()
	{
		return (super.isRemappingReady() && !isThreadRunning());
	}
	
	@Override
	boolean isPreviewReady()
	{
		return (super.isPreviewReady() && !isThreadRunning());
	}
	
	@Override
	boolean isUploadReady()
	{
		return (super.isUploadReady() && !isThreadRunning());
	}
	
	@Override
	public EntitiesRemapping getRemappingSchema() throws IOException
	{
		if (!isRemappingReady())
		{
			if (!isThreadRunning())
			{
				if (exc != null)
					throw exc;
				threadedRemappingSchema();
			}
			return null;
		}
		return super.getRemappingSchema();
	}
	
	@Override
	public UploadPreview getUploadPreview() throws IOException
	{
		if (!isPreviewReady())
		{
			if (!isThreadRunning())
			{
				if (exc != null)
					throw exc;
				threadedUploadPreview();
			}
			return null;
		}
		return super.getUploadPreview();
	}
	
	@Override
	public UploadPreview upload() throws IOException
	{
		if (!isUploadReady())
		{	
			if (!isThreadRunning())
			{
				if (exc != null)
					throw exc;
				threadedUpload();
			}
			return null;
		}
		return super.upload();
	}
	
	
	void threadedRemappingSchema() throws IOException
	{
		if (currentThread != null && currentThread.isRunning())
			return;
		currentThread = new WrapperThread()
		{
			public void wrapped() throws Exception
			{
				logger.info("Starting processRemappingSchema");
				try
				{
					processRemappingSchema();
				} catch (Exception e)
				{
					exc = new IOException(e);
					e.printStackTrace();
				}
			}			
		};
		currentThread.setName("processRemappingSchema");
		currentThread.userSession = Globals.userSession();
//		currentThread.exclusiveLock = true;
		currentThread.updateSessionTime = false;
		currentThread.start();
	}
	
	void threadedUploadPreview() throws IOException
	{
		if (currentThread != null && currentThread.isRunning())
			return;
		
		currentThread = new WrapperThread()
		{
			public void wrapped() throws Exception
			{
				try
				{
					logger.info("Starting processUploadPreview");
					processUploadPreview();
				} catch (Exception e)
				{
					exc = new IOException(e);
					e.printStackTrace();
				}
			}			
		};
		currentThread.setName("processUploadPreview");
		currentThread.userSession = Globals.userSession();
//		currentThread.exclusiveLock = true;
		currentThread.updateSessionTime = false;
		currentThread.start();
	}
	
	void threadedUpload() throws IOException
	{
		if (currentThread != null && currentThread.isRunning())
			return;
		
		currentThread = new WrapperThread()
		{
			public void wrapped() throws Exception
			{
				try
				{
					logger.info("Starting processUpload");
					processUpload();
				} catch (Exception e)
				{
					exc = new IOException(e);
					e.printStackTrace();
				}
			}			
		};
		currentThread.setName("processUpload");
		currentThread.userSession = Globals.userSession();
//		currentThread.exclusiveLock = true;
		currentThread.updateSessionTime = false;
		currentThread.start();
	}
}