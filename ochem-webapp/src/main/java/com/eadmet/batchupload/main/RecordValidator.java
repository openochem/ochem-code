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

import qspr.Globals;
import qspr.entities.ExperimentalProperty;
import qspr.util.UploadContext;

import com.eadmet.batchupload.main.RecordPreview.PreviewRecordStatus;
import com.eadmet.batchupload.main.RecordPreview.PreviewUploadAction;

public class RecordValidator 
{
	RecordPreview preview;
	UploadContext context;
	ExperimentalProperty ep;
	
	public void validate(RecordPreview preview, UploadContext context)
	{
		this.preview = preview;
		this.context = context;
		this.ep = preview.ep;
		
		try
		{	
			if (ep == null)
			{
				preview.messages.add(BatchUploadMessage.newError("No experimental property found in preview"));
				return;
			}
		
			ep.initializeBasicCollections();
			
			rangesAndSizesValidation();
			
			basicHashCalculationValidation();
	
			internalDuplicateValidation();

			externalDuplicateValidation();
			
			someValueCompatibilityValidation();

	
			if (ep.md5 != null)
				context.recordHashCache.put(ep.md5, preview.getCoordinates());
			
		} catch (Exception e)
		{
			preview.messages.add(BatchUploadMessage.newError(e));
		} finally
		{
			if (preview.numErrorMessages() > 0)
			{
				if (preview.status == null)
				{
					if (ep != null && ep.molecule != null && ep.property != null && ep.article != null && (ep.property.isNumeric() || ep.option != null))
						preview.status = PreviewRecordStatus.error;
					else
					{
						preview.ep = null;
						preview.status = PreviewRecordStatus.fatal_error;
					}
				}
				if (preview.action == null)
					preview.action = PreviewUploadAction.skip;
			} else
			{
				if (preview.status == null)
					preview.status = PreviewRecordStatus.valid;
				
				if (preview.action == null)
					preview.action = PreviewUploadAction.save;
			}
		}
	}
	
	private void someValueCompatibilityValidation() throws Exception
	{
		if (ep.rights == Globals.RIGHTS_NONE && ep.owner == null)
		{
			preview.messages.add(BatchUploadMessage.newError("Anonymous users are not allowed to create hidden records"));
			ep.ep_status = ExperimentalProperty.STATUS_ERROR;
		}
	}
	
	private void rangesAndSizesValidation() throws Exception
	{
		if (ep.artMolId != null && ep.artMolId.length() > 50)
			preview.messages.add(BatchUploadMessage.newError("Label of the molecule in the article (\"N\" column) has a value longer than 50 characters. Please provide a shorter label."));
	}
	
	private void basicHashCalculationValidation()
	{
		try
		{
			ep.updateHash();
		} catch (Exception e)
		{
			ep.ep_status = ExperimentalProperty.STATUS_ERROR;
			ep.errorComment = e.getMessage();
			ep.updateHash();			
			preview.messages.add(BatchUploadMessage.newError(e));			
		}
	}
	
	private void internalDuplicateValidation() throws Exception
	{
		// "Internal duplicates" matter only in preview mode. On real upload they turn to normal duplicates.
		if (context.preview && ep.md5 != null && context.recordHashCache.get(ep.md5) != null && !context.recordHashCache.get(ep.md5).equals(preview.getCoordinates()))
		{
			String[] coordinates = context.recordHashCache.get(ep.md5).split("_");
			
			preview.messages.add(BatchUploadMessage.newError("Internal duplicate with sheet "+(Integer.valueOf(coordinates[0])+1)+", row "+(Integer.valueOf(coordinates[1])+1)+", record "+(Integer.valueOf(coordinates[2])+1)));
			ep.ep_status = ExperimentalProperty.STATUS_ERROR;
			ep.updateHash();
			preview.status = PreviewRecordStatus.duplicate_internal;
		}
	}
	
	private void externalDuplicateValidation()
	{
		if (ep.hasConflicts() && !ep.duplicate.id.equals(ep.id))
		{
			String error = "";
			
			if (ep.duplicate.rights.equals(Integer.valueOf(2)))
				error = "Duplicate of a public record (RecordID: R"+ep.duplicate.id+")";
			else
				if (ep.duplicate.owner != null)
				{
					if (ep.duplicate.owner.equals(ep.owner))
						error = "Duplicate of your private record (RecordID: R"+ep.duplicate.id+")";
					else
						error = "Duplicate of other user's private record (RecordID: R"+ep.duplicate.id+", User: "+ep.duplicate.owner.login+")";
				} else
					error = "Duplicate of anonymous private record (RecordID: R"+ep.duplicate.id+")";
			
			

			ep.ep_status = ExperimentalProperty.STATUS_ERROR;
			ep.errorComment = error;
			ep.updateHash();
			
			preview.duplicate = ep.duplicate;
			preview.status = PreviewRecordStatus.duplicate_external;
			if(preview.action != PreviewUploadAction.save) preview.action = PreviewUploadAction.put_original_to_basket;

			ep.duplicate = null;
		}
	}
}
