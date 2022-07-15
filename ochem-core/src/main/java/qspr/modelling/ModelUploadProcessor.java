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

package qspr.modelling;

import java.sql.Timestamp;
import java.util.Date;

import qspr.metaserver.protocol.Task;
import qspr.workflow.datatypes.WorkflowNodeData;

import com.eadmet.exceptions.UserFriendlyException;

public class ModelUploadProcessor extends CDSModelProcessor
{
	//	private static transient final Logger logger = Logger.getLogger(ModelUploadProcessor.class);

	public boolean methodNeedsStructures = false;

	WorkflowNodeData uploadedWnd;
	public boolean prepared = false;

	public ModelUploadProcessor()
	{
	}

	@Override
	protected Task createTeacherTask() throws Exception
	{
		throw new UserFriendlyException("Can not train the uploaded model "+model.name);
	}

	@Override
	public Task createApplierTask() throws Exception
	{
		throw new UserFriendlyException("Can not apply the uploaded model "+model.name);
	}

	public void setUploadedData(WorkflowNodeData uploadedWnd)
	{
		this.uploadedWnd = uploadedWnd;
	}

	public ModelUploadProcessor(WorkflowNodeData uploadedWnd)
	{
		this.uploadedWnd = uploadedWnd;
	}

	@Override
	protected void postTeacherTask() throws Exception
	{
		model.taskId = null;
	}

	@Override
	public Task waitForTask() throws Exception
	{
		Task t = new Task();
		t.setResult(uploadedWnd);
		t.timeAssigned = new Timestamp(new Date().getTime());
		t.timeCompleted = new Timestamp(new Date().getTime());
		t.status = "ready";
		return t;
	}

	@Override
	public void prepare() throws Exception
	{
		if (prepared)
			return;

		super.prepare();
		prepared = true;
	}
}
