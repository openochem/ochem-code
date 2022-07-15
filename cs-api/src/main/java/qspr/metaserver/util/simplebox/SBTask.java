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

package qspr.metaserver.util.simplebox;

import java.io.IOException;
import java.io.Writer;
import java.util.HashMap;
import java.util.LinkedHashMap;

import org.json.simple.JSONValue;

/**
 * 
 * 
 * 
 * @author Pantelis Sopasakis
 */

public class SBTask extends AbstractComponent {

	/**
	 * 
	 * Status of the Task
	 */

	public enum TaskStatus {
		COMPLETED, RUNNING, QUEUED, REJECTED, ERROR, CANCELLED
	}

	// Private Fields:
	private String resultURI;
	private int httpStatus = -1;
	private long duration = -1;
	private long expectedDuration = -1;
	private long creationTimestamp = System.currentTimeMillis();
	private ErrorReport error;
	private TaskStatus taskStatus = TaskStatus.QUEUED;

	public SBTask() {
		super();
	}

	// <editor-fold defaultstate="collapsed" desc="Getters and Setters">
	public long getExpectedDuration() {
		return expectedDuration;
	}

	public void setExpectedDuration(long expectedDuration) {
		this.expectedDuration = expectedDuration;
	}

	public long getCreationTimestamp() {
		return creationTimestamp;
	}

	public long getDuration() {
		return duration;
	}

	public void setDuration(long duration) {
		this.duration = duration;
	}

	public ErrorReport getError() {
		return error;
	}

	public void setError(ErrorReport error) {
		this.error = error;
	}

	public int getHttpStatus() {
		return httpStatus;
	}

	public void setHttpStatus(int httpStatus) {
		this.httpStatus = httpStatus;
	}

	public String getResultURI() {
		return resultURI;
	}

	public void setResultURI(String resultURI) {
		this.resultURI = resultURI;
	}

	public TaskStatus getTaskStatus() {
		return taskStatus;
	}

	public void setTaskStatus(TaskStatus taskStatus) {
		this.taskStatus = taskStatus;
	}

	// </editor-fold>

	public void writeJSONString(Writer writer) throws IOException {
		HashMap<String, Object> obj = new LinkedHashMap<String, Object>();
		obj.put("object_type", "Task");
		if (taskStatus != null) {
			obj.put("taskStatus", taskStatus.toString());
		}

		if (duration > 0) {
			obj.put("duration", Long.valueOf(duration));
		}

		if (expectedDuration > 0) {
			obj.put("exp_duration", Long.valueOf(expectedDuration));
		}

		if (httpStatus > 0) {
			obj.put("httpStatus", Integer.valueOf(httpStatus));
		}

		if (resultURI != null) {
			obj.put("resultURI", resultURI);
		}

		if (error != null) {
			obj.put("error", error);
		}

		if (getMeta() != null) {
			obj.put("metadata", getMeta());
		}

		obj.put("creationTimestamp", Long.valueOf(creationTimestamp));
		JSONValue.writeJSONString(obj, writer);
	}

}
