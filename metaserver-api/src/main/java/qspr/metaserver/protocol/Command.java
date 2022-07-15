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

package qspr.metaserver.protocol;

import java.io.Serializable;

public class Command implements Serializable
{
	public static final long serialVersionUID = 2;
	public Integer id;
	public Serializable data;
	public String senderId;

	// Common commands
	public static final int FAILURE = -1;
	public static final int DENIAL = -2;
	public static final int KILL_TASK = -3;

	// Client >>> MetaServer
	public static final int  CL_SUBMIT_TASK = 1;
	public static final int  CL_QUERY_TASK = 2;
	public static final int  CL_QUERY_TASK_STATUS = 3;
	public static final int  CL_QUERY_READY_TASK = 4;
	public static final int  CL_DELETE_TASK = 5;
	public static final int  CL_DELETE_CHILDREN = 6;
	public static final int  CL_REGISTER_ADMIN_IP = 7;
	public static final int  CL_SET_PRIORITY = 8;
	public static final int  CL_GET_SUPPORTED_TASKS = 9;
	public static final int  CL_GET_TASK_BY_MD5 = 10;
	public static final int  CL_SET_PARENT = 11;
	public static final int  CL_GET_TASKS_SUMMARY = 12;
	public static final int  CL_GET_PARENT_ID = 13;

	// MetaServer >>> Client & CS
	public static final int  MS_UNKNOWN_TASK = 0x10;
	public static final int  MS_TASK_STATUS = 0x20;
	public static final int  MS_TERMINATE = 0x30;
	public static final int  MS_RESTART = 0x40;
	public static final int  MS_UPDATEREQUIRED = 0x50;
	public static final int  MS_TASKFOUND = 0x60;
	public static final int  MS_COMMAND = 0x70;
	public static final int  MS_GET_LOGS = 0x80;
	public static final int  MS_STOP = 0x90;

	// CS >>> MetaServer
	public static final int  CS_TASK_CALCULATED = 0x100;
	public static final int  CS_REGISTER = 0x300;
	public static final int  CS_QUERYINFO = 0x200;


	// MetaServer >> CS
	public static final int  MS_ASSIGN_TASK = 0x2000;
	public static final int  MS_OK = 0;
	public static final String MEASERVER_DOWN = "Can't fetch a task since the metaserver is down. Retry in a while...";

	public static final String UNKNOWN_TASK = "Unknown task ID:";

	static public String LOCAL = "local";
	static public String FIXED = "mongofix";
	public static final String VERSIONTEMPLATEXML = "version-template.xml";

	public Command(Integer id, Serializable data)
	{
		this.id = id;
		this.data = data;
	}

	public Command sid(String sid)
	{
		this.senderId = sid;
		return this;
	}

	public String toString()
	{
		switch (this.id)
		{
		case KILL_TASK:
			return "KILL_TASK";
		case CL_SUBMIT_TASK:
			return "CL_SUBMIT_TASK";
		case CL_QUERY_TASK:
			return "CL_QUERY_TASK";
		case CL_QUERY_TASK_STATUS:
			return "CL_QUERY_TASK_STATUS";
		case CL_QUERY_READY_TASK:
			return "CL_QUERY_READY_TASK";
		case CL_SET_PARENT:
			return "CL_SET_PARENT";
		case CL_REGISTER_ADMIN_IP:
			return "CL_REGISTER_ADMIN_IP";
		case CS_TASK_CALCULATED:
			return "CS_TASK_CALCULATED";
		case CS_REGISTER:
			return "CS_REGISTER";
		case CS_QUERYINFO:
			return "CS_QUERYINFO";
		case MS_ASSIGN_TASK:
			return "MS_ASSIGN_TASK";
		case MS_TASK_STATUS:
			return "MS_TASK_STATUS";
		case FAILURE:
			return "FAILURE";
		case DENIAL:
			return "DENIAL";
		case CL_GET_SUPPORTED_TASKS:
			return "CL_GET_SUPPORTED_TASKS";
		case MS_UNKNOWN_TASK:
			return "MS_UNKNOWN_TASK";
		case CL_SET_PRIORITY:
			return "CL_SET_PRIORITY";
		case MS_TASKFOUND:
			return "MS_TASKFOUND";
		case CL_GET_TASK_BY_MD5:
			return "CL_GET_TASK_BY_MD5";
		case CL_GET_TASKS_SUMMARY:
			return "CL_GET_TASKS_SUMMARY";
		case MS_COMMAND:
			return "MS_COMMAND";

		default:
			return "Unknown command "+this.id;
		}
	}
}
