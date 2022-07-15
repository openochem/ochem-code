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

package qspr.metaserver.serv;

/**
 * A set of constants that define the exchange protocol between MultiServer and ServerRunner
 * 
 * @author midnighter
 *
 */
public class ExchangeCommands 
{
	public static final String MULTISERVER_COMMAND_PREFIX = ">>";
	public static final String MULTISERVER_PING = MULTISERVER_COMMAND_PREFIX +"ping";
	public static final String MULTISERVER_TERMINATE = MULTISERVER_COMMAND_PREFIX +"terminate";
	public static final String MULTISERVER_STARTED = MULTISERVER_COMMAND_PREFIX +"task started";
	public static final String MULTISERVER_RESTART = MULTISERVER_COMMAND_PREFIX +"restart me";
	public static final String MULTISERVER_RETEST = MULTISERVER_COMMAND_PREFIX +"retest me";
	public static final String MULTISERVER_UPDATEREQUIRED = MULTISERVER_COMMAND_PREFIX + "update required";
	public static final String MULTISERVER_TESTS_FILE_NAME="/tests.txt";
	
	public static final String SERVER_RESTART_COMMAND="RESTART";
	public static final String SERVER_TERMINATE_COMMAND="TERMINATE";
	public static final long   SERVER_NO_ACTIVITY_MAX_TIME = 30 * 60 * 1000; // 30 minutes
}
