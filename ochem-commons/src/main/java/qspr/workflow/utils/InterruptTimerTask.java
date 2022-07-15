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

package qspr.workflow.utils;
import java.util.TimerTask;

/**
 * Task will interrupt unless there was a message from Server that it is still alive
 * @author itetko
 *
 */

public class InterruptTimerTask extends TimerTask
{
	private Thread thread;

	private boolean skipInterruption = false;

	public InterruptTimerTask(Thread t)
	{
		this.thread = t;
	}

	public void run()
	{
		if(skipInterruption)skipInterruption = false;
		else
			thread.interrupt();
	}

	public void resetInterruption() {
		skipInterruption = true;
	}
}