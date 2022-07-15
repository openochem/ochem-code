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

package qspr.schedule;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.OCHEMConfiguration;
import qspr.util.WrapperThread;

import com.eadmet.mmpa.MMPIndexingService;

@DatabaseMaintenanceJob
public class MMPFetchIndexTask extends OchemCronjobTask
{

	private static final Logger logger = LogManager.getLogger(MMPFetchIndexTask.class);


	@Override
	public void executeTask() throws Exception
	{
		new WrapperThread()
		{
			@Override
			public void wrapped() throws Exception
			{
				logger.info(" fetching started");
				MMPIndexingService.getInstance().fetchIndexingTasks();
			}
		}.run();
	}

	@Override
	protected boolean shouldRun()
	{
		if (!OCHEMConfiguration.mmpSimilarityIndexingTask) return false;

		return super.shouldRun();
	}

	public static void main(String[] args) throws Exception
	{
		new MMPFetchIndexTask().executeTask();
	}
}
