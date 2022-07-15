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

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.util.WrapperThread;


@DatabaseMaintenanceJob
public class MMPSimilarityIndexingTask extends OchemCronjobTask
{
	private static final Logger logger = LogManager.getLogger(MMPSimilarityIndexingTask.class);
	@Override
	public void executeTask() throws Exception
	{
		new WrapperThread()
		{

			@Override
			public void wrapped() throws Exception
			{

				Globals.session().createSQLQuery("update MMPair SET similarity = (select similarity(a.screen,b.screen) from  StructureQuery as a, StructureQuery as b where similarity is NULL and a.mapping2_id=mol1 and b.mapping2_id=mol2)").list(); 
				logger.info("MMPair were updated");
				Globals.restartAllTransactions(true);
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
		new MMPSimilarityIndexingTask().executeTask();
	}
}
