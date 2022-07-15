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

package qspr.modelling.ad;

import java.util.ArrayList;
import java.util.List;

import qspr.Globals;
import qspr.entities.ModelMapping;
import qspr.metaserver.configurations.ADConfiguration;
import qspr.modelling.SetStatistics;

public abstract class ADFactory
{
	protected SetStatistics statsWithDM;
	protected SetStatistics statsWithPredictions;
	protected ADConfiguration adConfiguration;
	protected ADConfiguration initialConf;
	protected String dm;
	protected ModelMapping mm;
	public boolean considerPredicates;

	public static ThreadLocal<Integer> thresholdPercentage = new ThreadLocal<Integer>();

	public ADFactory(SetStatistics statistics, SetStatistics dmSource, ModelMapping mm, String dm, ADConfiguration initialConf)
	{
		if (initialConf == null)
			initialConf = new ADConfiguration(); // Defaults
		if (dmSource == null)
			dmSource = statistics;
		this.initialConf = initialConf;
		adConfiguration = new ADConfiguration();
		adConfiguration.intervals = new ArrayList<Double>();
		adConfiguration.errors = new ArrayList<Double>();
		adConfiguration.percents = new ArrayList<Double>();
		adConfiguration.epIds = new ArrayList<Long>();
		adConfiguration.averagingType = initialConf.averagingType;
		adConfiguration.dmName = dm;

		int percentage = thresholdPercentage.get() != null ? thresholdPercentage.get() : 95;

		// Identify the DM threshold for AD as the one embracing 95% of the set compounds
		List<Integer> order = dmSource.getDMOrder(dm);
		if (order.size() > 0)
			adConfiguration.threshold = dmSource.getDM(order.get((order.size() - 1) * percentage / 100), dm);
		else
			adConfiguration.threshold = 0D;

		this.statsWithDM = dmSource;
		this.statsWithPredictions = statistics;
		this.mm = mm;
		this.dm = dm;

		considerPredicates = Globals.considerPredicates();
	}

	public static ADFactory getInstance(SetStatistics statistics, SetStatistics dmSource, ModelMapping mm, String dm, ADConfiguration initialConf)
	{
		if (initialConf == null)
			initialConf = new ADConfiguration();
		if (initialConf.averagingType == null)
			initialConf.averagingType = mm.property.isNumeric() ? "mgd" : "ma";
		if (dmSource == null)
			dmSource = statistics;

		if ("cumulative".equals(initialConf.averagingType))
			return new CumulativeADFactory(statistics, dmSource, mm, dm, initialConf);
		else if ("mgd".equals(initialConf.averagingType))
			return new StepsADFactory(statistics, dmSource, mm, dm, initialConf);
		else
			return new SlidingWindowADFactory(statistics, dmSource, mm, dm, initialConf);
	}

	abstract public ADConfiguration getAD();
}