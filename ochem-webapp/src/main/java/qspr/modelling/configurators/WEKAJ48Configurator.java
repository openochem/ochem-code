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

package qspr.modelling.configurators;

import javax.servlet.http.HttpServletRequest;

import qspr.entities.Model;
import qspr.frontend.WebModel;
import qspr.metaserver.configurations.J48Configuration;
import qspr.modelling.configurations.CDSConfiguration;

public class WEKAJ48Configurator extends DescriptorsConfigurator 
{
	@Override
	public void setModel(Model model) 
	{
		super.setModel(model);
		((CDSConfiguration) model.attachment.getObject().configuration).modelConfiguration = new J48Configuration();
		this.stepAfterDescriptors = "j48";
	}

	public WebModel j48() 
	{
		return new WebModel(getDefaultTemplate()).setTemplate("modeller/configurators/j48");
	}

	public void j48Submit(HttpServletRequest request) 
	{
		J48Configuration conf = (J48Configuration) ((CDSConfiguration) model.attachment.getObject().configuration).modelConfiguration;
		conf.useUnpruned = (request.getParameter("use-unpruned") != null);
		if (conf.useUnpruned)
		{
			conf.confidence = null;
			conf.instancesPerLeaf = null;
			conf.reducedPruning = null;
			conf.numReducedPruningFolds = null;
			conf.useBinarySplits = null;
		} else
		{
			conf.instancesPerLeaf = Integer.valueOf(request.getParameter("instances-per-leaf"));
			conf.reducedPruning = (request.getParameter("reduced-pruning") != null);

			if (conf.reducedPruning)
			{
				conf.numReducedPruningFolds = Integer.valueOf(request.getParameter("num-reduced-pruning-folds"));
				conf.confidence = null;
			} else {
				conf.confidence = Double.valueOf(request.getParameter("confidence"));
				if(conf.confidence>0.999) conf.confidence = 0.999;
				else
					if(conf.confidence<0.001) conf.confidence = 0.001;
			}

		}
		conf.useBinarySplits = (request.getParameter("use-binary-splits") != null);
		conf.dontPerformRaising = (request.getParameter("dont-perform-raising") != null);
		conf.noCleanup = (request.getParameter("no-cleanup") != null);
		conf.useLaplaseSmoothing = (request.getParameter("use-laplase-smoothing") != null);
		conf.seed = Integer.valueOf(request.getParameter("seed"));
		currentPage = startStep;
	}
}
