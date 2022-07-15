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

import qspr.Globals;
import qspr.entities.Model;
import qspr.frontend.WebModel;
import qspr.metaserver.configurations.HAMNETConfiguration;
import qspr.modelling.configurations.CDSConfiguration;

public class HAMNETConfigurator extends DescriptorsConfigurator {

	@Override
	public void setModel(Model model)
	{
		super.setModel(model);
		((CDSConfiguration)model.attachment.getObject().configuration).modelConfiguration = new HAMNETConfiguration();
		((CDSConfiguration)model.attachment.getObject().configuration).selection = null;
		this.stepAfterDescriptors = "hamnet";
	}

	public HAMNETConfigurator()
	{
		firstConfigurationStep = "hamnet";
	}

	public WebModel hamnet()
	{
		return new WebModel(getDefaultTemplate()).setTemplate("modeller/configurators/hamnet");
	}

	public void hamnetSubmit(HttpServletRequest request)
	{
		HAMNETConfiguration conf = (HAMNETConfiguration) 
				((CDSConfiguration)model.attachment.getObject().configuration).modelConfiguration;

		conf.learningRate = Double.valueOf(request.getParameter("learningRate"));
		conf.batch_size= Integer.valueOf(request.getParameter("batch_size"));
		conf.refine = request.getParameter("refine") != null;
		conf.dropout = Double.valueOf(request.getParameter("dropout"));

		addStandard(conf, request, Globals.isOCHEMDeveloper()?1000:100);

		currentPage = startStep;
	}

}
