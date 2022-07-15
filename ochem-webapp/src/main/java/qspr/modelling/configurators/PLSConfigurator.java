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
import qspr.metaserver.configurations.PLSConfiguration;
import qspr.modelling.configurations.CDSConfiguration;

public class PLSConfigurator extends DescriptorsConfigurator 
{
	@Override
	public void setModel(Model model)
	{
		super.setModel(model);
		((CDSConfiguration)model.attachment.getObject().configuration).modelConfiguration = new PLSConfiguration();
		this.stepAfterDescriptors = "pls";
	}

	public WebModel pls()
	{
		return new WebModel(getDefaultTemplate()).setTemplate("modeller/configurators/pls");
	}

	public void plsSubmit(HttpServletRequest request)
	{
		PLSConfiguration conf = (PLSConfiguration) 
				((CDSConfiguration)model.attachment.getObject().configuration).modelConfiguration;

		String numLatentVariables = request.getParameter("num-latent-vars");
		if (numLatentVariables!=null)
			conf.numLatentVariables = Integer.valueOf(numLatentVariables);
		else
			conf.numLatentVariables = 0;

		if (request.getParameter("limitrange") != null)
			conf.limitRange = true;

		currentPage = startStep;
	}
}
