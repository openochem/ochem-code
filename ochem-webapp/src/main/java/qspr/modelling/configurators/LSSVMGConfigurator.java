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
import qspr.metaserver.configurations.LSSVMGConfiguration;
import qspr.modelling.configurations.CDSConfiguration;

public class LSSVMGConfigurator extends DescriptorsConfigurator{

	@Override
	public void setModel(Model model)
	{
		super.setModel(model);
		((CDSConfiguration)model.attachment.getObject().configuration).modelConfiguration = new LSSVMGConfiguration();
		this.stepAfterDescriptors = "lssvmg";
	}

	public WebModel lssvmg()
	{
		return new WebModel(getDefaultTemplate()).setTemplate("modeller/configurators/lssvmg");
	}

	public void lssvmgSubmit(HttpServletRequest request)
	{
		LSSVMGConfiguration conf = (LSSVMGConfiguration) 
				((CDSConfiguration)model.attachment.getObject().configuration).modelConfiguration;
		conf.kernel = request.getParameter("kernel-lssvmg");
		conf.cv = Integer.parseInt(request.getParameter("cv-lssvmg"));
		conf.useGlobaL = (request.getParameter("use-global") != null);
		if(!conf.useGlobaL)conf.useGlobaL = null;
		if (request.getParameter("additionalParam") != null) {
			conf.additionalParam = request.getParameter("additionalParam");
			if(conf.additionalParam.length()<3)
				conf.additionalParam = null;
		}

		currentPage = startStep;
	}
}
