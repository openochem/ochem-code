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
import qspr.metaserver.configurations.IMGQSARConfiguration;
import qspr.modelling.configurations.CDSConfiguration;

public class IMGQSARConfigurator extends DescriptorsConfigurator {

	@Override
	public void setModel(Model model)
	{
		super.setModel(model);
		((CDSConfiguration)model.attachment.getObject().configuration).modelConfiguration = new IMGQSARConfiguration();
		((CDSConfiguration)model.attachment.getObject().configuration).selection = null;
		this.stepAfterDescriptors = "imgqsar";
	}

	public IMGQSARConfigurator()
	{
		firstConfigurationStep = "imgqsar";
	}

	public WebModel imgqsar()
	{
		return new WebModel(getDefaultTemplate()).setTemplate("modeller/configurators/imgqsar");
	}

	public void imgqsarSubmit(HttpServletRequest request)
	{
		IMGQSARConfiguration conf = (IMGQSARConfiguration) 
				((CDSConfiguration)model.attachment.getObject().configuration).modelConfiguration;

		conf.batch_size = Integer.valueOf(request.getParameter("batch_size"));
		conf.learning_rate = Double.valueOf(request.getParameter("learning_rate"));
		//conf.backbone = request.getParameter("backbone");
		
		addStandard(conf, request, 1000);
		
		currentPage = startStep;
	}

}
