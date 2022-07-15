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
import qspr.metaserver.configurations.DLCAConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration.ScalingType;
import qspr.modelling.configurations.CDSConfiguration;

public class DLCAConfigurator extends DescriptorsConfigurator 
{

	public DLCAConfigurator()
	{
		dataPreprocessingStep = true;
		firstConfigurationStep = "dlca";
	}

	@Override
	public void setModel(Model model)
	{
		super.setModel(model);
		((CDSConfiguration)model.attachment.getObject().configuration).modelConfiguration = new DLCAConfiguration();
		((CDSConfiguration)model.attachment.getObject().configuration).selection = null;
		this.stepAfterDescriptors = "dlca";
	}

	public WebModel dlca()
	{
		return new WebModel(getDefaultTemplate()).setTemplate("modeller/configurators/dlca");
	}

	public void dlcaSubmit(HttpServletRequest request)
	{
		DLCAConfiguration conf = (DLCAConfiguration) 
				((CDSConfiguration)model.attachment.getObject().configuration).modelConfiguration;

		conf.epochs = Integer.valueOf(request.getParameter("epochs"));
		conf.batch = Integer.valueOf(request.getParameter("batch"));
		conf.scaleTypeY = ScalingType.RANGE;

		currentPage = startStep;
	}

}