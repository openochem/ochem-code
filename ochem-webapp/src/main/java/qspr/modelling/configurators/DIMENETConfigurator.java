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

import java.io.IOException;

import javax.servlet.http.HttpServletRequest;

import qspr.entities.Model;
import qspr.frontend.WebModel;
import qspr.metaserver.configurations.DIMENETConfiguration;
import qspr.modelling.configurations.CDSConfiguration;

public class DIMENETConfigurator extends DescriptorsConfigurator{

	@Override
	public void setModel(Model model)
	{
		super.setModel(model);
		CDSConfiguration cds = (CDSConfiguration)model.attachment.getObject().configuration;
		cds.modelConfiguration = new DIMENETConfiguration();
		this.stepAfterDescriptors = "dimenet";
	}

	public DIMENETConfigurator()
	{
		dataPreprocessingStep = true;
		firstConfigurationStep = "dimenet";
	}

	public WebModel dimenet()
	{
		return new WebModel(getDefaultTemplate()).setTemplate("modeller/configurators/dimenet");
	}

	public void dimenetSubmit(HttpServletRequest request) throws IOException
	{
		CDSConfiguration cds = (CDSConfiguration)model.attachment.getObject().configuration;

		DIMENETConfiguration conf = (DIMENETConfiguration) 
				cds.modelConfiguration;

		conf.nbepochs = Integer.valueOf(request.getParameter("nbepochs"));
		conf.batch = Integer.valueOf(request.getParameter("batch"));
		conf.early = request.getParameter("early") != null;
		if(conf.batch > 64) conf.batch = 64;

		cds.optimisationConfiguration=DescriptorsConfigurator.configureStructureOptimisation(request);
		if(cds.optimisationConfiguration != null) conf.external3D = true;

		currentPage = startStep;
	}
}