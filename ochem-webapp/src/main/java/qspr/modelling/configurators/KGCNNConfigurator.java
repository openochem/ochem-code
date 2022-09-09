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
import qspr.metaserver.configurations.KGCNNConfiguration;
import qspr.metaserver.configurations.KGCNNConfiguration.KGCNN;
import qspr.modelling.configurations.CDSConfiguration;

public class KGCNNConfigurator extends DescriptorsConfigurator {

	@Override
	public void setModel(Model model)
	{
		super.setModel(model);
		((CDSConfiguration)model.attachment.getObject().configuration).modelConfiguration = new KGCNNConfiguration();
		((CDSConfiguration)model.attachment.getObject().configuration).selection = null;
		this.stepAfterDescriptors = "kgcnn";
	}

	public KGCNNConfigurator()
	{
		firstConfigurationStep = "kgcnn";
	}

	public WebModel kgcnn()
	{
		return new WebModel(getDefaultTemplate()).setTemplate("modeller/configurators/kgcnn");
	}

	public void kgcnnSubmit(HttpServletRequest request)
	{
	
		CDSConfiguration cds = (CDSConfiguration)model.attachment.getObject().configuration;

		KGCNNConfiguration conf = (KGCNNConfiguration) 
				((CDSConfiguration)model.attachment.getObject().configuration).modelConfiguration;

		if(request.getParameter("sanitize") != null)conf.setSanitize();

		conf.method = KGCNN.valueOf(request.getParameter("method"));
		conf.nepochs = Integer.valueOf(request.getParameter("nepochs"));
		conf.batch = Integer.valueOf(request.getParameter("batch"));
		conf.nepochs = conf.nepochs == null || conf.nepochs > 10000 ? 10000: conf.nepochs <1? 1: conf.nepochs;

		cds.optimisationConfiguration=DescriptorsConfigurator.configureStructureOptimisation(request);
		conf.setUse3D(cds.optimisationConfiguration != null);

		currentPage = startStep;
	}

}
