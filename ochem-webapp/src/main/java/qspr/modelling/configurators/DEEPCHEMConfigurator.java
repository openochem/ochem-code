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
import qspr.metaserver.configurations.DeepChemConfiguration;
import qspr.metaserver.configurations.DeepChemConfiguration.DeepChemMethod;
import qspr.modelling.configurations.CDSConfiguration;

public class DEEPCHEMConfigurator extends DescriptorsConfigurator 
{

	public DEEPCHEMConfigurator()
	{
		dataPreprocessingStep = true;
		firstConfigurationStep = "deepchem";
	}

	@Override
	public void setModel(Model model)
	{
		super.setModel(model);
		((CDSConfiguration)model.attachment.getObject().configuration).modelConfiguration = new DeepChemConfiguration();
		((CDSConfiguration)model.attachment.getObject().configuration).selection = null;
		this.stepAfterDescriptors = "deepchem";
	}

	public WebModel deepchem()
	{
		return new WebModel(getDefaultTemplate()).setTemplate("modeller/configurators/deepchem");
	}

	public void deepchemSubmit(HttpServletRequest request)
	{
		DeepChemConfiguration conf = (DeepChemConfiguration) 
				((CDSConfiguration)model.attachment.getObject().configuration).modelConfiguration;

		conf.method = DeepChemMethod.valueOf(request.getParameter("method"));
		//conf.learning_rate = Double.valueOf(request.getParameter("learning_rate"));
		//conf.dropout = Double.valueOf(request.getParameter("dropout"));
		//conf.dense_layer_size = Integer.valueOf(request.getParameter("dense_layer_size"));
		//conf.M = Integer.valueOf(request.getParameter("M"));
		//conf.T = Integer.valueOf(request.getParameter("T"));
		//conf.n_hidden = Integer.valueOf(request.getParameter("n_hidden"));
		//conf.n_embedding = Integer.valueOf(request.getParameter("n_embedding"));
		//conf.graph_conv_layers = request.getParameter("graph_conv_layers");

		addStandard(conf, request, 200);
		if(conf.method  != DeepChemMethod.TEXTCNN)
			conf.setAugmentations(1, 1, false);

		currentPage = startStep;
	}

}