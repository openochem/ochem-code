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
import qspr.metaserver.configurations.ATTFPConfiguration;
import qspr.modelling.configurations.CDSConfiguration;

public class ATTFPConfigurator extends DescriptorsConfigurator {

	@Override
	public void setModel(Model model)
	{
		super.setModel(model);
		((CDSConfiguration)model.attachment.getObject().configuration).modelConfiguration = new ATTFPConfiguration();
		((CDSConfiguration)model.attachment.getObject().configuration).selection = null;
		this.stepAfterDescriptors = "attfp";
	}

	public ATTFPConfigurator()
	{
		firstConfigurationStep = "attfp";
	}

	public WebModel attfp()
	{
		return new WebModel(getDefaultTemplate()).setTemplate("modeller/configurators/attfp");
	}

	public void attfpSubmit(HttpServletRequest request)
	{
		ATTFPConfiguration conf = (ATTFPConfiguration) 
				((CDSConfiguration)model.attachment.getObject().configuration).modelConfiguration;

		conf.nepochs = Integer.valueOf(request.getParameter("nepochs"));
		conf.batch = Integer.valueOf(request.getParameter("batch"));
		conf.nepochs = conf.nepochs == null || conf.nepochs > 1000 ? 1000: conf.nepochs <1? 1: conf.nepochs;
		conf.patience_reduce = Integer.valueOf(request.getParameter("patience_reduce"));
		conf.patience_early = Integer.valueOf(request.getParameter("patience_early"));
		conf.radius = Integer.valueOf(request.getParameter("radius"));
		conf.T = Integer.valueOf(request.getParameter("T"));
		conf.fp_dim = Integer.valueOf(request.getParameter("fp_dim"));
		conf.cosineT = Integer.valueOf(request.getParameter("cosineT"));
		conf.dropout = Double.valueOf(request.getParameter("dropout"));
		conf.lr = Double.valueOf(request.getParameter("lr"));
		conf.weight_decay = Double.valueOf(request.getParameter("weight_decay"));
		conf.cosine = request.getParameter("cosine") != null;
		conf.early = request.getParameter("early") != null;
		conf.singleT = request.getParameter("singleT") != null?true:null;
		conf.simpleO = request.getParameter("simpleO") != null?true:null;
		conf.lngru = request.getParameter("lngru") != null?true:null;
		conf.molsize = Integer.valueOf(request.getParameter("molsize"));
		currentPage = startStep;
	}

}
