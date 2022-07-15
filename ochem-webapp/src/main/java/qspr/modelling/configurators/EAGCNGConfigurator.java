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

import com.eadmet.exceptions.UserFriendlyException;

import qspr.entities.Model;
import qspr.frontend.WebModel;
import qspr.metaserver.configurations.CNFConfiguration;
import qspr.metaserver.configurations.EAGCNGConfiguration;
import qspr.modelling.configurations.CDSConfiguration;

public class EAGCNGConfigurator extends DescriptorsConfigurator{

	public EAGCNGConfigurator()
	{
		dataPreprocessingStep = true;
		firstConfigurationStep = "eagcng";
	}

	@Override
	public void setModel(Model model)
	{
		super.setModel(model);
		((CDSConfiguration)model.attachment.getObject().configuration).modelConfiguration = new EAGCNGConfiguration();
		this.stepAfterDescriptors = "eagcng";
	}

	public WebModel eagcng()
	{
		return new WebModel(getDefaultTemplate()).setTemplate("modeller/configurators/eagcng");
	}

	public void eagcngSubmit(HttpServletRequest request)
	{
		EAGCNGConfiguration conf = (EAGCNGConfiguration) 
				((CDSConfiguration)model.attachment.getObject().configuration).modelConfiguration;

		addStandard(conf, request, 1000);

		//conf.seed = Integer.valueOf(request.getParameter("seed"));
		conf.batchsize = Integer.valueOf(request.getParameter("batchsize"));
		conf.molsize = Integer.valueOf(request.getParameter("molsize"));
		conf.dropout = Double.valueOf(request.getParameter("dropout"));
		conf.learningRate = Double.valueOf(request.getParameter("learningRate"));
		conf.weightDecay = Double.valueOf(request.getParameter("weightDecay"));
		conf.method = request.getParameter("method");
		conf.n_sgc1 = request.getParameter("n_sgc1"); if(conf.n_sgc1!=null)conf.n_sgc1=conf.n_sgc1.trim();
		conf.n_sgc2 = request.getParameter("n_sgc2"); if(conf.n_sgc2!=null)conf.n_sgc2=conf.n_sgc2.trim();
		if(conf.n_sgc2 != null && conf.n_sgc2.length() == 0)conf.n_sgc2 = null;
		conf.n_sgc3 = request.getParameter("n_sgc3"); if(conf.n_sgc3!=null)conf.n_sgc3=conf.n_sgc3.trim();
		if(conf.n_sgc3 != null && conf.n_sgc3.length() == 0)conf.n_sgc3 = null;
		conf.n_den = request.getParameter("n_den");
		conf.normalisation = request.getParameter("normalisation") == null ? null : 
			CNFConfiguration.NORMALIZE.valueOf(request.getParameter("normalisation"));
		conf.gate = request.getParameter("gate") != null;

		if(conf.gate) {
			if(!conf.n_sgc1.equalsIgnoreCase(conf.n_sgc2) && conf.n_sgc3 == null ) 
				throw new UserFriendlyException("For gate use n_sgc1 == n_sgc2");

			if(conf.n_sgc3 != null && (!conf.n_sgc1.equalsIgnoreCase(conf.n_sgc2) || !conf.n_sgc1.equalsIgnoreCase(conf.n_sgc3))) 
				throw new UserFriendlyException("For gate use n_sgc1 == n_sgc2 == n_sgc3");
		}

		//addAugmentation(conf, request);

		currentPage = startStep;
	}
}
