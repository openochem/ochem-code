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
import qspr.metaserver.configurations.CNFConfiguration;
import qspr.metaserver.configurations.CNFConfiguration.TOKENIZER;
import qspr.modelling.configurations.CDSConfiguration;

public class CNFConfigurator extends DescriptorsConfigurator{

	@Override
	public void setModel(Model model)
	{
		super.setModel(model);
		((CDSConfiguration)model.attachment.getObject().configuration).modelConfiguration = new CNFConfiguration();
		this.stepAfterDescriptors = "cnf";
	}


	// firstConfigurationStep indicates to skip descriptors
	public CNFConfigurator()
	{
		dataPreprocessingStep = true;
		firstConfigurationStep = "cnf";
	}

	public WebModel cnf()
	{
		return new WebModel(getDefaultTemplate()).setTemplate("modeller/configurators/cnf");
	}

	public void cnfSubmit(HttpServletRequest request) throws IOException
	{
		CNFConfiguration conf = (CNFConfiguration) 
				((CDSConfiguration)model.attachment.getObject().configuration).modelConfiguration;

		redValues(conf, request);
		currentPage = startStep;
	}

	static public void redValues(CNFConfiguration conf, HttpServletRequest request) throws IOException {
		addStandard(conf, request,10000);
		if(conf.augmentation != null && conf.augmentation < -1) conf.augmentation = -1;
		//conf.setOriginalTaxonomy(request.getParameter("taxonomy"));

		conf.FPDim = Integer.valueOf(request.getParameter("FPDim"));
		conf.nLayers = Integer.valueOf(request.getParameter("nLayers"));
		if (request.getParameter("FFNET_dim") != null)
			conf.FFNET_dim = (request.getParameter("FFNET_dim"));
		conf.rate = Double.valueOf(request.getParameter("rate"));
		conf.batch = Integer.valueOf(request.getParameter("batch"));

		conf.type = request.getParameter("fingerprint_function") == null ? null : 
			CNFConfiguration.CNFTYPE.valueOf(request.getParameter("fingerprint_function"));

		/*
		conf.normalisation = request.getParameter("normalisation") == null ? null : 
			CNFConfiguration.NORMALIZE.valueOf(request.getParameter("normalisation"));

		if (request.getParameter("normy") != null) {
			String normY = request.getParameter("normy");
			if (!normY.isEmpty())
				conf.scaleTypeY = ScalingType.valueOf(normY);
		}

		conf.errorFunction = request.getParameter("error_function") == null ? null : 
			CNFConfiguration.ERROR_FUNCTION.valueOf(request.getParameter("error_function"));

		//if(conf.isUsesFilter()) conf.filterSize = Integer.valueOf(request.getParameter("filter_size"));
		//else
		//	conf.filterSize = null;

		if(conf.isUsesAlpha()) {
			conf.nfp_alpha = Double.valueOf(request.getParameter("nfp_alpha"));
			if(conf.nfp_alpha == null || conf.nfp_alpha < 0 || conf.nfp_alpha > 1 )throw new IOException("Alpha should be in the range [0,1] but it is :" + conf.nfp_alpha);			
		}else
			conf.nfp_alpha = null;
		 */

		//conf.highway =request.getParameter("highway") == null ? null : 
		//	Integer.valueOf(request.getParameter("highway"));
		//if(conf.highway != null && conf.highway < 1) conf.highway = null;

		conf.tokenizer = request.getParameter("tokenizer") == null ? null : 
			CNFConfiguration.TOKENIZER.valueOf(request.getParameter("tokenizer"));

		if(conf.tokenizer != null && conf.tokenizer == TOKENIZER.CHAR) conf.tokenizer = null; // the default type

		if(request.getParameter("dropout") != null) {
			conf.dropout = Double.valueOf(request.getParameter("dropout"));
			if(conf.dropout<0)conf.dropout=-1.;
			if(conf.dropout>0.5)conf.dropout=0.5;
		}		

	}


}
