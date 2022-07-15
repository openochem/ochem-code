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

import java.util.ArrayList;
import java.util.List;

import javax.servlet.http.HttpServletRequest;

import qspr.Globals;
import qspr.entities.Basket;
import qspr.entities.Model;
import qspr.entities.PropertyOption;
import qspr.frontend.WebModel;
import qspr.metaserver.configurations.LibSvmConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration.ScalingType;
import qspr.modelling.OptionsEnumeration;
import qspr.modelling.configurations.CDSConfiguration;
import qspr.util.BasicRecordMapper;
import qspr.workflow.utils.QSPRConstants;

public class LibSVMConfigurator extends DescriptorsConfigurator 
{
	OptionsEnumeration en;

	@Override
	public void setModel(Model model)
	{
		super.setModel(model);
		LibSvmConfiguration conf = new LibSvmConfiguration();
		conf.scaleTypeX = ScalingType.RANGE;
		conf.scaleTypeY = ScalingType.NONE;
		((CDSConfiguration)model.attachment.getObject().configuration).modelConfiguration = conf;	
		this.stepAfterDescriptors = "libsvm";
	}

	public WebModel libsvm()
	{
		model.trainingSet = (Basket)Globals.session().get(Basket.class, model.trainingSet.id);
		en = OptionsEnumeration.enumeratePropertyOptions(model.getFilteredSet(QSPRConstants.TRAINING), new BasicRecordMapper(model.getFilteredSet(QSPRConstants.TRAINING)));
		List<PropertyOption> list = new ArrayList<PropertyOption>();
		for (Long option_id : en.optionsMapping.keySet())
			list.add((PropertyOption)Globals.session().get(PropertyOption.class, option_id));
		return new WebModel(getDefaultTemplate()).setList(list).setTemplate("modeller/configurators/libsvm");
	}

	public void libsvmSubmit(HttpServletRequest request)
	{
		LibSvmConfiguration conf = (LibSvmConfiguration) 
				((CDSConfiguration)model.attachment.getObject().configuration).modelConfiguration;
		//svm option
		Integer type = Integer.valueOf(request.getParameter("use-nu"));
		if (type == 0)
		{
			conf.useNu = false;
			conf.oneClass = false;
		} else if (type == 1)
		{
			conf.useNu = true;
			conf.oneClass = false;
		} else if (type == 2)
		{
			conf.useNu = true;
			conf.oneClass = true;
			Long option_id =  Long.valueOf(request.getParameter("one-class-class"));
			conf.oneClassLabel = en.optionsMapping.get(option_id).toString();
		}

		if (request.getParameter("weight") != null)
			conf.useWeighting = true;

		conf.kernel_type = Integer.parseInt(request.getParameter("kernel-type"));

		if (request.getParameter("advanced") != null)
		{
			//advanced/kernel options
			conf.cost = Double.parseDouble(request.getParameter("cost"));


			if (request.getParameter("gamma").equals(""))
				conf.gamma = null;
			else
				conf.gamma = Double.parseDouble(request.getParameter("gamma"));

			conf.svrEpsilon = Double.parseDouble(request.getParameter("svr-epsilon"));

			//advanced-advanced options
			conf.degree = Integer.parseInt(request.getParameter("degree"));
			conf.coef0 = Integer.parseInt(request.getParameter("coef0"));
			conf.nu = Double.parseDouble(request.getParameter("nu"));
			conf.epsilon = Double.parseDouble(request.getParameter("epsilon"));
			conf.classWeightRatio = Double.parseDouble(request.getParameter("class-ratio"));
		}


		//grid option
		//add grid option
		if(request.getParameter("grid") != null)
		{
			conf.gridSearchSetSize = Double.parseDouble(request.getParameter("grid-search-set-size"));
			if(request.getParameter("grid-search-parallel") !=null && request.getParameter("grid-search-parallel").length()>0)
				conf.PARALLEL = Integer.parseInt(request.getParameter("grid-search-parallel"));			

			conf.costMin = Double.parseDouble(request.getParameter("cost-min"));
			conf.costMax = Double.parseDouble(request.getParameter("cost-max"));
			conf.costStep = Double.parseDouble(request.getParameter("cost-step"));
			conf.gammaMin = Double.parseDouble(request.getParameter("gamma-min"));
			conf.gammaMax = Double.parseDouble(request.getParameter("gamma-max"));
			conf.gammaStep = Double.parseDouble(request.getParameter("gamma-step"));
			conf.svrEpsilonMin = Double.parseDouble(request.getParameter("svr-epsilon-min"));
			conf.svrEpsilonMax = Double.parseDouble(request.getParameter("svr-epsilon-max"));
			conf.svrEpsilonStep = Double.parseDouble(request.getParameter("svr-epsilon-step"));

			conf.classWeightRatioMin = Double.parseDouble(request.getParameter("class-ratio-min"));
			conf.classWeightRatioMax = Double.parseDouble(request.getParameter("class-ratio-max"));
			conf.classWeightRatioStep = Double.parseDouble(request.getParameter("class-ratio-step"));
		}
		else {
			conf.gridSearch = false;
			conf.PARALLEL = null;
			conf.gridSearchSetSize = conf.costMin = conf.costMax = conf.costStep = conf.gammaMin = conf.gammaMax = null;
			conf.gammaStep = conf.svrEpsilonMin = conf.svrEpsilonMax = conf.svrEpsilonStep = null;
			conf.classWeightRatioMin = conf.classWeightRatioMax = conf.classWeightRatioStep = null;
		}

		currentPage = startStep;
	}
}
