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

import qspr.dao.Repository;
import qspr.entities.Attachment.AttachmentType;
import qspr.entities.Model;
import qspr.entities.ModelProtocol;
import qspr.frontend.WebModel;
import qspr.metaserver.configurations.ASNNConfiguration;
import qspr.metaserver.configurations.BaggingConfiguration;
import qspr.metaserver.configurations.SelectionConfiguration;
import qspr.modelling.configurations.CDSConfiguration;
import qspr.modelling.configurations.CDSModelData;

import com.eadmet.exceptions.UserFriendlyException;

public class LibraryModelConfigurator extends BasicModelConfigurator 
{

	@Override
	public void setModel(Model model)
	{
		super.setModel(model);
	}

	public LibraryModelConfigurator()
	{
		firstConfigurationStep = "chooseModel";
		dataPreprocessingStep = false;
	}

	public WebModel chooseModel()
	{
		return new WebModel(getDefaultTemplate()).setTemplate("modeller/configurators/choose-library-model");
	}

	public void chooseModelSubmit(HttpServletRequest request)
	{
		Model baseModel = Repository.model.getById(Long.valueOf(request.getParameter("basemodel-id")));

		if (baseModel.attachment.getObject().protocol.validationConfiguration instanceof BaggingConfiguration)
			throw new UserFriendlyException("The base model created using bagging are not yet supported. Please, choose another base model.");

		ModelProtocol chosenProtocol = model.attachment.getObject().protocol;

		// Copy the model attachment completely form the base model
		model.attachment.setObject(baseModel.attachment.getObject(), AttachmentType.MARSHALABLE); // Copy the model protocol completely

		// .. but keep the model protocol selected by the user
		model.attachment.getObject().protocol = chosenProtocol;

		// Copy the base model data
		ASNNConfiguration annConf = (ASNNConfiguration) ((CDSConfiguration)model.attachment.getObject().configuration).modelConfiguration;

		if(((CDSModelData)baseModel.readyModelAttachment.getObject().modelData).methodSpecificData==null)throw new UserFriendlyException(" This model was not saved and cannot be applied to new data");

		ASNNConfiguration annConfTrained = (ASNNConfiguration) ((CDSModelData)baseModel.readyModelAttachment.getObject().modelData).methodSpecificData;
		if(annConfTrained.savedmodel==null)throw new UserFriendlyException("Library model can be only used with new ASNN model which were built after 1.08.2012");

		annConf.libraryModel = annConfTrained.savedmodel;
		annConf.savedmodel = null;

		if (request.getParameter("property-num") != null && !"".equals(request.getParameter("property-num")))
			annConf.libraryOutput = Integer.valueOf(request.getParameter("property-num"));

		// Take the descriptors filters from the base model
		SelectionConfiguration descConfTrained = ((CDSModelData)baseModel.readyModelAttachment.getObject().modelData).selectionConfiguration;
		((CDSConfiguration)model.attachment.getObject().configuration).selection = descConfTrained;

		currentPage = "start";
	}
}
