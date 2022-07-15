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

package qspr.util;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.JAXBException;

import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.ThreadScope;
import qspr.dao.Repository;
import qspr.entities.ModelConfigurationTemplate;
import qspr.entities.Session;
import qspr.export.ExportableModel;
import qspr.export.ModelTemplatesList;
import qspr.modelling.configurators.ModelConfigurationTemplateFactory;


public class ModelTemplateTranser 
{
	@SuppressWarnings("unchecked")
	public void exportTemplates() throws FileNotFoundException, JAXBException
	{
		List<ModelConfigurationTemplate> mcTemplates = Globals.session()
				.createCriteria(ModelConfigurationTemplate.class)
				.add(Restrictions.eq("isPublic", true))
				.list();
		List<ExportableModel> eModels = new ArrayList<ExportableModel>();

		for (ModelConfigurationTemplate mcTemplate : mcTemplates) {
			try {
				eModels.add((ExportableModel) mcTemplate.configuration.getObject());
			} catch (Exception e) {
				System.out.println("Warning: could not export template " + mcTemplate.name + ": " + e.getMessage());
			}
		}
		ModelTemplatesList list = new ModelTemplatesList();
		list.mcTemplates = mcTemplates;
		list.eModels = eModels;

		Globals.jaxbContext.createMarshaller().marshal(list, new FileOutputStream("/Users/midnighter/templates.xml"));
	}

	public void importTemplates() throws JAXBException
	{
		ModelTemplatesList templates = (ModelTemplatesList) Globals.jaxbContext.createUnmarshaller().unmarshal(new File("/Users/midnighter/templates.xml"));
		for (int i = 0; i < templates.mcTemplates.size(); i++)
		{
			ExportableModel eModel = templates.eModels.get(i);
			ModelConfigurationTemplate mcTemplateOrig = templates.mcTemplates.get(i);

			ModelConfigurationTemplate mcTemplate = ModelConfigurationTemplateFactory.create(eModel, mcTemplateOrig.type);
			mcTemplate.modelTemplate = Repository.modelTemplate.getByName(mcTemplateOrig.modelTemplate.name);

			if (mcTemplate.modelTemplate == null)
				System.out.println("WARNING: No method with name " + mcTemplateOrig.modelTemplate.name);

			mcTemplate.isPublic = true;
			mcTemplate.name = mcTemplateOrig.name;
			Globals.session().save(mcTemplate);
		}
	}

	public static void main(String[] args) throws JAXBException 
	{
		Globals.jaxbContext = Globals.createJAXBContext();
		Globals.startAllTransactions();
		try
		{
			ThreadScope.get().userSession = Session.getSuperuserSession();
			new ModelTemplateTranser().importTemplates();
			Globals.commitAllTransactions();
		}
		catch (Exception e)
		{
			e.printStackTrace();
			Globals.rollbackAllTransactions();
		}
		finally
		{
			ThreadScope.reset();
		}
	}
}
