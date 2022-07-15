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

package qspr.business;

import qspr.Globals;
import qspr.entities.Model;
import qspr.metaserver.configurations.BaggingConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.configurations.NoDescriptors;
import qspr.modelling.configurations.CDSConfiguration;

public class ModelPeer 
{

	public static final int NAME_LENGTH = 30;

	public static String getDescriptorsCode(ModelAbstractConfiguration conf)
	{
		return (conf == null || !(conf  instanceof NoDescriptors))?"": ((NoDescriptors)conf).augemenationString();
	}

	public static String getModelName(Model model, boolean randomize)
	{
		if (model.template.isDescriptorCalculationOnly())
			return ((CDSConfiguration) model.attachment.getObject().configuration).descriptors.toString();
		StringBuilder name = new StringBuilder();

		for (int i = 0; i < model.modelMappings.size(); i++)
			name.append((i == 0 ? "" : "+") + model.modelMappings.get(i).property.getName());

		if(name.length() > NAME_LENGTH)
			name.delete(NAME_LENGTH, name.length());

		if (model.attachment.getObject().configuration instanceof CDSConfiguration) {
			CDSConfiguration cds = (CDSConfiguration) model.attachment.getObject().configuration;
			String descName = cds.descriptors.getInformativeName(); descName = descName.length() >0? "_" + descName:"";
			String augmentation = getDescriptorsCode(cds.modelConfiguration); augmentation = augmentation.length() >0? "_" + augmentation:"";
			name.append( "_" + cds.modelConfiguration.getInformativeName() + descName + augmentation);
		}

		if (model.attachment.getObject().protocol.validationConfiguration instanceof BaggingConfiguration)
			name.append("_Bag");

		String result = (name.length() > (NAME_LENGTH * 3)) ? name.substring(0, NAME_LENGTH * 2) + "..." : name.toString();

		if (randomize)
		{
			Long id = (Long) Globals.session().createQuery("select max(id) from Model").uniqueResult() + 1;
			result = result + (" - " + id);
		}

		return result;
	}
}
