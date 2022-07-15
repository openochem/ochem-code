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

package qspr.modelling.configurations;

import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.metaserver.configurations.LinearConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.configurations.NoDescriptors;
import qspr.metaserver.configurations.SelectionConfiguration;
import qspr.metaserver.configurations.StructureOptimisationConfiguration;


@XmlRootElement(name = "cds-model-data")
public class CDSModelData 
{
	//	public boolean useCorina;
	public StructureOptimisationConfiguration optimisationConfiguration;
	public DescriptorsConfiguration descriptors;
	public SelectionConfiguration selectionConfiguration;


	/**
	 * These data are null for model that was not saved
	 */

	@XmlElement(name = "method-data")
	public ModelAbstractConfiguration methodSpecificData;

	public String toString()
	{
		String res = 
				selectionConfiguration == null || methodSpecificData instanceof NoDescriptors? "no descriptors" : 
					(selectionConfiguration.getDescriptorsSize() == 0)?"descriptors list is not avaialble - model will not work"
							:(selectionConfiguration.getDescriptorsSize() + " pre-filtered descriptors");
		if (methodSpecificData != null)
		{

			res += "\n" + methodSpecificData;

			String equation = getFullEquation();
			if (equation != null)
				res += "\n" + equation;
		}else
			res += "\nModel was not saved.";

		return res;
	}

	@XmlTransient
	public String getFullEquation()
	{
		if (methodSpecificData!=null && methodSpecificData instanceof LinearConfiguration)
		{
			LinearConfiguration mlraConf = (LinearConfiguration) methodSpecificData;			
			String equation = mlraConf.theEquation(selectionConfiguration.descriptorAsStrings(),false);
			if(equation==null)return "\nModel was not saved.";			
			equation += "\n\nNormalised "+mlraConf.theEquation(selectionConfiguration.descriptorAsStrings(),true);
			return equation.toString();
		}

		return null;
	}


}
