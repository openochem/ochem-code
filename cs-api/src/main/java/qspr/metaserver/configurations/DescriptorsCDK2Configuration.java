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

package qspr.metaserver.configurations;

import java.util.ArrayList;
import java.util.Arrays;

import javax.servlet.http.HttpServletRequest;
import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement(name = "cdk2-configuration")
public class DescriptorsCDK2Configuration extends DescriptorsCDKConfiguration {
	private static final long serialVersionUID = 1L;

	public DescriptorsCDK2Configuration(){
	}

	public DescriptorsCDK2Configuration(String[] classesList){
		this.descriptorTypes.addAll(Arrays.asList(classesList));
	}

	@Override
	public String getDefaultTypeName() {
		return DescriptorsConfiguration.CDK2;
	}

	@Override
	DescriptorsCDK2Configuration setConfiguration(HttpServletRequest request) {
		descriptorTypes = new ArrayList<String>();
		if (request.getParameter("cdk21") != null)
			descriptorTypes.add(constitutionalDescriptor);
		if (request.getParameter("cdk22") != null)
			descriptorTypes.add(topologicalDescriptor);
		if (request.getParameter("cdk23") != null)
			descriptorTypes.add(geometricalDescriptor);
		if (request.getParameter("cdk24") != null)
			descriptorTypes.add(electronicDescriptor);
		if (request.getParameter("cdk25") != null)
			descriptorTypes.add(hybridDescriptor);

		if (request.getParameter("cdk2Timeout") != null)
			setTimeout(Integer.valueOf(request.getParameter("cdk2Timeout")));
		return this;	
	}
}
