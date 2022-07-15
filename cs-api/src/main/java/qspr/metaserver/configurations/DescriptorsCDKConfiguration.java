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
import java.util.List;

import javax.servlet.http.HttpServletRequest;

import qspr.dao.Aromatisation;

abstract public class DescriptorsCDKConfiguration extends DescriptorsAbstractConfiguration {
	private static final long serialVersionUID = 1L;

	static final public String constitutionalDescriptor = "constitutionalDescriptor";
	static final public String topologicalDescriptor = "topologicalDescriptor";
	static final public String geometricalDescriptor = "geometricalDescriptor";
	static final public String electronicDescriptor = "electronicDescriptor";
	static final public String hybridDescriptor = "hybridDescriptor";

	public List<String> descriptorTypes = new ArrayList<String>();

	public DescriptorsCDKConfiguration(){
		addAll();
	}

	private void addAll(){
		descriptorTypes = new ArrayList<String>();
		descriptorTypes.add(constitutionalDescriptor);
		descriptorTypes.add(topologicalDescriptor);
		descriptorTypes.add(geometricalDescriptor);
		descriptorTypes.add(electronicDescriptor);
		descriptorTypes.add(hybridDescriptor);
	}

	@Override
	public boolean requires3D() {
		return  true;
	}

	public DescriptorsCDKConfiguration(String[] classesList){
		this.descriptorTypes.addAll(Arrays.asList(classesList));
	}

	@Override
	public String toString() {
		String str = "";
		for(String s: descriptorTypes)
			switch(s) {
			case constitutionalDescriptor: str += "cons,"; break;
			case topologicalDescriptor: str += "topol,"; break;
			case geometricalDescriptor: str += "geom,"; break;
			case electronicDescriptor: str += "elect,"; break;
			default: str += "hybr";
			}

		return str.replaceAll(",$", "");
	}

	@Override
	public boolean isLongCalculation() {
		return true;
	}

	@Override
	public DescriptorsAbstractConfiguration setAllOn() {
		addAll();
		return this;
	}

	public Aromatisation getCDKAromatization() {
		return null;
	}

	@Override
	public void cleanXML(){
		super.cleanXML();
	}

	@Override
	public boolean isCachable() {
		return true;
	}


	@Override
	DescriptorsCDKConfiguration setConfiguration(HttpServletRequest request) {
		descriptorTypes = new ArrayList<String>();
		if (request.getParameter("cdk1") != null)
			descriptorTypes.add(constitutionalDescriptor);
		if (request.getParameter("cdk2") != null)
			descriptorTypes.add(topologicalDescriptor);
		if (request.getParameter("cdk3") != null)
			descriptorTypes.add(geometricalDescriptor);
		if (request.getParameter("cdk4") != null)
			descriptorTypes.add(electronicDescriptor);
		if (request.getParameter("cdk5") != null)
			descriptorTypes.add(hybridDescriptor);

		if (request.getParameter("cdkTimeout") != null)
			setTimeout(Integer.valueOf(request.getParameter("cdkTimeout")));

		return this;	
	}
}
