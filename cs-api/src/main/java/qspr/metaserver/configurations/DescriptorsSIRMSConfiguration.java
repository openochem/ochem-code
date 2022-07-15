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
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.servlet.http.HttpServletRequest;
import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement(name = "sirms-configuration")
public class DescriptorsSIRMSConfiguration extends DescriptorsAbstractConfiguration{
	private static final long serialVersionUID = 1L;
	/**
	 * Do not use H (even explicit)
	 */
	public Boolean noH = true;
	/**
	 * Use mixtures for each compound (including quasi mixtures)
	 */
	public Boolean	mixtures = null;

	public Integer minFragment = 1;
	public Integer maxFragment = 4;

	public List<Labeling> descriptorTypes = new ArrayList<Labeling>();

	public enum Labeling {
		none,
		elm,
		CHARGE,
		LOGP,
		HB,
		REFRACTIVITY
	};

	public DescriptorsSIRMSConfiguration(int min, int max){
		minFragment = min;
		maxFragment = max;
		addDefault();
	}

	public DescriptorsSIRMSConfiguration(){
		addDefault();
	}

	void addDefault(){
		descriptorTypes.add(Labeling.CHARGE);
		descriptorTypes.add(Labeling.LOGP);
		descriptorTypes.add(Labeling.HB);
		descriptorTypes.add(Labeling.REFRACTIVITY);
	}

	public boolean ignoreH(){
		return noH == null || noH;
	}

	public boolean asMixtures(){
		return mixtures != null && mixtures;
	}

	public String toString()
	{
		String lab = "";
		for(Labeling l: getUniqueAndSorted())
			lab += l + " ";

		lab = lab.trim().replace(" ", "+");

		return "" + ("labels: "+lab.toLowerCase()); 
	}

	@Override
	public boolean requires3D() {
		return false;
	}

	@Override
	public boolean isCachable() {
		return true;
	}

	@Override
	public String getDefaultTypeName() {
		return DescriptorsConfiguration.SIRMS;
	}

	public List<Labeling> getUniqueAndSorted(){
		Set<Labeling> unique = new HashSet<Labeling>(descriptorTypes);
		descriptorTypes.clear();
		unique.remove(Labeling.LOGP); // not supported by CDK
		unique.remove(Labeling.REFRACTIVITY);  // not supported by CDK
		descriptorTypes.addAll(unique);
		Collections.sort(descriptorTypes);
		return descriptorTypes;
	}

	@Override
	public void cleanXML(){
		super.cleanXML();
		getUniqueAndSorted();
	}

	@Override
	DescriptorsAbstractConfiguration setConfiguration(HttpServletRequest request) {
		descriptorTypes = new ArrayList<Labeling>();

		if(request.getParameter("sirms_noh")!=null)
			noH = true;
		if(request.getParameter("sirms_mixtures")!=null)
			mixtures = true;
		minFragment = Integer.valueOf(request.getParameter("sirms-min-frag"));
		maxFragment = Integer.valueOf(request.getParameter("sirms-max-frag"));

		if(request.getParameter("sirms_none")!=null)
			descriptorTypes.add(Labeling.none);
		if(request.getParameter("sirms_elm")!=null)
			descriptorTypes.add(Labeling.elm);
		if(request.getParameter("sirms_charge")!=null)
			descriptorTypes.add(Labeling.CHARGE);
		if(request.getParameter("sirms_logp")!=null)
			descriptorTypes.add(Labeling.LOGP);
		if(request.getParameter("sirms_hb")!=null)
			descriptorTypes.add(Labeling.HB);
		if(request.getParameter("sirms_refractivity")!=null)
			descriptorTypes.add(Labeling.REFRACTIVITY);

		return this;
	}

}
