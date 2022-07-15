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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.servlet.http.HttpServletRequest;
import javax.xml.bind.annotation.XmlRootElement;

import qspr.interfaces.Descriptable;

@XmlRootElement(name = "descriptors-configuration")
public class DescriptorsConfiguration  extends DescriptorsAbstractConfiguration implements Descriptable
{
	private static final long serialVersionUID = 4L;

	public static final String RDKIT = "RDKIT";

	public static final String CDK2 = "CDK2";
	public static final String JPLOGP = "JPlogP";
	public static final String CDDD = "CDDD";
	public static final String MORDRED = "MORDRED";
	public static final String SilicosItScaffold = "SilicosItScaffold";
	public static final String ExternalDescriptors = "ExternalDescriptors";
	public static final String FRAGMENTS = "Fragmentor";
	public static final String GSFrag = "GSFrag";
	public static final String INDUCTIVE = "InductiveDescriptors";
	public static final String Mera = "Mera";
	public static final String Mersy = "Mersy";	
	public static final String MolPrint = "MolPrint";
	public static final String MOLD2 = "MOLD2";
	public static final String MOPAC2016 = "MOPAC2016";
	public static final String KRAKENX = "KrakenX";
	public static final String SIGMA = "SIGMA";
	public static final String MOPAC = "MOPAC";
	public static final String OEstate = "OEstate";
	public static final String PyDescriptor = "PyDescriptor";
	public static final String QNPR = "QNPR";
	public static final String SIRMS = "SIRMS";
	public static final String MAP4 = "MAP4";
	public static final String Spectrophores = "Spectrophores";
	public static final String StructuralAlerts = "StructuralAlerts";
	public static final String PADEL2 = "PaDEL2";
	public static final String ECFP = "ECFP";

	//Special
	public static final String ApplyModel = "ApplyModel"; // a special type, starts a model via web-service
	public static final String ALogPS = "ALogPS";
	public static final String RANDOM = "Random";
	public static final String ExpValues = "ExpValues";
	public static final String EPA = "EPA";


	public List<DescriptorType> types = new ArrayList<DescriptorType>();

	public MixturesProcessing mixtures;

	public Boolean allowMerge;

	public Boolean forceUpdateDescriptorCache;


	/** 
	 * Indicates descriptors in the order as they will be processed 
	 * All processed descriptors should be specified in this file
	 */

	final public static String descriptors[] = {
			ALogPS,
			ApplyModel,
			CDDD,
			CDK2,
			EPA,
			ExternalDescriptors,
			FRAGMENTS,
			GSFrag,
			INDUCTIVE,
			JPLOGP,
			KRAKENX,
			MAP4,
			Mera,
			Mersy,
			MOLD2,
			MolPrint,
			MOPAC,
			MOPAC2016,
			MORDRED,
			OEstate,
			ECFP,
			PADEL2,
			PyDescriptor,
			QNPR,
			RDKIT,
			SIGMA,
			SilicosItScaffold,
			SIRMS,
			Spectrophores,
			StructuralAlerts,

	};



	public DescriptorsConfiguration() {
	}

	public String getInformativeName() {
		return types.size()>0?types.toString():"";
	}

	@Override
	public DescriptorsAbstractConfiguration setConfiguration(HttpServletRequest request) throws Exception {
		types.clear();

		for(String descriptor: descriptors)
			if (request.getParameterValues(descriptor) != null)
				addByName(descriptor, request);
		return this;
	}

	void addByName(String name, HttpServletRequest request) throws Exception {
		switch(name) {

		case RANDOM:
			addDescriptorType(RANDOM,null); break;

		case ALogPS:
			addDescriptorType(name,new DescriptorsEmptyConfiguration()); break;

		default: 
			@SuppressWarnings("rawtypes")
			Class myClass = Class.forName("qspr.metaserver.configurations.Descriptors"+name+"Configuration");
			DescriptorsAbstractConfiguration o = (DescriptorsAbstractConfiguration) myClass.newInstance();
			o = o.setConfiguration(request);
			if(o != null)
				addDescriptorType(name,o);
		}
	}

	public DescriptorsConfiguration(DescriptorsConfiguration old){
		types = old.types;
		mixtures = old.mixtures;
		allowMerge = old.allowMerge;
		forceUpdateDescriptorCache = old.forceUpdateDescriptorCache;
		types = old.types;
	}

	public DescriptorType addDescriptorType(String taskName, DescriptorsAbstractConfiguration configuration)
	{
		return addDescriptorType(taskName, configuration, null, null, false);
	}

	public DescriptorType addDescriptorType(String taskName, DescriptorsAbstractConfiguration configuration, String columnTitle, Boolean writeToCache, boolean atFront)
	{
		DescriptorType type = new DescriptorType();
		type.type = taskName;
		type.configuration = configuration;
		type.title = columnTitle;
		type.writeToCache = writeToCache;

		if (types.contains(type)) { // replaces at the same position of the list
			int index = types.indexOf(type);
			types.remove(type);
			types.add(index, type);
		}
		else
			if(atFront)
				types.add(0, type);
			else
				types.add(type);

		return type;
	}

	public String toString()
	{
		String type; 
		if (types.isEmpty())
			type = "";
		else
			type = types.toString();

		if (mixtures != null)
		{
			type += "\n Mixture processing : ";
			switch (mixtures) 
			{
			case FRACTION: type += "Use molar fractions as descriptors";break;
			case CONCATENATE: type += "Concatenate descriptors";break;
			case WSUMDIFF: type += "Sum and absolute differences of weighted (by molar fraction) descriptors";break;
			case WSUM: type += "Descriptor sum (weighted by molar fraction)";break;
			default:
				type += "Unknown";
				break;
			}

		}

		return type;
	}

	public boolean requires3D()
	{

		for (DescriptorType descType : types)
			if (descType.requires3D())
				return true;
		return false;
	}

	public enum MixturesProcessing
	{
		NONE, FRACTION, CONCATENATE, WSUM, WSUMDIFF, SUMONLY; ///NONE will correspond to null
		;

		public static MixturesProcessing fromString(String arg)
		{
			if(arg == null) return null;

			switch (arg) {
			case "sumonly":
				return MixturesProcessing.SUMONLY;
			case "fraction":
				return MixturesProcessing.FRACTION;
			case "concatenate":
				return MixturesProcessing.CONCATENATE;
			case "wsumdiff":
				return MixturesProcessing.WSUMDIFF;
			case "wsum":
				return MixturesProcessing.WSUM;
			default:
				return null;
			}
		}
	}

	public Map<String, Object> getParameters() 
	{
		Map<String, Object> parameters = new HashMap<String, Object>();

		List<String> dTypes = new ArrayList<String>();
		for (DescriptorType type : types) 
			dTypes.add(type.type);

		parameters.put("Descriptors", dTypes.toString());

		return parameters;
	}

	public boolean allowsMissedValues() {
		return allowMerge == null ? false :  allowMerge;
	}

	public boolean supportSplicing() {
		for(DescriptorType type: types)
			if(type.supportSplicing())return true;
		return false;
	}

	@Override
	public String getDefaultTypeName() {
		return "DESCRIPTORS";
	}

}
