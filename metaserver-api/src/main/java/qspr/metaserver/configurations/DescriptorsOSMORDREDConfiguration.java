package qspr.metaserver.configurations;

import javax.servlet.http.HttpServletRequest;
import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement(name="osmo-configuration")
public class DescriptorsOSMORDREDConfiguration extends DescriptorsAbstractConfiguration
{
	private static final long serialVersionUID = 1L;

	public DescriptorsOSMORDREDConfiguration()
	{

	}

	public String toString()
	{
		return ""; 
	}

	@Override
	public boolean requires3D() {
		return false;
	}

	@Override
	public String getDefaultTypeName() {
		return DescriptorsConfiguration.OSMORDRED;
	}

	@Override
	DescriptorsOSMORDREDConfiguration setConfiguration(HttpServletRequest request) {
		return this;
	}
}
