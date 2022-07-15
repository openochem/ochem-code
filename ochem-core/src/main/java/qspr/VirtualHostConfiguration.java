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

package qspr;

import java.util.HashMap;
import java.util.Map;

import com.eadmet.utils.config.ConfigurableClass;
import com.eadmet.utils.config.ConfigurableProperty;


@ConfigurableClass(name = "vhosts", comment = "Configurations specific to particular access URLs ('virtual hosts')")
public class VirtualHostConfiguration
{
	/**
	 * This map allows to define different themes for different access URLs.
	 * A typical use case - OCHEM trade is physically the same server as iPRIOR, but is accessed via ochemtrade.com and has an own theme 
	 */
	@ConfigurableProperty(name = "outer_template", comment = "Allows to bind domain names to particular 'themes'")
	public static Map<String, String> outerTemplate = new HashMap<String,String>();
	
	/**
	 * This map allows to define different themes for different access URLs.
	 * A typical use case - OCHEM trade is physically the same server as iPRIOR, but is accessed via ochemtrade.com and has an own theme 
	 */
	@ConfigurableProperty(name = "home_page", comment = "Allows to bind domain names to particular 'themes'")
	public static Map<String, String> defaultHomePage = new HashMap<String,String>();
	
	
	public static String getOuterTemplate()
	{
		return outerTemplate.get(getHostName());
	}
	
	public static String getHomePage()
	{
		return defaultHomePage.get(getHostName());
	}
	
	private static  String getHostName()
	{
		return OCHEMConfiguration.getRootHost().substring(7).replace(":8080", "");
	}
	

}
