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

package com.eadmet.useractions;

import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement(name = "login-action")
public class LoginAction extends AbstractUserAction
{
	public String userAgent;
	public String ipAddress;
	public String provider;
	
	public LoginAction()
	{
		
	}
	
	public LoginAction(String userAgent, String ipAddress, String provider)
	{
		this.userAgent = userAgent;
		this.ipAddress = ipAddress;
		this.provider = provider;
	}
	
	@Override
	public String getLogLine()
	{
		return "has logged in " + (provider != null ? " with his " + provider + " account" : "") + " using " + getBrowserName(userAgent);
	}
	
	private static String getBrowserName(String userAgent) {
		
		if (userAgent == null)
			return "Unknown browser";
		
	    if(userAgent.contains("MSIE")){ 
	        return "Internet Explorer";
	    }
	    if(userAgent.contains("Firefox")){ 
	        return "Firefox";
	    }
	    if(userAgent.contains("Chrome")){ 
	        return "Chrome";
	    }
	    if(userAgent.contains("Opera")){ 
	        return "Opera";
	    }
	    if(userAgent.contains("Safari")){ 
	        return "Safari";
	    }
	    return userAgent;
	}
	
}
