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

package qspr.controllers.login.providers;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlRootElement;

import java.util.Map;

import qspr.OCHEMConfiguration;
import qspr.controllers.login.providers.facebook.FacebookLoginHandler;
import qspr.controllers.login.providers.facebook.FacebookProviderConfig;
import qspr.controllers.login.providers.ldap.LDAPLoginHandler;
import qspr.controllers.login.providers.ochem.OCHEMLoginConfig;
import qspr.controllers.login.providers.ochem.OCHEMLoginHandler;

public class LoginHandlers {
	
	@XmlRootElement
	public static class HandlerInfo {
		
		@XmlAttribute
		String name;
		
		public HandlerInfo() {
			// no action
		}
		
		public HandlerInfo(String name) {
			this.name = name;
		}
	}
	
	static final Map<String,LoginByProviderHandler> handlerMap = new HashMap<String, LoginByProviderHandler>();
	static {	
		handlerMap.put(
				"Facebook", 
				new FacebookLoginHandler(
						new FacebookProviderConfig(
								"Facebook",
								OCHEMConfiguration.loginInfo.get("facebook.api_key"),
								OCHEMConfiguration.loginInfo.get("facebook.secret")
						)
				)
		);
		handlerMap.put(
				"OCHEM",
				new OCHEMLoginHandler(
						new OCHEMLoginConfig(
								"OCHEM",
								OCHEMConfiguration.loginInfo.get("ochem.url"),
								OCHEMConfiguration.loginInfo.get("ochem.api_key"),
								OCHEMConfiguration.loginInfo.get("ochem.secret")
						)
				)
		);
		handlerMap.put(
				"LDAP", 
				LDAPLoginHandler.createFromOCHEMConfig()
		);
	};
	
	public static LoginByProviderHandler get(String urlName) {
		return handlerMap.get(urlName);
	}
	
	public static List<HandlerInfo> getInfos() {
		List<HandlerInfo> ret = new ArrayList<>();
		for (String key : handlerMap.keySet()) {
			ret.add(new HandlerInfo(handlerMap.get(key).getConfig().getURLName()));
		}
		return ret;
	}
}
