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

package qspr.controllers.login.providers.ldap;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import qspr.controllers.login.providers.ProviderConfig;

public class LDAPProviderConfig implements ProviderConfig {
	
	Map<String,LDAPServer> servers = new HashMap<>();
	
	public LDAPProviderConfig(List<LDAPServer> servers) {
		for (LDAPServer server : servers) {
			this.servers.put(server.name, server);
		}
	}

	@Override
	public String getURLName() {
		return "LDAP";
	}
	
	public Collection<LDAPServer> getServers() {
		return servers.values();
	}

	public LDAPServer getServer(String serverName) {
		return servers.get(serverName);
	}

}
