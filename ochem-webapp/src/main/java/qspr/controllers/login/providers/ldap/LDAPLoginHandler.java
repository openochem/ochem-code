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

import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.naming.AuthenticationException;
import javax.naming.AuthenticationNotSupportedException;
import javax.naming.NamingException;
import javax.naming.directory.SearchResult;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.json.JSONObject;
import org.springframework.web.servlet.ModelAndView;

import com.eadmet.exceptions.UserFriendlyException;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.controllers.LoginController;
import qspr.controllers.login.providers.LoginByProviderHandler;
import qspr.controllers.login.providers.ProviderConfig;
import qspr.entities.User;
import qspr.util.DynaWrap;

public class LDAPLoginHandler extends LoginByProviderHandler {

	public LDAPLoginHandler(ProviderConfig config) {
		super(config);
	}

	@Override
	public String getProviderLoginURL() {
		return OCHEMConfiguration.getRootHost() + "/login/ldap.do?render-mode=popup";
	}

	@Override
	protected ProviderCallbackData parseCallbackRequest(HttpServletRequest req) {
		// TODO Auto-generated method stub
		return null;
	}
	
	@Override
	protected String getAuthURL(ProviderCallbackData callbackData) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	protected ProviderAccessToken parseAuthResponse(JSONObject data) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	protected User getUser(ProviderAccessToken token) {
		// TODO Auto-generated method stub
		return null;
	}
	
	@Override
	protected User createUser(ProviderAccessToken token) throws UserFriendlyException {
		// TODO Auto-generated method stub
		return null;
	}
	
	@Override
	public ModelAndView callbackController(HttpServletRequest req, HttpServletResponse res, LoginController caller) {
		String username = req.getParameter("username");
		String password = req.getParameter("password");
		String serverName = req.getParameter("ldapServer");
		
		LDAPProviderConfig config = (LDAPProviderConfig) getConfig();
		LDAPServer server = config.getServer(serverName);
		SearchResult result;
		try {
			result = server.authenticate(username, password);
		} catch (AuthenticationNotSupportedException e) {
			e.printStackTrace();
			return caller.redirect("login/ldap.do?render-mode=popup&error=unsupported");
		} catch (AuthenticationException e) {
			return caller.redirect("login/ldap.do?render-mode=popup&error=invalid_credentials");
		} catch (NamingException e) {
			return caller.redirect("login/ldap.do?render-mode=popup&error=unknown");
		}
		
		User user = getUser(result);
		if (user != null) {
			caller.loginUser(user);
			logUserLogin(user, req);
		} else {
			user = createUser(result);
			Globals.session().save(user);
			caller.loginUser(user);
			logUserCreate(user, req);
			logUserLogin(user, req);
		}
		
		return caller.redirect("login/chooseProvider.do?login=" + user.login);
	}
	
	protected User getUser(SearchResult result) {
		return User.getByLogin(encodeUsername(result.getNameInNamespace()));
	}
	
	protected User createUser(SearchResult result) throws UserFriendlyException {
		String nsName = result.getNameInNamespace();
		String login = encodeUsername(nsName);
		User sameLoginUser = User.getByLogin(login);
		if (sameLoginUser != null)
			throw new UserFriendlyException("Cannot create a new account from profile, since userID " + sameLoginUser.login + " (LDAP namespace name: "+ nsName + ") already exists. If its you, please login or request a forgotten password.");
		User user = User.getNewInstance();
		user.login = login;
		user.rank = User.NOTVALIDATED;

		String cn = "";
		String username = "";
		try {
			username = result.getAttributes().get("uid").get(0).toString();
			cn = result.getAttributes().get("cn").get(0).toString();
		} catch (NamingException e) {
			throw new UserFriendlyException(e);
		}
		String[] names = cn.split("\\s+");
		if (names.length > 1 && (user.isExtended())) {
			DynaWrap extended =  user.dynaWrapped();
			extended.setField("firstName", names[0]);
			String last = names[1];
			if (names.length > 2) {
				for (int i = 2; i < names.length; i++) {
					last += " " + names[i];
				}
			}
			extended.setField("lastName", last);
			try {
				User.getExtended().getMethod("setPassword", int.class).invoke(user, 8);
			} catch (Exception e) {
				throw new UserFriendlyException(e); 
			}
		} else {
			user.login = login;
		}
		
		return user;
	}

	private static String bytesToHex(byte[] hash) {
	    StringBuilder hexString = new StringBuilder(2 * hash.length);
	    for (int i = 0; i < hash.length; i++) {
	        String hex = Integer.toHexString(0xff & hash[i]);
	        if(hex.length() == 1) {
	            hexString.append('0');
	        }
	        hexString.append(hex);
	    }
	    return hexString.toString();
	}
	
	private static String encodeUsername(String username) {
		MessageDigest digest;
		try {
			digest = MessageDigest.getInstance("SHA-256");
		} catch (NoSuchAlgorithmException e) {
			throw new UserFriendlyException(e);
		}
		byte[] encodedhash = digest.digest(
				username.getBytes(StandardCharsets.UTF_8));
		return bytesToHex(encodedhash).substring(0, 16);
	}
	
	public static LDAPLoginHandler createFromOCHEMConfig() {
		Map<String, Map<String, String>> servers_config = new HashMap<>();
		for (Map.Entry<String,String> entry : OCHEMConfiguration.loginInfo.entrySet()) {
			String key = entry.getKey();
			if (key.startsWith("ldap")) {
				String[] config  = key.split("\\.");
				if (config.length < 3) {
					throw new UserFriendlyException("Incorrect LDAP server configuration: " + key);
				}
				String server_name = config[1];
				if (!servers_config.containsKey(server_name)) {
					servers_config.put(config[1], new HashMap<String,String>());
				}
				String config_prop = config[2];
				servers_config.get(server_name).put(config_prop, entry.getValue());
			}
		}
		List<LDAPServer> servers = new ArrayList<>();
		for (Map.Entry<String,Map<String, String>> entry : servers_config.entrySet()) {
			servers.add(LDAPServer.fromPropertyMap(entry.getKey(), entry.getValue()));
		}
		
		return new LDAPLoginHandler(new LDAPProviderConfig(servers));
	}
}
