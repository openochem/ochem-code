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

import java.util.Hashtable;
import java.util.Map;
import java.util.Properties;

import javax.naming.AuthenticationException;
import javax.naming.AuthenticationNotSupportedException;
import javax.naming.Context;
import javax.naming.NamingEnumeration;
import javax.naming.NamingException;
import javax.naming.directory.DirContext;
import javax.naming.directory.InitialDirContext;
import javax.naming.directory.SearchControls;
import javax.naming.directory.SearchResult;
import javax.xml.bind.annotation.XmlRootElement;

import com.eadmet.exceptions.UserFriendlyException;

@XmlRootElement
public class LDAPServer {
	public String name;
	private String url;
	private Hashtable<String, String> env;
	String rootDomain = "";
	
	public LDAPServer() {
		name = "default";
	}
	
	public LDAPServer(String name, Hashtable<String, String> env) {
		this.name = name;
		if (env.get(Context.PROVIDER_URL) == null) {
			throw new UserFriendlyException("LDAP server context must specify PROVIDER_URL. None supplied.");
		} else if (env.get(Context.SECURITY_PRINCIPAL) == null || env.get(Context.SECURITY_CREDENTIALS) == null) {
			throw new UserFriendlyException("LDAP server context must SECURITY_PRINCIPAL. None supplied.");
		}
		this.url = env.get(Context.PROVIDER_URL);
		
		if (env.get(Context.INITIAL_CONTEXT_FACTORY) == null) {
			env.put(Context.INITIAL_CONTEXT_FACTORY, "com.sun.jndi.ldap.LdapCtxFactory");
		}
		if (env.get(Context.SECURITY_AUTHENTICATION) == null) {
			env.put(Context.SECURITY_AUTHENTICATION, "simple");
		}
		this.env = env;
	}
	
	public LDAPServer(String name, Hashtable<String, String> env, String rootDomain) {
		this(name, env);
		this.rootDomain = rootDomain;
	}
	
	public static LDAPServer fromPropertyMap(String name, Map<String, String> props) {
		Hashtable<String, String> env = new Hashtable<String, String>();
		env.put(Context.PROVIDER_URL, props.get("url"));
		env.put(Context.SECURITY_PRINCIPAL, props.get("security_principal"));
		env.put(Context.SECURITY_CREDENTIALS, props.get("security_credentials"));
		if (props.containsKey("root_domain")) {
			return new LDAPServer(name, env, props.get("root_domain"));
		} else {

			return new LDAPServer(name, env);
		}
	}
	
	public String getURL() {
		return this.url;
	}
	
	public SearchResult searchUser(String username, String[] retAttributes) throws AuthenticationException, AuthenticationNotSupportedException, NamingException {
		SearchControls ctrls = new SearchControls();
	    ctrls.setReturningAttributes(retAttributes);
	    ctrls.setSearchScope(SearchControls.SUBTREE_SCOPE);
	    DirContext ctx = new InitialDirContext(env);
    	String[] parts = username.split("@");
    	String[] domainParts;
    	if (parts.length == 1) {
    		domainParts = rootDomain.split("\\.");
    	} else if (parts.length == 2) {
    		domainParts = parts[1].split("\\.");
    	} else {
    		throw new AuthenticationException("Authentication for user " + username + "failed: Malformed username.");
    	}
    	String searchCtx = "dc=" + String.join(",dc=", domainParts);
	    NamingEnumeration<SearchResult> answers = ctx.search(searchCtx, "(uid=" + parts[0]+ ")", ctrls);
    	if (answers.hasMoreElements()) {
    		SearchResult result = answers.nextElement();
    		ctx.close();
    		return result;
    	} else {
    		ctx.close();
    		throw new AuthenticationException("Authentication for user " + username + "failed: User does not exist.");	    
    	}
	}
	
	public SearchResult searchUser(String username) throws AuthenticationException, AuthenticationNotSupportedException, NamingException {
		return searchUser(username, new String[] { "uid", "cn"});
	}
	
	public SearchResult authenticate(String username, String password) throws AuthenticationException, AuthenticationNotSupportedException, NamingException {
		// search the user in the current context
		SearchResult searchResult = searchUser(username);
		String user = searchResult.getNameInNamespace();
	    Properties props = new Properties();
	    String url = env.get(Context.PROVIDER_URL);
        props.put(Context.INITIAL_CONTEXT_FACTORY, env.get(Context.INITIAL_CONTEXT_FACTORY));
        props.put(Context.PROVIDER_URL, url);
        props.put(Context.SECURITY_PRINCIPAL, user);
        props.put(Context.SECURITY_CREDENTIALS, password);
	    
        DirContext userCtx = new InitialDirContext(props);
	    System.out.println(user + " was authenticated with LDAP server: " + url);
	    userCtx.close();
	    return searchResult;
	}
}
