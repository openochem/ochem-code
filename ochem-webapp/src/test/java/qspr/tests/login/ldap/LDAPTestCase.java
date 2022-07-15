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

package qspr.tests.login.ldap;

import static org.junit.Assert.*;

import java.util.Hashtable;
import java.util.Properties;

import javax.naming.AuthenticationException;
import javax.naming.AuthenticationNotSupportedException;
import javax.naming.Context;
import javax.naming.NamingEnumeration;
import javax.naming.NamingException;
import javax.naming.directory.DirContext;
import javax.naming.directory.InitialDirContext;
import javax.naming.directory.SearchControls;

import org.junit.Test;

public class LDAPTestCase {

	@Test
	public void testServer() throws Exception {
//		String username = "joe";
//		String password = "joespassword";
//		
//		String url = "ldap://localhost:389";
//		Hashtable<String, String> env = new Hashtable<String, String>();
//		env.put(Context.INITIAL_CONTEXT_FACTORY, "com.sun.jndi.ldap.LdapCtxFactory");
//		env.put(Context.PROVIDER_URL, url);
//		env.put(Context.SECURITY_AUTHENTICATION, "simple");
//		env.put(Context.SECURITY_PRINCIPAL, "cn=admin,dc=ochem,dc=eu");
//		env.put(Context.SECURITY_CREDENTIALS, "demo");
//		
//		try {
//			DirContext ctx = new InitialDirContext(env);
//			System.out.println("connected");
////		    System.out.println(ctx.getEnvironment());
//		    
//		    SearchControls ctrls = new SearchControls();
//		    ctrls.setReturningAttributes(new String[] { "cn", "sn", "objectclass", "uid" });
//		    ctrls.setSearchScope(SearchControls.SUBTREE_SCOPE);
//		    
//		    NamingEnumeration<javax.naming.directory.SearchResult> answers = ctx.search("dc=ochem,dc=eu", "(uid=" + username + ")", ctrls);
//		    javax.naming.directory.SearchResult result = answers.nextElement();
//
//		    String user = result.getNameInNamespace();
//		    Properties props = new Properties();
//	        props.put(Context.INITIAL_CONTEXT_FACTORY, "com.sun.jndi.ldap.LdapCtxFactory");
//	        props.put(Context.PROVIDER_URL, url);
//	        props.put(Context.SECURITY_PRINCIPAL, user);
//	        props.put(Context.SECURITY_CREDENTIALS, password);
//		    
//		    ctx = new InitialDirContext(props);
//		    System.out.println(user + " was authenticated");
//	        ctx.close();
//			
//		} catch (AuthenticationNotSupportedException ex) {
//		    System.out.println("The authentication is not supported by the server");
//		    ex.printStackTrace();
//		    throw ex;
//		} catch (AuthenticationException ex) {
//		    System.out.println("incorrect password or username");
//		    ex.printStackTrace();
//		    throw ex;
//		} catch (NamingException ex) {
//		    System.out.println("error when trying to create the context");
//		    ex.printStackTrace();
//		    throw ex;
//		}
		
	}

}
