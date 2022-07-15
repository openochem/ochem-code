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

package qspr.controllers;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.xml.bind.annotation.XmlRootElement;

import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import com.eadmet.exceptions.UserFriendlyException;

import qspr.Globals;
import qspr.controllers.LoginController.LoginSettings;
import qspr.entities.OAuthApp;
import qspr.entities.OAuthCode;
import qspr.entities.OAuthToken;
import qspr.entities.Session;
import qspr.entities.User;
import qspr.frontend.WebModel;
import qspr.util.DynaWrap;

@Controller
public class OauthController extends ControllerWrapper {
	
	
	@XmlRootElement
	public static class AuthError {
		
		public String message;
	}
	
	private static ModelAndView error(String message, String template) {
//		TODO: make this more user friendly
		AuthError err = new AuthError();
		err.message = message;
		return new WebModel(err).setTemplate(template).getModelAndView();
	}

	public ModelAndView login(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		String client_id = req.getParameter("client_id");
		String redirect_uri = req.getParameter("redirect_uri");
		
		if (!OAuthApp.get(client_id).checkDomain(redirect_uri)) {
			return error("The domain of redirect URI is not whitelisted: " + redirect_uri, "oauth-authorize");
		}
		
		if (client_id != null && redirect_uri != null) {
			return redirect("login/show.do?oauth=true&" + req.getQueryString());
		} else {
			return error("Missing client ID or redirect URI.", "oauth-authorize");
		}
	}
	
	public ModelAndView authorize(HttpServletRequest req, HttpServletResponse res) throws Exception {
		Session session = Globals.userSession();
		if (session != null && session.user != null && req.getParameter("login") != null && req.getParameter("login").equals(session.user.login))
		{
			if (req.getMethod().equals("POST")) {
				String clientID = req.getParameter("clientID");
				String redirectUrl = req.getParameter("redirectURI");
				
				if (!OAuthApp.get(clientID).checkDomain(redirectUrl)) {
					return error("The domain of redirect URI is not whitelisted: " + redirectUrl, "oauth-authorize");
				}
				
				// FIXME: code should be short-lived -> add date of creation so that we can revoke it during authorization if used too late
				OAuthCode code = OAuthCode.create(OAuthApp.get(clientID), redirectUrl, session.user);
				code.save();
				return redirect(redirectUrl + "&code=" + code.toString());
			}
			LoginSettings settings = new LoginSettings();
			settings.oauth = new LoginSettings.OAuthSettings(req.getParameter("redirect_uri"), req.getParameter("client_id"));
			return new WebModel(settings).setTemplate("oauth-authorize").getModelAndView();
		} else {
			return error("You are not authorized to access this page.", "oauth-authorize");
		}
	}
	
	public ModelAndView token(HttpServletRequest req, HttpServletResponse res) throws Exception {
		String client_id = req.getParameter("client_id");
		String redirect_uri = req.getParameter("redirect_uri");
		
		if (!OAuthApp.get(client_id).checkDomain(redirect_uri)) {
			return error("The domain of redirect URI is not whitelisted: " + redirect_uri, "oauth-authorize");
		}
		
		String code = req.getParameter("code");
		if (req.getMethod().equals("GET") && client_id != null && redirect_uri != null && code != null) { // FIXME: should actually be POST, but even Facebook allows GET nowadays...
			OAuthCode code_ = OAuthCode.get(code);
			if (code_ == null) {
				return error("Wrong code.", "oauth-authorize");
			}
			if (!code_.app.clientID.equals(client_id)) {
				return error("Unknown client.", "oauth-authorize");
			}
			if (!code_.redirectURI.equals(redirect_uri)) {
				return error("Bad redirect URI.", "oauth-authorize");
			}
			
//			FIXME: add time of creation so we can revoke the token if too old
			OAuthToken token = OAuthToken.create(code_.app, code_.user);
			token.save();
			code_.delete();
			return new WebModel(token).getModelAndView();
		} else {
			return error("You are not authorized to access this page.", "oauth-authorize");
		}
	}
	
	@XmlRootElement
	public static class UserData {
		public String login;
		public String firstName;
		public String lastName;
	}
	
	public ModelAndView data(HttpServletRequest req, HttpServletResponse res) {
		String token = req.getParameter("access_token");
		String client_id = req.getParameter("client_id");
		if (token != null &&  client_id != null && isTokenValid(client_id, OAuthToken.get(token))) {
			OAuthToken token_ = OAuthToken.get(token);
			UserData data  = getUserData(token_.user);
			token_.delete();
			return new WebModel(data).getModelAndView();
		} else {
			return error("No access token received or token invalid.", "oauth-authorize");
		}
	}

	private static UserData getUserData(User user) {
		UserData data = new UserData();
		data.login = user.login;
		DynaWrap extended  = (user.isExtended())  ? user.dynaWrapped() :  null;
		data.firstName =  (extended !=  null)  && extended.getString("firstName") != null ? extended.getString("firstName") : "Not Set";
		data.lastName = (extended !=  null) && extended.getString("lastName") != null ? extended.getString("lastName") : "Not Set";
		return data;
	}

	private static boolean isTokenValid(String client_id, OAuthToken token) {
		if (token != null && token.app.clientID.equals(client_id)) {
			return true;
		}
		return false;
	}
	
	public ModelAndView create(HttpServletRequest req, HttpServletResponse res) {
		Session session = Globals.userSession();
		User user = session.user;
		if (session != null && user != null && req.getParameter("user") != null && req.getParameter("user").equals(user.login))
		{
			if (req.getMethod().equals("POST")) { // FIXME: add option to  delete
				String create = req.getParameter("create-new");
				String domainList = req.getParameter("domain-list");
				if (create != null && create.equals("true")) {
					try {
						OAuthApp new_app = OAuthApp.create(session.user);
						Set<String> domainSet = new HashSet<String>();
						domainSet.add("localhost");
						if (!domainList.equals("")) {
							domainSet.addAll(Arrays.asList(domainList.split(",")));	
						}
						new_app.domainList = String.join(",", domainSet);
						new_app.save();
					} catch (Exception e) {
						throw new UserFriendlyException(e);
					}
					return redirect("/user/profile.do?login=" + user.login);
				} else {
					throw new UserFriendlyException("Invalid data.");
				}
			}
			List<OAuthApp> list = OAuthApp.getAppsForUser(session.user);
			return new WebModel().setList(list).setTemplate("create-app").getModelAndView();
		} else {
			return error("You are not authorized to access this page.", "oauth-authorize");
		}
	}
	
	public ModelAndView delete(HttpServletRequest req, HttpServletResponse res) {
		Session session = Globals.userSession();
		User user = session.user;
		if (req.getMethod().equals("POST") && session != null && user != null && req.getParameter("clientID") != null)
		{
			String delete = req.getParameter("clientID");			
			OAuthApp to_delete = OAuthApp.get(delete);
			if (to_delete.user.id == user.id) {
				to_delete.delete();
				return redirect("/user/profile.do?login=" + user.login);
			} else {
				throw new UserFriendlyException("You can only delete your apps!");
			}
		} else {
			return error("You are not authorized to access this page.", "oauth-authorize");
		}
	}
}
