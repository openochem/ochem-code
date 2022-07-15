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

import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLEncoder;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.json.JSONException;
import org.json.JSONObject;
import org.springframework.web.servlet.ModelAndView;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.useractions.EventFactory;
import com.eadmet.useractions.LoginAction;
import com.eadmet.utils.mailer.Email;
import com.eadmet.utils.mailer.Mailer;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.ThreadScope;
import qspr.controllers.LoginController;
import qspr.entities.User;
import qspr.util.DynaWrap;
import qspr.util.URLUtil;

public abstract class LoginByProviderHandler {
	public abstract String getProviderLoginURL();
	protected abstract ProviderCallbackData parseCallbackRequest(HttpServletRequest req);
	protected abstract String getAuthURL(ProviderCallbackData callbackData);
	protected abstract ProviderAccessToken parseAuthResponse(JSONObject data);
	protected abstract User getUser(ProviderAccessToken token);
	protected abstract User createUser(ProviderAccessToken token) throws UserFriendlyException;
	
	public ProviderConfig config;
	private static transient final Logger logger = LogManager.getLogger(LoginController.class);
	
	public LoginByProviderHandler(ProviderConfig config) {
		this.config = config;
	}
	
	protected class ProviderCallbackData {
		private List<String> errors = new ArrayList<>();
		private String code = null;
		
		public ProviderCallbackData() {
			// no action
		};
		
		public void addError(String error) {
			errors.add(error);
		}

		public List<String> getErrors() {
			return errors;
		}

		public String getCode() {
			return code;
		}

		public void setCode(String code) {
			this.code = code;
		};
		
		public boolean isValid() {
			return errors.isEmpty() && code != null;
		}
		
	}
	
	protected class ProviderAccessToken {
		private String token = null;
		
		public ProviderAccessToken(String token) {
			this.token = token;
		}

		public String get() {
			return token;
		}		
	}
	
	private ProviderAccessToken getAccessToken(ProviderCallbackData callbackData) {
		String authURL = getAuthURL(callbackData);
		URL url;
		try {
			url = new URL(authURL);
		} catch (MalformedURLException e) {
			throw new UserFriendlyException(e);
		}
		String result;
		try {
			result = URLUtil.readURL(url);
		} catch (IOException e) {
			throw new UserFriendlyException("Error getting access token.");
		}
		JSONObject jsonObj;
		try {
			jsonObj = new JSONObject(result);
		} catch (JSONException e) {
			throw new UserFriendlyException(e);
		}
		return parseAuthResponse(jsonObj);
	};
	
	public ModelAndView callbackController(HttpServletRequest req, HttpServletResponse res, LoginController caller) {
		ProviderCallbackData data = parseCallbackRequest(req);
		if (!data.getErrors().isEmpty()) {
			return caller.redirect(getOChemStartURL() + "&errors=" + String.join(",", data.getErrors()));
		} else if (data.getCode() == null) {
			return caller.redirect(getOChemStartURL() + "&errors=No authorization code obtained in the request.");
		}
		
		ProviderAccessToken accessToken = getAccessToken(data);
		if (accessToken == null) {
			return caller.redirect(getOChemStartURL() + "&errors=No token obtained. Log in failed.");
		}
		
		User user = getUser(accessToken);
		if (user != null) {
			caller.loginUser(user);
			logUserLogin(user, req);
		} else {
			user = createUser(accessToken);
			Globals.session().save(user);
			caller.loginUser(user);
			logUserCreate(user, req);
			logUserLogin(user, req);
		}
		
		return caller.redirect(getOChemStartURL() + "&login=" + user.login);
	};
	
	public String getOCHemRedirectURL() {
		return OCHEMConfiguration.getRootHost() + "/login/providerCallback.do?provider=" + config.getURLName() + "&render-mode=popup";
	}
	
	public String getOCHemRedirectURL(boolean encoded) {
		String redirectURI = getOCHemRedirectURL();
		if (!encoded) return redirectURI;
		
		try {
			redirectURI = URLEncoder.encode(redirectURI, StandardCharsets.UTF_8.toString());
		} catch (UnsupportedEncodingException e) {
			throw new UserFriendlyException(e);
		}
		return redirectURI;
	}
	
	public String getOChemStartURL() {
		return OCHEMConfiguration.getRootHost() + "/login/chooseProvider.do?provider=" + config.getURLName() + "&render-mode=popup";
	}
	
	public ProviderConfig getConfig() {
		return config;
	}
	
	protected void logUserLogin(User user, HttpServletRequest req) {
		logger.info("User " + user.login + " has logged in with provider: " + config.getURLName());
		EventFactory.document("User login", new LoginAction(req.getHeader("User-agent"), ThreadScope.resolveRemoteAddr(req), config.getURLName()));
	}
	
	protected void logUserCreate(User user, HttpServletRequest req) {
		logger.info("User " + user.login + " created from " + config.getURLName());
		Mailer.notifyAdmins("New user has just registered with an identity provider:" + config.getURLName(), "A new user, " + user.login + ", has just registered with their" + config.getURLName() + " account!\nThe lucky person's new login is " + user.login + ".\n\nI will send them a confirmation email.");
		// Send the user a warm welcome
		if (user.isExtended()) {
			try {
					DynaWrap extended = user.dynaWrapped();
					Mailer.postMailAsync(new Email(extended.getString("email"), "Welcome to OCHEM!", String.format("Dear %s %s, \n\nThank you for  registering at OCHEM. You account has beed automatically created using your " + config.getURLName() + " data.\n\n"+
							"Your login data is:\nUsername: %s\nPassword (autogenerated): %s\n\nYou can always login into OCHEM via Facebook, without entering your OCHEM password.\n\nKind regards,\nOCHEM Team",
							extended.getString("firstName"), extended.getString("lastName"), user.login, extended.getString("password"))));
				
			} catch (Exception e) {
				throw new UserFriendlyException(e);
			}
		}
	}
}
