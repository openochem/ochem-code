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

package qspr.controllers.login.providers.ochem;

import java.io.IOException;
import java.net.URL;

import javax.servlet.http.HttpServletRequest;

import org.json.JSONException;
import org.json.JSONObject;

import com.eadmet.exceptions.UserFriendlyException;

import qspr.Globals;
import qspr.controllers.login.providers.LoginByProviderHandler;
import qspr.entities.User;
import qspr.util.DynaWrap;
import qspr.util.URLUtil;

public class OCHEMLoginHandler extends LoginByProviderHandler {
	
	private String apiKey;
	private String secret;
	private String mainOCHEMURL;

	public OCHEMLoginHandler(OCHEMLoginConfig config) {
		super(config);
		apiKey = config.apiKey;
		secret = config.secret;
		mainOCHEMURL = config.mainOCHEMURL.replaceAll("/+$", "");
	}

	@Override
	public String getProviderLoginURL() {
		return mainOCHEMURL + "/oauth/login.do?client_id=" + apiKey + "&redirect_uri=" + getOCHemRedirectURL(true);
	}

	@Override
	protected ProviderCallbackData parseCallbackRequest(HttpServletRequest req) {
		ProviderCallbackData data = new ProviderCallbackData();
		String error = req.getParameter("error");
		String code = req.getParameter("code");
		if (error != null) {
			data.addError(error);
		}
		if (code != null) {
			data.setCode(code);
		}
		return data;
	}

	@Override
	protected String getAuthURL(ProviderCallbackData callbackData) {
		return mainOCHEMURL + "/oauth/token.do?client_id=" +
				apiKey+"&redirect_uri=" + getOCHemRedirectURL(true) +
				"&client_secret="+secret+"&code=" + callbackData.getCode() + "&out=json";
	}

	@Override
	protected ProviderAccessToken parseAuthResponse(JSONObject data) {
		JSONObject obj;
		try {
			obj = data.getJSONObject("token");
			String accessToken = obj.getString("value");
			return new ProviderAccessToken(accessToken);
		} catch (JSONException e) {
			throw new UserFriendlyException(e);
		}
	}
	
	private static class UserData {
		public String login;
		public String firstName;
		public String lastName;
	}
	
	private String getUserURL(ProviderAccessToken token) {
		return mainOCHEMURL + "/oauth/data.do?access_token=" +
				token.get() + "&client_id=" + apiKey +  "&out=json";
	}
	
	private UserData fetchUserData(ProviderAccessToken token) {
		try {
			URL url = new URL(getUserURL(token));
			String result = URLUtil.readURL(url);
			JSONObject resp = new JSONObject(result);
			// TODO: process errors
			JSONObject user = resp.getJSONObject("userData");
			UserData data = new UserData();
			data.login = user.getString("login");
			data.firstName = user.getString("firstName");
			data.lastName = user.getString("lastName");
			return data;
		} catch (IOException | JSONException e) {
			throw new UserFriendlyException(e);
		}
	}

	@Override
	protected User getUser(ProviderAccessToken token) {
		UserData data = fetchUserData(token);
		User user;
		try {
			user = User.getByLogin(data.login);
			if (user == null) {
				throw new Exception("User not found.");
			}
			
			DynaWrap extended = null;
			if (user.isExtended()) {
				extended = user.dynaWrapped();
			}
			if (extended != null && (!extended.getString("firstName").equals(data.firstName) || !extended.getString("lastName").equals(data.lastName))) {
				extended.setField("firstName", data.firstName);
				extended.setField("lastName", data.lastName);
				Globals.session().save(user);
			}
			return user;
		} catch (Exception exp) {
			return createUser(data);
		}
	}
	
	protected User createUser(UserData data) {
		User sameLoginUser = User.getByLogin(data.login);
		if (sameLoginUser != null)
			throw new UserFriendlyException("Cannot create a new account from your main OCHEM profile, since userID " + sameLoginUser.login + "already exists.");
		
		User user = User.getNewInstance();
		user.login = data.login;
		
		DynaWrap extended = null;
		if (user.isExtended()) {
			extended = user.dynaWrapped();
		}
		if  (extended != null) {
			extended.setField("firstName", data.firstName);
			extended.setField("lastName", data.lastName);
			try {
				User.getExtended().getMethod("setPassword", int.class).invoke(user, 8);
			} catch (Exception e) {
				throw new UserFriendlyException(e);
			}
		}
		user.rank = User.NOTVALIDATED;
		
		return user;
	}

	@Override
	protected User createUser(ProviderAccessToken token) throws UserFriendlyException {
		UserData data = fetchUserData(token);
		return createUser(data);
	}

}
