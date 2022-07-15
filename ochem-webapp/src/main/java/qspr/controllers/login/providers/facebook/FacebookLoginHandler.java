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

package qspr.controllers.login.providers.facebook;

import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;

import javax.servlet.http.HttpServletRequest;

import org.json.JSONException;
import org.json.JSONObject;

import com.eadmet.exceptions.UserFriendlyException;

import qspr.controllers.login.providers.LoginByProviderHandler;
import qspr.entities.User;
import qspr.util.StringUtils;
import qspr.util.URLUtil;

public class FacebookLoginHandler extends LoginByProviderHandler {
	
	public static class FacebookUser
	{
		public String id;
		public String userName;
		public String firstName;
		public String lastName;

		public FacebookUser(ProviderAccessToken accessToken) throws MalformedURLException, JSONException, IOException
		{
			JSONObject resp = new JSONObject(URLUtil.readURL(new URL("https://graph.facebook.com/me?access_token=" + accessToken.get())));
			id = resp.getString("id");
			userName = resp.getString("id");
			firstName = "Facebook User"; // FIXME: get the real name
			lastName = resp.getString("id");
		}
	}
	
	String apiKey;
	String secret;
	private static final String[] permissions = new String[] {"email"};
	
	public FacebookLoginHandler(FacebookProviderConfig config) {
		super(config);
		apiKey = config.getKey();
		secret = config.getSecret();
	}

	@Override
	public String getProviderLoginURL() {
		return "https://www.facebook.com/v12.0/dialog/oauth?client_id=" +
				apiKey + "&state=123456&scope=" + StringUtils.arrayToString(permissions) + "&redirect_uri=" + getOCHemRedirectURL(true);
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
		return "https://graph.facebook.com/v12.0/oauth/access_token?client_id=" +
				apiKey+"&redirect_uri=" + getOCHemRedirectURL(true) +
				"&client_secret="+secret+"&code=" + callbackData.getCode();
	}

	@Override
	protected ProviderAccessToken parseAuthResponse(JSONObject data) {
		try {
			String accessToken = data.getString("access_token");
//			Integer expires = data.getInt("expires_in");
			return new ProviderAccessToken(accessToken);
		} catch (JSONException e) {
			throw new UserFriendlyException(e);
		}
	}
	
	protected FacebookUser getFBUser(ProviderAccessToken token) {
		try {
			return new FacebookUser(token);
		} catch (JSONException | IOException e) {
			throw new UserFriendlyException(e);
		}
	}

	@Override
	protected User getUser(ProviderAccessToken token) {
		// FIXME: update first and last name
		return User.getByFacebookId(getFBUser(token).userName); //FIXME: change this to get by provider ID -> each provider will have a user ID
	}

	@Override
	protected User createUser(ProviderAccessToken token) throws UserFriendlyException {
		// FIXME: save first and last name
		FacebookUser fbUser = getFBUser(token);
		User sameLoginUser = User.getByLogin(fbUser.userName);
		if (sameLoginUser != null)
			throw new UserFriendlyException("Cannot create a new account from your facebook profile, since userID " + sameLoginUser.login + " already exists. If its you, please login or request a forgotten password.");
		
		User user = User.getNewInstance();
		user.login = fbUser.userName;
		user.facebookId = fbUser.userName; // FIXME: change to providerID and add a field for provider itself
		if (user.isExtended()) {
			try {
				User.getExtended().getMethod("setPassword", int.class).invoke(user, 8);
			} catch (Exception e) {
				throw new UserFriendlyException(e);
			}
		}
		user.rank = User.NOTVALIDATED;
		
		return user;
	}

	

}
