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

import java.io.IOException;
import java.net.URLEncoder;
import java.util.ArrayList;
import java.util.List;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.servlet.http.HttpSession;
import javax.xml.bind.JAXBException;
import javax.xml.bind.annotation.XmlRootElement;

import org.hibernate.criterion.Restrictions;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.ThreadScope;
import qspr.controllers.login.providers.LoginByProviderHandler;
import qspr.controllers.login.providers.LoginHandlers;
import qspr.controllers.login.providers.ldap.LDAPLoginHandler;
import qspr.controllers.login.providers.ldap.LDAPProviderConfig;
import qspr.entities.Basket;
import qspr.entities.Session;
import qspr.entities.User;
import qspr.frontend.WebModel;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.useractions.EventFactory;
import com.eadmet.useractions.LoginAction;
import com.eadmet.utils.mailer.Mailer;

@Controller
public class LoginController extends ControllerWrapper 
{
	@SuppressWarnings("unchecked")
	public ModelAndView login(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		String login = req.getParameter("login");
		String pwd   = req.getParameter("pwd");
		
		HttpSession httpSession = req.getSession(true);
		httpSession.setMaxInactiveInterval(12*60*60);
		String targetPage = (String) req.getSession().getAttribute("requested-page");
		if (targetPage == null)
			targetPage = "/home/show.do";
		targetPage =  targetPage.replaceFirst(OCHEMConfiguration.rootDir, "").replace("//", "/");

		if (assertParam(QSPRConstants.ANONYMOUS))
		{
			checkAnonymousUsers();
			loginUser(null);
			return redirect("login/acceptLicenseAgreement.do");
		}

		List<User> User = Globals.session().createQuery(
				//"from User where login=:login or email=:login").setString("login", login)
				"from User where login=:login").setString("login", login) // There is no email unless ExtendedUser
				.list();
		if (User.size() > 0)
		{
			User user = User.get(0);
			if  (!(user.isExtended())) {
				return redirect("static/pwd-error.do?wrong-login");
			}
			
			if (user.authorizeWithSecret(pwd)) 
			{
				loginUser(user);
				EventFactory.document("User login", new LoginAction(req.getHeader("User-agent"), ThreadScope.resolveRemoteAddr(req), null));
				if (req.getParameter("oauth") != null && req.getParameter("oauth").equals("true")) {
					String clientID = req.getParameter("clientID");
					String redirectURI = req.getParameter("redirectURI");
					return redirect("oauth/authorize.do?client_id=" + clientID + "&redirect_uri=" + URLEncoder.encode(redirectURI) + "&login=" + user.login);
				}
				return redirect(targetPage.replace("render-mode=popup", "render-mode=full"));
			}
			else
				return redirect("static/pwd-error.do?wrong-login");
		} 
		return redirect("static/pwd-error.do?wrong-password");
	}

	void checkAnonymousUsers(){
		if(OCHEMConfiguration.disableAnonymousUsers) {
			Globals.userSession(null); // Go away
			throw new UserFriendlyException("Anonymous users are disabled for this installation. Please, register.");
		}
	}
	
	public ModelAndView ldap(HttpServletRequest req, HttpServletResponse res) {
		LDAPLoginHandler handler = (LDAPLoginHandler) LoginHandlers.get("LDAP");
		LDAPProviderConfig config = (LDAPProviderConfig) handler.getConfig();
		if ( req.getMethod().equals("POST") ) {
			return handler.callbackController(req, res, this);
		} else {
			return new WebModel().setList(new ArrayList<>(config.getServers())).setTemplate("ldap").getModelAndView();
		}
	}
	
	public ModelAndView chooseProvider(HttpServletRequest req, HttpServletResponse res) {
		return new WebModel().setList(LoginHandlers.getInfos()).setTemplate("login-provider").getModelAndView();
	}
	
	public ModelAndView withProvider(HttpServletRequest req, HttpServletResponse res)
	{
		
		if (!assertParam("provider"))
			throw new UserFriendlyException("You did not select an identity provider.");
		
		String provider = req.getParameter("provider");
		LoginByProviderHandler handler = LoginHandlers.get(provider);
		if (handler != null) {
			return redirect(handler.getProviderLoginURL());
		} else {
			throw new UserFriendlyException("Unknown identity provider: " + provider);
		}
	}

	public ModelAndView providerCallback(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		if (!assertParam("provider"))
			throw new UserFriendlyException("No identity provider selected.");
		
		String provider = req.getParameter("provider");
		LoginByProviderHandler handler = LoginHandlers.get(provider);
		if (handler != null) {
			return handler.callbackController(req, res, this);
		} else {
			throw new UserFriendlyException("Unknown identity provider: " + provider);
		}
	}

	@SuppressWarnings("unchecked")
	public ModelAndView logout(HttpServletRequest req, HttpServletResponse res) throws IOException, JAXBException
	{
		if (Globals.userSession() != null && Globals.userSession().user == null)
		{
			// It is an anonymous user who is logging out
			// Delete his basket entries
			List<Basket> bes = Globals.session()
					.createCriteria(Basket.class)
					.add(Restrictions.eq("session", Globals.userSession()))
					.list();
			for (Basket basket : bes) 
				Globals.session().delete(basket);
		}

		EventFactory.document("Logout", null, "has logged out");

		Globals.userSession(null);
		req.getSession().invalidate();

		return redirect("login/show.do");
	}

	public ModelAndView acceptLicenseAgreement(HttpServletRequest req, HttpServletResponse res)
	{
		Session session = Globals.userSession();
		if(session.user == null || assertParam(QSPRConstants.ANONYMOUS) || assertParam("guest"))checkAnonymousUsers();
		if (assertParam("accepted") || session.license || OCHEMConfiguration.inhouseInstallation)
		{
			session.license = true;
			if (session.user != null && session.user.licenseVersion < Globals.CURRENT_LICENSE_VERSION)
			{
				// Update the accepted license version for this user
				session.user.licenseVersion = Globals.CURRENT_LICENSE_VERSION;
				Globals.session().saveOrUpdate(session.user);
			}
			String targetPage = (String) req.getSession().getAttribute("requested-page");
			if (targetPage == null)
				targetPage = "home/show.do";
			targetPage =  targetPage.replaceFirst(OCHEMConfiguration.rootDir, "");

			// After login, we load the target page in the global frame, so use the full render mode
			targetPage =  targetPage.replace("render-mode=popup", "render-mode=full");

			// create the selectionBasket for this user when (s)he logs in. --> selected records are shown in basket filter
			Globals.selectionBasket(false); 

			return redirect(targetPage);
		}
		else if (assertParam("denied"))
		{
			if (session.user != null)
				Mailer.notifyAdmins("User " + session.user.login + " has rejected the license agreement", "");
			Globals.userSession(null); // Go away
			return redirect("home/show.do");
		}
		else
			return new WebModel().setTemplate("accept-license").getModelAndView();
	}

	public ModelAndView anonymous(HttpServletRequest req, HttpServletResponse res) throws IOException
	{
//		checkAnonymousUsers();
		return new WebModel().setTemplate("login").getModelAndView();
	}
	
	@XmlRootElement
	public static class LoginSettings {
		public boolean guestLogin = OCHEMConfiguration.guestLogin;
		public boolean providerLogin = OCHEMConfiguration.providerLogin;
		public boolean registerLogin = OCHEMConfiguration.registerLogin;
	
		public OAuthSettings oauth;
		
		@XmlRootElement
		public static class OAuthSettings {
			public String redirectURI;
			public String clientID;
			
			public OAuthSettings() {
				// no action
			}
			
			public OAuthSettings(String uri, String client) {
				this.redirectURI = uri;
				this.clientID = client;
			}
		}
		
	}

	public ModelAndView show(HttpServletRequest req, HttpServletResponse res) throws IOException
	{
		LoginSettings settings = new LoginSettings();
		if (req.getParameter("oauth") != null && req.getParameter("oauth").equals("true")) {
			settings.guestLogin = false;
			settings.providerLogin = false;
			settings.oauth = new LoginSettings.OAuthSettings(
					req.getParameter("redirect_uri"), req.getParameter("client_id")
			);
		}
		return new WebModel(settings).setTemplate("login").getModelAndView();
	}

}
