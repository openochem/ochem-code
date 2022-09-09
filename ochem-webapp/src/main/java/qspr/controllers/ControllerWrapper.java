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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.Filter;
import org.hibernate.criterion.Order;
import org.springframework.web.multipart.MultipartHttpServletRequest;
import org.springframework.web.servlet.ModelAndView;
import org.springframework.web.servlet.mvc.multiaction.MultiActionController;
import org.springframework.web.servlet.mvc.multiaction.NoSuchRequestHandlingMethodException;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.ThreadScope;
import qspr.entities.Alert;
import qspr.entities.Announcement;
import qspr.entities.Session;
import qspr.entities.User;
import qspr.export.ExportableModel;
import qspr.frontend.WebModel;
import qspr.util.Operation;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.FileUtils;
import com.eadmet.utils.mailer.Mailer;

public class ControllerWrapper extends MultiActionController {
	protected boolean sessionRequired;
	protected boolean useIframe = true;
	public String outerTemplate;

	static Map<String, Integer> runningRequest = new HashMap<String, Integer>();

	public void cleanRequest(String method) {
		String request = getProfile(method);
		if(request == null) return;
		int num = runningRequest.get(request);
		runningRequest.put(request,num-1);
	}

	boolean startNewRequest(String request, int maxrequests){
		if(!runningRequest.containsKey(request))
			runningRequest.put(request,0);
		int num = runningRequest.get(request);
		if(num >= maxrequests) return false;
		runningRequest.put(request,num+1);
		return true;
	}

	String getProfile(String method) {
		if (Globals.userSession().user == null) return null;
		return  Globals.userSession().user + "_"+ method;
	}

	public String canExecuteLongRequest(String method,int maxrequests) {
		String request = getProfile(method);
		if(request == null) return null;

		boolean waited = false;
		int i, maxn = 100;
		for(i =0; !startNewRequest(request,maxrequests) && i<maxn; i++) {
			try {
				Thread.sleep(1000);
			}catch (Exception e) {
			}
			waited = true;
			logger.info("waiting to start: " + request + " " + i);
		}
		if (i >= maxn) {
			logger.info("cancelling: " + request);
			return null;
		}
		if(waited)
			logger.info("starting after waiting: " + request);
		return request;
	}

	public void setSessionRequired(boolean something) {
		sessionRequired = something;
	}

	protected String getRedirectCheckConsensus(ExportableModel template){
		return template.method.equals(QSPRConstants.CONSENSUS)?"modelconfigurator/choose.do?render-mode=popup":"modelconfigurator/choose.do?render-mode=popup&skip-configuration=1";
	}

	@SuppressWarnings("unchecked")
	protected ModelAndView handleRequestInternal(HttpServletRequest request, HttpServletResponse response) throws Exception {
		ThreadScope.get().localRequest = request;
		ThreadScope.get().localMpRequest = null;

		try {
			MultipartHttpServletRequest mpr = (MultipartHttpServletRequest) request;
			ThreadScope.get().localMpRequest = mpr;
			logger.debug("Multipart request received");
		} catch (Exception e) {
			logger.debug("Non-multipart request received");
		}

		Session session = null;
		boolean nodb = request.getParameter("nodb") != null;

		if (assertParam("referral"))
			request.getSession().setAttribute("referral", getParam("referral"));

		if (!nodb) {
			// Create a guest session it it has not been yet created, for anonymous
			session = Globals.userSession();

			if (session != null && session.user != null && Globals.uniqueLoginRequired) {
				Session ls = Session.getLastSession(session.user.login);
				if (ls.id > session.id) // If we have a logged in session after the current one, invalidate the current one
				{
					Globals.userSession(null);
					request.getSession().invalidate();
				}
			}

			if (sessionRequired && session == null && !handlerAnnotatedWith(NoLoginRequired.class, request)) {
				// Ease for Developers: Login automatically / Midnighter
				if (OCHEMConfiguration.autoLoginUser != null && !OCHEMConfiguration.autoLoginUser.isEmpty())
					loginUser(User.getByLogin(OCHEMConfiguration.autoLoginUser));
				else {
					if (assertParam("out")) // Ajax requests in browsers
						return new WebModel(Alert.Error("friendly",
								"You are not logged in. Probably you have been logged out by timeout. Please log in again to proceed.")).getModelAndView();

					setTargetPage(request);
					return redirect("login/show.do");
				}
			}
		}

		Filter filter;
		if (session != null) {

			if (!session.license && !request.getRequestURI().contains("acceptLicenseAgreement")) {
				setTargetPage(request);
				return redirect("login/acceptLicenseAgreement.do");
			}
//			if (session.user != null && session.user.organisation == null && !request.getRequestURI().contains("user") && !request.getRequestURI().contains("login")) {
//				setTargetPage(request);
//				return redirect("user/show.do");
//			}
			if (session.user == null) {
				filter = Globals.session().enableFilter("sessionFilter");
				filter.setParameter("sessionId", session.id.longValue());
				// logger.info("Session filter applied "+session.id);
			} else {
				filter = Globals.session().enableFilter("userFilter");
				filter.setParameter("userId", session.user.id.longValue());
				// logger.info("User filter has been apllied "+session.user.id);
			}

			if (assertParam("developer"))
				Globals.userSession().isDeveloperSession = true; // Enable the features visible only to the developers

			// logger.info("[CW] Entered session id "+session.id);
		}

		String[] pieces = this.getClass().getName().split("\\.");
		ThreadScope.get().controller = pieces[pieces.length - 1].toLowerCase().replace("controller", "");

		ModelAndView mav;
		try {
			if (session != null && session.user != null && session.user.isSuspended())
				throw new UserFriendlyException("The user (" + session.user.login + ") has been suspended. Please, contact the administrator: " + QSPRConstants.INFOEMAIL);
			if (isIFrame() || "redirect".equals(getParam("render-mode"))) {
				// Show content of iframe directly

				if (request.getParameter("operation-id") != null && (request.getParameter("out") == null || request.getParameter("start-operation") != null)
						&& Operation.getOperation(getLongParam("operation-id")) == null && !request.getRequestURI().contains("operationWaitingScreen"))
					ThreadScope.get().operation = new Operation(getLongParam("operation-id"));

				if (session != null && session.user != null && session.user.isTestUser())
					Mailer.enableForThread(false);
				mav = super.handleRequestInternal(request, response);
			} else {
				// Show global template with iframe
				String requestedUrl = request.getRequestURL().toString();
				if (request.getQueryString() != null)
					requestedUrl += "?" + request.getQueryString();
				requestedUrl = requestedUrl.replace("render-mode=full", "");

				WebModel wmOuter = new WebModel(requestedUrl);
				wmOuter.outerTemplate = outerTemplate;
				wmOuter.setTemplate("faq");

				mav = (wmOuter).getModelAndView();
			}
			if (!nodb && Globals.isMainTransactionRunning())
				Globals.session().flush();
		}
		catch (UserFriendlyException ee) {
			logger.info("UserFriendlyException: " + ee.getMessage());
			if (!nodb)
				Globals.rollbackMainTransaction();
			return new WebModel(Alert.Exception(ee)).setTemplate("exception").setRenderMode("popup").getModelAndView();
		}
		catch (Exception e) {
			logger.error("Managed exception", e);
			if (!nodb)
				Globals.rollbackMainTransaction();
			return new WebModel(Alert.Exception(e)).setTemplate("exception").setRenderMode("popup").getModelAndView();
		}
		finally
		{
			Mailer.enableForThread(true);
		}

		if (mav != null) {
			if (mav.getModel().get("object") != null && mav.getModel().get("object") instanceof WebModel)
			{
				WebModel wm = (WebModel) mav.getModel().get("object");
				if (wm != null) {
					wm.session = session;
					if (wm.templateName == null)
						wm.templateName = this.getClass().getName().toLowerCase().replaceAll(".*\\.([^.]*)controller", "$1");
					wm.versionInfo = getVersionInfo();

					if (!nodb) {
						Criteria c = Globals.session().createCriteria(Announcement.class).addOrder(Order.desc("id")).setMaxResults(1);
						List<Announcement> l = c.list();
						if (l.size() > 0)
							wm.announcement = l.get(0);
					}

					if (isIFrame())
						wm.setRenderMode("popup");
				}
			}
		}
		return mav;
	}

	protected void setTargetPage(HttpServletRequest request) {
		if (!request.getRequestURI().contains("user") && !request.getRequestURI().contains("login") && !"json".equals(request.getParameter("out"))) {
			request.getSession().setAttribute("requested-page", getRequestURI(request));
		}
	}

	protected boolean isIFrame() {
		// This method is a set of work-arounds
		// Subject for refactoring: find a robust way of determination
		// whether current page is in IFRAME or not
		// Midnighter

		if (!useIframe)
			return true;

		HttpServletRequest request = ThreadScope.get().localRequest;
		if (ThreadScope.get().controller != null)
		{
			if (ThreadScope.get().controller.toLowerCase().equals("physprop"))
				return true;

			if (ThreadScope.get().controller.toLowerCase().equals("marvinintegration"))
				return true;



			try
			{
				if (handlerAnnotatedWith(NoFrame.class, request))
					return true;
				if (ThreadScope.get().controller.toLowerCase().equals("home") && "index".equals(getMethodNameResolver().getHandlerMethodName(request)))
					return true;
			}
			catch (NoSuchRequestHandlingMethodException e)
			{
				// ...
			} catch (SecurityException e)
			{
				// TODO Auto-generated catch block
				e.printStackTrace();
			} 
		}


		if ("full".equals(getParam("render-mode")) || "redirect".equals(getParam("render-mode")))
			return false;

		return ((request.getHeader("Referer") != null && request.getHeader("Referer").contains("render-mode=popup")) || "popup".equals(getParam("render-mode"))
				|| "iprior".equals(getParam("render-mode")) || "json".equals(getParam("out")) || "xml".equals(getParam("out")) || request.getMethod().equals(
						"POST"));
	}

	@SuppressWarnings({ "rawtypes", "unchecked" })
	private boolean handlerAnnotatedWith(Class clazz, HttpServletRequest request) {
		try
		{
			return this.getClass().getMethod(getMethodNameResolver().getHandlerMethodName(request), HttpServletRequest.class, HttpServletResponse.class).getAnnotation(clazz) != null;
		} catch (Exception e)
		{
			e.printStackTrace();
			return false;
		}
	}

	protected String getVersionInfo() throws IOException {
		// Do we have the release information file
		File releaseInfo = new File(getServletContext().getRealPath("/") + "release.info");
		if (releaseInfo.exists())
			return "v." + FileUtils.getFileAsString(releaseInfo.getAbsolutePath());

		// Do we have the version.info file?
		File versionInfo = new File(getServletContext().getRealPath("/") + "version.info");
		if (versionInfo.exists()) {
			BufferedReader br = new BufferedReader(new FileReader(versionInfo));
			String revision = br.readLine();
			String owner = br.readLine();
			String checkInDate = br.readLine();
			String buildDate = br.readLine();
			String buildIp = br.readLine();
			br.close();
			return "Revision " + revision + " by " + owner + " checked in on " + checkInDate + ". Built from " + buildIp + " on " + buildDate;
		}

		return null;

	}

	public boolean assertParam(String name) {
		MultipartHttpServletRequest mpRequest = ThreadScope.get().localMpRequest;
		if (mpRequest == null)
			return (request().getParameter(name) != null && !request().getParameter(name).equals("") && !request().getParameter(name).equals("undefined"));
		else
			return (mpRequest.getParameter(name) != null && !mpRequest.getParameter(name).equals("") && !mpRequest.getParameter(name).equals("undefined"));
	}

	public Long getLongParam(String name) {
		if (!assertParam(name))
			return null;
		if (ThreadScope.get().localMpRequest == null) {
			String param = ThreadScope.get().localRequest.getParameter(name);
			if (param != null)
				return Long.valueOf(param.trim());
			else
				return null;
		} else {
			String param = ThreadScope.get().localMpRequest.getParameter(name);
			if (param != null)
				return Long.valueOf(param.trim());
			else
				return null;
		}
	}

	public Integer getIntParam(String name) {
		if (!assertParam(name))
			return null;
		if (ThreadScope.get().localMpRequest == null)
			//FIXME It is work around to avoid problem with Integer exception
			// If there is more than 1 provider, two  or more int values are provided for each provider, but there is only one return value!
			// Exception because Integer cannot parse string with ore than one int value
			// Fix is to parse the first value only.
			return Integer.valueOf(ThreadScope.get().localRequest.getParameter(name.trim()).split("\\s+")[0]);
		else
			return Integer.valueOf(ThreadScope.get().localMpRequest.getParameter(name.trim()));
	}

	public Double getDoubleParam(String name) {
		if (!assertParam(name))
			return null;
		if (ThreadScope.get().localMpRequest == null)
			return Double.valueOf(ThreadScope.get().localRequest.getParameter(name.trim()));
		else
			return Double.valueOf(ThreadScope.get().localMpRequest.getParameter(name.trim()));
	}

	public String getParamTrim(String name) {
		String s = getParam(name);
		return s == null?null:s.trim();
	}

	public String getParam(String name) {
		if (ThreadScope.get().localMpRequest == null)
			return ThreadScope.get().localRequest.getParameter(name);
		else
			return ThreadScope.get().localMpRequest.getParameter(name);
	}

	public void loginUser(User user) {
		if (Globals.userSession() == null || Globals.userSession().user != user) {
			if (user != null)
				logger.info("Logging in " + user.login + "...");
			else
				logger.info("Logging in as guest...");

			Session newUserSession = Session.createNewSession(user);
			newUserSession.setIpAddress(ThreadScope.resolveRemoteAddr(ThreadScope.get().localRequest));
			newUserSession.userAgent = ThreadScope.get().localRequest.getHeader("user-agent");
			Globals.session().save(newUserSession);
			Globals.userSession(newUserSession);
		}
	}

	HttpServletRequest request() {
		return ThreadScope.get().localRequest;
	}

	public void assertNotNull(Object object, String name) throws Exception {
		if (object == null)
			throw new Exception("The field \"" + name + "\" should be defined!");
	}

	public ModelAndView redirect(String url) {
		String chr = url.contains("?") ? "&" : "?";
		if (!url.contains("render-mode"))
		{
			if (isIFrame())
				url += chr + "render-mode=popup";
			else
				url += chr + "render-mode=full";
		}

		if (url.startsWith("http"))
			return new ModelAndView("redirect:" + url);
		else
			return new ModelAndView("redirect:" + OCHEMConfiguration.getRootHost() + OCHEMConfiguration.rootDir + "/" + url);
	}

	public static String now() {
		Calendar cal = Calendar.getInstance();
		SimpleDateFormat sdf = new SimpleDateFormat("dd-MM-yyyy HH:mm:ss");
		return sdf.format(cal.getTime());
	}

	@SuppressWarnings("unchecked")
	protected String getRequestURI(HttpServletRequest request) {
		Enumeration<String> params = request.getParameterNames();
		StringBuilder data = new StringBuilder();
		while (params.hasMoreElements()) {
			String key = params.nextElement();
			String value = request.getParameter(key);
			data.append("&" + key + "=" + value);
		}

		if (data.length() > 0)
			return request.getRequestURI() + "?" + data.substring(1);
		else
			return request.getRequestURI();
	}

	protected String[] getURIParts(HttpServletRequest request)
	{
		String[] parts = request.getRequestURI().substring(1).split("/");
		for (int i = 0; i < parts.length; i++)
			if (parts[i].contains(".do"))
				parts[i] = parts[i].substring(0, parts[i].length() - 3);

		return parts;
	}

	public void setStatus(String status) {
		if (ThreadScope.get().operation != null)
			ThreadScope.get().operation.setStatus(status);
	}

	private static Logger logger = LogManager.getLogger(ControllerWrapper.class);
}
