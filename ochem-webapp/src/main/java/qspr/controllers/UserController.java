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

import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.List;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.hibernate.Criteria;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Restrictions;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.MarshallingOption;
import qspr.OCHEMConfiguration;
import qspr.business.UserProfile;
import qspr.business.WebFilters;
import qspr.entities.Invite;
import qspr.entities.Session;
import qspr.entities.User;
import qspr.frontend.BrowserModel;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;
import qspr.util.AccessChecker;
import qspr.util.DynaWrap;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.mailer.Email;
import com.eadmet.utils.mailer.Mailer;

@Controller
public class UserController extends BrowserWrapper
{
	User sessionUser;

	public ModelAndView pleaseregister(HttpServletRequest request, HttpServletResponse response) throws Exception 
	{
		return new WebModel().setTemplate("please-register").getModelAndView();
	}

	public ModelAndView show(HttpServletRequest request, HttpServletResponse response) throws Exception 
	{
		init(request);		

		if (sessionUser == null)
			return redirect("login/show.do");

		ModelAndView mav = new ModelAndView("xslt");
		mav.addObject("object", new WebModel(sessionUser));
		return mav;
	}	

	public ModelAndView newuser(HttpServletRequest request, HttpServletResponse response) throws Exception 
	{
		sessionUser = User.getNewInstance();
		sessionUser.suspended = true;
		sessionUser.referralId = getReferral(request);

		WebModel wm = new WebModel(sessionUser);

		if (OCHEMConfiguration.inhouseInstallation)
			wm.addParam("inhouse", "true");
		if (OCHEMConfiguration.disableAnonymousUsers)
			wm.addParam("noanonymous", "true");
		ModelAndView mav = new ModelAndView("xslt");
		mav.addObject("object", wm);
		return mav;
	}

	void init(HttpServletRequest req) 
	{
		Session session = Globals.userSession();
		if (session != null)
		{
			sessionUser = session.user;
		}
		if (sessionUser == null)
			sessionUser = User.getNewInstance();
	}

	@Override
	public ModelAndView list(HttpServletRequest req, HttpServletResponse res) throws Exception 
	{
		AccessChecker.requestSuperuserPrivileges();
		WebList list = new WebList();
		WebFilters filters = formFilters(req);

		Criteria criteria = Globals.session().createCriteria(User.getCurrentClass());
		if (assertParam("login"))
			criteria.add(Restrictions.like("login", "%"+getParam("login")+"%"));
		if(assertParam("username"))
			criteria.add(Restrictions.or(Restrictions.like("login", "%"+getParam("username")+"%"), Restrictions.or(Restrictions.like("firstName", "%"+getParam("username")+"%") ,Restrictions.like("lastName", "%"+getParam("username")+"%"))));
		if (assertParam("group"))
			criteria.add(Restrictions.eq("group.id", getLongParam("group")));
		if (assertParam("organisation"))
		{
			if ("empty".equals(getParam("organisation")))
				criteria.add(Restrictions.isNull("organisation"));
			else
				criteria.add(Restrictions.like("organisation", getParam("organisation")));
		}

		//criteria.add(Restrictions.isNotNull(getParam("sort")));
		criteria.addOrder(Order.desc(getParam("sort")));

		//criteria.addOrder(Order.desc("rank"));
		list.useEntity(User.getCurrentClass()).loadDistinctFromCriteria(criteria, getPageNum(), getPageSize(10));

		return new BrowserModel().setFilters(filters).setObject(list).setTemplate("user-browser").getModelAndView();
	}

	private Long getReferral(HttpServletRequest request) {
		if (request.getSession().getAttribute("referral") != null)
			return User.getByString("" + request.getSession().getAttribute("referral")).id;
		else
			return null;
	}

	private void updateActivityData() {
		GregorianCalendar cal = new GregorianCalendar();
		cal.add(java.util.Calendar.MONTH, -1);
		Date monthAgo = cal.getTime();
		Globals.session().createSQLQuery("update User set activities_count_total=(select count(*) from UserEvent natural join Session s where s.user_id=User.user_id)").executeUpdate();
		Globals.session().createSQLQuery("update User set activities_count_month=(select count(*) from UserEvent ue natural join Session s where s.user_id=User.user_id and ue.time >= :monthAgo)")
		.setParameter("monthAgo", monthAgo)
		.executeUpdate();

		Globals.session().createSQLQuery("update User set latest_activity_time=(select max(time) from UserEvent natural join Session s where s.user_id=User.user_id)").executeUpdate();
	}

	public ModelAndView all(HttpServletRequest request, HttpServletResponse response) throws Exception 
	{
		AccessChecker.requestSuperuserPrivileges();
		updateActivityData();
		return new WebModel().setTemplate("user-browser").getModelAndView();
	}

	@SuppressWarnings("unchecked")
	public ModelAndView action(HttpServletRequest request, HttpServletResponse response) throws Exception 
	{
		User user = null;
		if(assertParam("id"))
			user = (User) Globals.session().get(User.getCurrentClass(), getLongParam("id"));

		DynaWrap extended = null;
		if (user != null && user.isExtended()) {
			extended = user.dynaWrapped();
		}
		

		if(assertParam("action"))
		{
			if ("suspend".equals(getParam("action")) && user != null)
			{
				AccessChecker.requestSuperuserPrivileges();
				user.suspended = true;
			}
			else if ("unsuspend".equals(getParam("action")) && user != null)
			{
				AccessChecker.requestSuperuserPrivileges();
				user.suspended = false;
			}
			else if ("changerank".equals(getParam("action")) && user != null)
			{
				AccessChecker.requestSuperuserPrivileges();
				Integer newRank = getIntParam("rank");
				if (extended  != null  && user.rank < User.VALIDATED && newRank == User.VALIDATED )
					Mailer.postMailAsync(new Email(extended.getString("email"), "Your OCHEM account has been validated!", "Dear " + extended.getString("firstName") + " " + extended.getString("lastName") + ",\n\nWe are happy to inform you that your OCHEM account \"" + user.login + "\" has been validated by the administration team!\nWe wish you a lot of fun. Thank you for trying OCHEM.\n\nBest regards,\nYour OCHEM team"));
				user.rank = newRank;
			}
			else if ("delete".equals(getParam("action")) && user != null)
			{
				AccessChecker.requestSuperuserPrivileges();
				List<Session> sessions = Globals.session().createCriteria(Session.class).add(Restrictions.eq("user", user)).list();
				for (Session session : sessions) {
					session.user = null;
				}
				Globals.session().delete(user);
			}
			else if(getParam("action").equals("check"))
			{
				String theuser = request.getParameter("login"); if(theuser != null) theuser = theuser.trim();
				User tuser = User.getByLogin(theuser);
				if(tuser == null)
					throw new UserFriendlyException("user "+theuser == null?"":theuser+" not found");
			}
			else{
				String login = "";
				String email = "";
				String passwd = "";

				if(getParam("action").equals("submit"))
				{
					login = getParamTrim("login");
					email = getParamTrim("emailId");
					passwd = getParam("passwd");

					user = User.getByLogin(login);
					if (user != null)
						throw new UserFriendlyException("This login name is already exit in our database. Try another login name");

					user = User.getByEmail(email);
					if (user != null)
						throw new UserFriendlyException("This E-mail is already exit in our database.");

					user = User.getNewInstance();
					if (!isAllowedLogin(login))
						throw new UserFriendlyException("Not allowed login " + login);
					user.login = login;
					user.referralId = getReferral(request);
					
					if  (user.isExtended())  {
						extended = user.dynaWrapped();
						extended.setField("email", email);
						User.getExtended().getMethod("setPassword", String.class).invoke(user, passwd);
					}

					if (user.referralId != null) {
						// FIXME midnighter: Add a bonus for inviting a new user
					}
				}
				else if(getParam("action").equals("update"))
				{
					email = request.getParameter("emailId"); if(email!=null)email = email.trim();
					passwd = request.getParameter("passwd");
					user = User.getByEmail(email);
					user.message = "your profile has been successfully updated";
					if(extended != null && passwd != null && !extended.getString("password").equals(passwd))
					{
						User.getExtended().getMethod("setPassword", String.class).invoke(user, passwd);
						extended.setField("message", "Password has been sucessfully updated");
					}
				}

				if (extended != null) {
					if (assertParam("title")) {
						extended.setField("title", request.getParameter("title"));
						if(extended.getString("title") == null || extended.getString("title").length() == 0 || extended.getString("title").length() >= 8)
							throw new UserFriendlyException("Provide your Title.");
					}
					if (assertParam("firstname"))
						extended.setField("firstName", request.getParameter("firstname"));
					if (assertParam("lastname"))
						extended.setField("lastName", request.getParameter("lastname"));
					if (assertParam("organisation"))
						extended.setField("organisation", request.getParameter("organisation"));

					if (assertParam("city"))
						extended.setField("city", request.getParameter("city"));
					if (assertParam("state"))
						extended.setField("state", request.getParameter("state"));
					if (assertParam("country"))
						extended.setField("country", request.getParameter("country"));
					if (assertParam("zip"))
						extended.setField("zip", request.getParameter("zip"));
					if (assertParam("telephone"))
						extended.setField("telephone", request.getParameter("telephone"));

					if (assertParam("affiliation"))
						extended.setField("affil", request.getParameter("affiliation"));		

					if (assertParam("company"))
						extended.setField("company", request.getParameter("company"));
					if (assertParam("occupation"))
						extended.setField("occupation", request.getParameter("occupation"));

					if (extended != null) {
						User.getExtended().getMethod("setLink", String.class).invoke(user, request.getParameter("website"));
					}
				}
				

				// new user
				if (user.id == null)
				{
					user.rank = User.NOTVALIDATED; // Eto loh.

					// Automatically validate users for inhouse installations
					if (OCHEMConfiguration.inhouseInstallation)
						user.rank = User.VALIDATED;

					if (assertParam("invite"))
					{
						Invite invite = Invite.getFree(request.getParameter("invite"));
						if (invite != null)
						{
							invite.user = user;
							user.suspended = false;
						}
					}
					
					if (extended != null) {
						extended.setField("registrationTime", Calendar.getInstance().getTime());
					}

					Globals.session().saveOrUpdate(user);
					loginUser(user);

					String message;
					if (extended != null) {
						message = "A new user, "+extended.getString("firstName") + " "+extended.getString("lastName") + " "+ 
								extended.getString("email") + ", has just registered!\n The lucky person's new login is " + user.login + ".\nHis IP address is " + request.getRemoteAddr() + "\n";
						if (extended.getString("affil") != null)
							message += "\nAffiliation: " + extended.getString("affil") + "\n";
						message += "You can view his profile at " + OCHEMConfiguration.getRootHost() + "/user/profile.do?login=" +  user.login;
					} else {
						message = "A new user, "+  user.login + ", has just registered!\n The lucky person's new login is " + user.login + ".\nHis IP address is " + request.getRemoteAddr() + "\n";
						message += "You can view his profile at " + OCHEMConfiguration.getRootHost() + "/user/profile.do?login=" +  user.login;
					}

					Mailer.notifyAdmins("New user has just registered", message);
				}
				else
				{
					Globals.session().merge(user);
				}
			}
		}
		return new WebModel().getModelAndView();
	}

	private boolean isAllowedLogin(String login)
	{
		String lower = login.toLowerCase().trim();
		if (lower.endsWith(".biz") || lower.endsWith(".com") || lower.endsWith(".info") || lower.endsWith(".net"))
			return false;
		return true;
	}

	public ModelAndView profile(HttpServletRequest request, HttpServletResponse response) throws Exception 
	{
		if (Globals.userSession() == null)
		{
			request.getSession().setAttribute("requested-page", getRequestURI(request));
			return redirect("login/show.do");
		}
		//Globals.setMarshallingOption(MarshallingOption.USER_RECORDS);
		Globals.setMarshallingOption(MarshallingOption.USER_PUBLIC_MODELS);
		User user = User.getByLogin(getParamTrim("login"));
		UserProfile profile = UserProfile.create(user);
		if (request.getParameter("first-login") != null && request.getParameter("first-login").equals("true")) {
			profile.firstLogin = true;
		}

		return new WebModel(profile).setTemplate("user-profile").getModelAndView();
	}

	public ModelAndView getSessionData(HttpServletRequest request,
			HttpServletResponse response)
	{
		return new WebModel().getModelAndView();
	}

	public ModelAndView contributions(HttpServletRequest request, HttpServletResponse response) throws Exception 
	{
		if (Globals.userSession() == null)
			return redirect("login/show.do");
		return new WebModel().setTemplate("user-contributions").getModelAndView();
	}


}
