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
import java.util.List;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.entities.User;
import qspr.util.DynaWrap;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.mailer.Email;
import com.eadmet.utils.mailer.Mailer;

@Controller
public class MailController extends ControllerWrapper {
	//	private static transient final Logger logger = Logger.getLogger(MailController.class);

	@SuppressWarnings("unchecked")
	public ModelAndView action(HttpServletRequest request,
			HttpServletResponse response) throws Exception {
		if ("reminder".equals(getParam("action")))
		{
			// Thats a password reminder!

			List<User> users;
			if (assertParam("login"))
				users = Globals.session().createQuery("from User where login=:login").setString("login", request.getParameter("login")).list();
			else if (assertParam("email"))
				users = Globals.session().createQuery("from User where email=:email").setString("email", request.getParameter("email")).list();
			else
				throw new UserFriendlyException("Please specify either email or login!");

			if (users.size() > 0 && (users.get(0).isExtended())) 
			{
				DynaWrap _user = users.get(0).dynaWrapped();
				String message = setMessage(null, _user.getString("lastName"), users.get(0).login, _user.getString("passwordEncoded"), null, null);
				Mailer.postMailSafely(new Email(_user.getString("email"), "Password reminder", message));
				ModelAndView mav = redirect("static/pwd-reminder-success.do");
				return mav;
			}
			else
				throw new UserFriendlyException("We are sorry, but we could not find a user with the data that you provided");
		}

		return null;
	}

	private String setMessage(String s_lastName, String r_lastName, String loginName, String pwd, String msg, String subj) {
		String message = "";
		if(pwd != null)
		{
			//message for password reminder
			message = "Dear "+r_lastName+"\n\nYou receive this email because you had requested a password be sent for your account at the OCHEM site. If you did not request this email then please ignore it, if you keep receiving it please contact the OCHEM Team.";
			message = message + "\n\n\nLogin Name\t-->\t"+loginName+"\nPassword\t-->\t"+pwd+"\n\nYou can change this password via the profile page. If you have any difficulties please contact the OCHEM Team.\n\n\nBest regards,\nYour OCHEM Team";
			return message;
		}
		else
		{
			//message between different user
			message = "Dear "+r_lastName+"\n\n\t"+s_lastName+" has sent you message\n\n" +
					"--------------------------------------------------------------------" +
					"\n"+"Dear "+r_lastName+",\n"+msg;
			message = message + "\nSincerely,\n"+s_lastName;
			message = message + "\n--------------------------------------------------------------------" +
					"\n\nTo reply this message, please visit your message box \n"+OCHEMConfiguration.getRootHost()+"/mail/show.do";
			message = message + "\n\n\nThank you\nOCHEM Team";
			return message;
		}
	}
}
