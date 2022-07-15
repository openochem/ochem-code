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

import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.entities.User;
import qspr.frontend.WebModel;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.messaging.Dialogue;
import com.eadmet.messaging.DialogueService;

@Controller
public class DialogueController extends ControllerWrapper
{
	@Autowired
	DialogueService dialogueService;
	
	public DialogueController() {
		sessionRequired = true;
	}
	
	public ModelAndView dialogue(HttpServletRequest request,
			HttpServletResponse response) {
		if (Globals.isGuestUser())
			throw new UserFriendlyException("Guest users cannot use the messaging and chat. Please, register a free account to use messaging.");
		return new WebModel().setTemplate("messaging/dialogue").getModelAndView();
	}
	
	public ModelAndView messages(HttpServletRequest request,
			HttpServletResponse response) {
		User opponent = User.getByLogin(getParam("user"));
		Dialogue dialogue = dialogueService.getDialogue(Globals.myself(), opponent);
		return new WebModel(dialogue).getModelAndView();
	}
	
	public ModelAndView dialogues(HttpServletRequest request,
			HttpServletResponse response) {
		List<Dialogue> dialogues = dialogueService.getDialogues(Globals.myself());
		return new WebModel().setList(dialogues).getModelAndView();
	}
	
	public ModelAndView markRead(HttpServletRequest request,
			HttpServletResponse response) {
		Globals.session().createQuery("update Message set isRead=1 where sender=:sender and receiver=:myself")
			.setParameter("sender", User.getByString(getParam("user")))
			.setParameter("myself", Globals.myself())
			.executeUpdate();
		
		return new WebModel().getModelAndView();
	}
	
	public ModelAndView send(HttpServletRequest request,
			HttpServletResponse response) {
		
		dialogueService.sendMessage(Globals.userSession().user, User.getByLogin(getParam("receiver")), null, getParam("text"), false);
		
		return new WebModel().setTemplate("messaging/dialogue").getModelAndView();
	}
}
