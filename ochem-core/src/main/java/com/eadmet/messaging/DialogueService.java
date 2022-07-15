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

package com.eadmet.messaging;

import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.hibernate.Criteria;
import org.hibernate.criterion.Disjunction;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;
import org.springframework.stereotype.Service;

import com.eadmet.utils.MAILERConstants;

import qspr.Globals;
import qspr.ThreadScope;
import qspr.entities.Message;
import qspr.entities.Session;
import qspr.entities.User;
import qspr.util.WrapperThread;

@Service
@SuppressWarnings("unchecked")
public class DialogueService
{
	public List<Dialogue> getDialogues(User user)
	{
		Set<Integer> opponents = new HashSet<Integer>();
		List<Integer> users = Globals.session().createSQLQuery("select sender_id from Message where receiver_id=:user group by sender_id")
				.setParameter("user", user.id).list();

		if (!users.isEmpty())
			opponents.addAll(users);
		users = Globals.session().createSQLQuery("select receiver_id from Message where sender_id=:user group by receiver_id")
				.setParameter("user", user.id).list();
		if (!users.isEmpty())
			opponents.addAll(users);

		List<Dialogue> dialogues = new ArrayList<Dialogue>();
		for (Integer userID : opponents)
		{
			User opponent = User.getById(userID);
			Dialogue dialogue = getDialogue(user, opponent);
			dialogue.opponent = opponent;
			dialogue.messages = null; // Save traffic
			dialogue.unread = Globals.session().createQuery("from Message where isRead=0 and receiver=:user and sender=:opponent")
					.setParameter("user", user)
					.setParameter("opponent", opponent)
					.list().size() > 0;
					dialogues.add(dialogue);
		}

		// Sort by the date of the last message
		Collections.sort(dialogues, new Comparator<Dialogue>(){
			@Override
			public int compare(Dialogue arg0, Dialogue arg1)
			{
				return 1 - arg0.lastMessage.time.compareTo(arg1.lastMessage.time);
			}
		});

		return dialogues;
	}

	public Dialogue getDialogue(User u1, User u2) {
		Criteria c = Globals.session().createCriteria(Message.class);
		Disjunction disj = Restrictions.disjunction();
		disj.add(Restrictions.and(Restrictions.eq("sender", u1), Restrictions.eq("receiver", u2)));
		disj.add(Restrictions.and(Restrictions.eq("sender", u2), Restrictions.eq("receiver", u1)));

		c.add(disj);
		c.addOrder(Order.asc("time"));

		Dialogue dialogue = new Dialogue();
		dialogue.messages = c.list();
		if (!dialogue.messages.isEmpty())
			dialogue.lastMessage = dialogue.messages.get(dialogue.messages.size() - 1);

		return dialogue;
	}

	public long getUnreadMessagesCount(User user) {
		if (user == null || user.id == null || user.id < 0)
			return 0;
		Criteria criteria = Globals.session().createCriteria(Message.class)
				.add(Restrictions.eq("receiver", user)).add(Restrictions.eq("isRead", false));
		criteria.setProjection(Projections.rowCount());
		return ((Long)criteria.list().get(0));
	}

	public void sendMessage(User sender, User receiver, String subject, String messageText, boolean systemMessage) {
		Message message = new Message();
		message.sender = sender;
		message.subject = subject;
		message.receiver = receiver;
		message.body = message.text = messageText;
		message.time = new Timestamp(Calendar.getInstance().getTimeInMillis());
		message.isRead = false;
		Globals.session().save(message);
	}

	public static void main(String[] args)
	{
		new WrapperThread()
		{
			@Override
			public void wrapped() throws Exception
			{
				ThreadScope.get().userSession = Session.getFirstSession(MAILERConstants.ADMIN);
				DialogueService service = new DialogueService();
				service.getDialogues(Globals.userSession().user);
			}
		}.run();
	}
}

