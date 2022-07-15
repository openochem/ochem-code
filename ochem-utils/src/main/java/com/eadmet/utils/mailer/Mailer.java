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

package com.eadmet.utils.mailer;

import java.util.Properties;

import javax.activation.DataHandler;
import javax.activation.DataSource;
import javax.activation.FileDataSource;
import javax.mail.Message;
import javax.mail.MessagingException;
import javax.mail.Multipart;
import javax.mail.Session;
import javax.mail.Transport;
import javax.mail.internet.InternetAddress;
import javax.mail.internet.MimeBodyPart;
import javax.mail.internet.MimeMessage;
import javax.mail.internet.MimeMultipart;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.eadmet.utils.MAILERConstants;
import com.eadmet.utils.OCHEMUtils;
import com.eadmet.utils.config.ConfigurableClass;
import com.eadmet.utils.config.ConfigurableProperty;

@ConfigurableClass(name = "mailer", comment = "Email settings")
public class Mailer
{

	public static String ProjectName = "OCHEM";

	@ConfigurableProperty(name = "server_signature")
	public static String serverSignature = "";

	@ConfigurableProperty(name = "enable")
	public static boolean enable = true;

	@ConfigurableProperty(name = "from_email")
	public static String from = "OCHEM Team <"+MAILERConstants.EMAIL_OCHEM+">";

	@ConfigurableProperty(name = "reply_to_email")
	public static String reply = MAILERConstants.EMAIL_OCHEM;

	@ConfigurableProperty(name = "mail_host", comment = "SMTP server address")
	public static String mailHost;
	// TODO: Authorization, SSL, TLS, Login/Pass, Port number and so on.

	@ConfigurableProperty(name = "admin_emails")
	public static String admins = MAILERConstants.EMAIL_ADMIN;

	@ConfigurableProperty(name = "developer_emails")
	public static String developers = MAILERConstants.EMAIL_ADMIN;

	/**
	 * Disable the mailer for this thread
	 */
	private static ThreadLocal<Boolean> disable = new ThreadLocal<Boolean>();

	public static void notifyAdmins(String subject, String message)
	{
		Mailer.postMailAsync(new Email(admins, subject, "Dear " + ProjectName + " administartor,\n\n" + ProjectName + " got some important news for you.\n" + message
				+ "\n\nSincerely yours,\n" + ProjectName + " server\n\nP.S.\n" + serverSignature));
	}

	public static void notifyDevelopers(Exception e, String callingClass)
	{
		Mailer.postMailAsync(new Email(developers, "An exception occured in " + callingClass,
				"Dear developer,\n\nregretfully, an exception occured:\n\n" + OCHEMUtils.exceptionToString(e) + "\n\nSincerely yours,\n" + ProjectName
				+ " server\n\n" + "P.S.\n" + serverSignature));
	}

	public static void notifyDevelopers(String subject, String message)
	{
		Mailer.postMailAsync(new Email(developers, subject, "Dear " + ProjectName + " developer,\n\n" + message + "\n\nSincerely yours,\n" + ProjectName + " server\n\n"
				+ "P.S.\n" + serverSignature));
	}

	public static void postMail(Email email) throws MessagingException
	{
		if (!enable)
			return;

		if (disable.get() != null && disable.get())
			return;

		if (email.recepients == null || email.recepients.trim().equals("") || mailHost == null)
			return;

		boolean debug = false;

		// Set the host smtp address
		Properties props = new Properties();
		props.put("mail.smtp.host", mailHost);

		// create some properties and get the default Session
		Session session = Session.getDefaultInstance(props, null);
		session.setDebug(debug);

		// create a message
		Message msg = new MimeMessage(session);

		// set from address
		msg.setFrom(new InternetAddress(from));
		// set reply address
		InternetAddress[] replyTo = new InternetAddress[1];

		if (email.replyEmail == null)
			email.replyEmail = email.replyEmail;

		replyTo[0] = new InternetAddress(reply);
		msg.setReplyTo(replyTo);
		// set to address
		InternetAddress[] addressTo = InternetAddress.parse(email.recepients);
		msg.addRecipients(Message.RecipientType.TO, addressTo);
		// set subject
		msg.setSubject(email.subject);

		if (email.attachmentPath == null && !email.html)
			msg.setText(email.getMessage());
		else
		{
			MimeBodyPart messageBodyPart = new MimeBodyPart();

			if (email.html)
				messageBodyPart.setContent(email.getMessage(), "text/html");
			else
				messageBodyPart.setText(email.getMessage());

			Multipart multipart = new MimeMultipart();
			multipart.addBodyPart(messageBodyPart);

			if (email.attachmentPath != null)
			{
				// Part two is attachment
				messageBodyPart = new MimeBodyPart();
				DataSource source = new FileDataSource(email.attachmentPath);
				messageBodyPart.setDataHandler(new DataHandler(source));
				messageBodyPart.setFileName(email.attachmentPath.substring(email.attachmentPath.lastIndexOf("/") + 1));
				multipart.addBodyPart(messageBodyPart);
			}

			if(mailHost != null && mailHost.length() > 0) msg.setContent(multipart);
		}

		logger.info( (mailHost != null && mailHost.length() > 0? "Sending" : "NOT Sending ") + " an email to " +
				email.recepients + " with subject \"" + email.subject + "\"");

		// Send the message
		Transport.send(msg);
	}

	public static void postMailSafely(Email email)
	{
		try
		{
			postMail(email);
		} catch (MessagingException e)
		{
			e.printStackTrace();
		}
	}

	public static void postMailAsync(final Email email)
	{
		new Thread()
		{
			public void run()
			{
				postMailSafely(email);
			}
		}.start();
	}

	public static void enableForThread(boolean enable)
	{
		disable.set(!enable);
	}

	public static void main(String[] args) throws MessagingException
	{
		postMail(new Email("test@gmail.com", "A test HTML message", "<head></head><body><h1>Nice!</h1>Very nice...</body>").useHTML());
	}

	private static transient final Logger logger = LoggerFactory.getLogger(Mailer.class);
}