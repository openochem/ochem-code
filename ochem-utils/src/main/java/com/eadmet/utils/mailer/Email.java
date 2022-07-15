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

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.io.IOUtils;

/**
 * Instead of passing a bunch of arguments to Mailer, pass a neat and extensible Email object.
 * @author midnighter
 *
 */
public class Email
{
	private String message;
	public String subject;
	public String recepients;
	public String attachmentPath;
	
	public String replyEmail;
	
	/**
	 * A flag to indicate that this email will be sent as HTML
	 */
	public boolean html;
	
	/**
	 * A classpath reference to the email template file (to be implemented)
	 */
	public String template;
	
	public Map<String, String> templateVars = new HashMap<String, String>();
	
	public Email(String recepients, String subject, String message) {
		this.recepients = recepients;
		this.subject = subject;
		this.message = message;
	}
	
	public Email setAttachmentPath(String attachmentPath) {
		this.attachmentPath = attachmentPath;
		
		return this;
	}
	
	public Email useHTML() {
		template = "email/basic.txt";
		html = true;
		return this;
	}
	
	public Email useTemplate(String template) {
		this.template = template;
		return this;
	}
	
	public String getMessage() {
		String result = null;
		if (html)
			result = message.replaceAll("\n", "<br/>");
		else
			result = message;
		
		if (template != null)
			try
			{
				result = loadTemplate(result);
			} catch (IOException e)
			{
				e.printStackTrace();
				return result;
			}
		
		return result;
	}
	
	public String loadTemplate(String body) throws IOException {
		String tpl = IOUtils.toString(this.getClass().getClassLoader().getResourceAsStream(template), "utf-8");
		
		tpl = tpl.replace("%body%", body);
		for (String var : templateVars.keySet())
			tpl.replaceAll("%" + var + "%", templateVars.get(var));
		return tpl;
	}
}
