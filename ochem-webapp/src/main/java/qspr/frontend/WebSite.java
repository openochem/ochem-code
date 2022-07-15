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

package qspr.frontend;

import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.annotation.XmlAnyElement;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlValue;

import qspr.entities.Alert;

@XmlRootElement(name = "web-site")
public class WebSite 
{
	@XmlElement
	public String title;
	
	@XmlAttribute
	public boolean ajax;
	
	@XmlElement(name = "main-menu")
	MainMenu mainMenu;
	
	@XmlElement(name = "content")
	ContentPanel contentPanel;
	
	@XmlElement(name = "commons")
	Commons commons;
	
	public void addContent(Object object)
	{
		if (object != null)
			contentPanel.contents.add(object);
		else
			contentPanel.contents.add(new Alert("Trying to add empty content!"));
	}

	
	public WebSite()
	{
		mainMenu = new MainMenu();
		contentPanel = new ContentPanel();
	}
}

class Link
{
	@XmlAttribute
	String link;
	
	@XmlValue
	String value;
}

class ContentPanel
{
	@XmlAnyElement(lax = true)
	public List<Object> contents = new ArrayList<Object>();
}

class Commons
{
	@XmlAnyElement(lax = true)
	public List<Object> objects = new ArrayList<Object>();
}

class MainMenu
{
	@XmlElement(name = "item")
	List<Link> items = new ArrayList<Link>();
}



