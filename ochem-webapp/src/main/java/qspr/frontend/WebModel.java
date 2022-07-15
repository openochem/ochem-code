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
import java.util.Set;

import javax.xml.bind.annotation.XmlAnyElement;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlElementWrapper;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlValue;

import org.springframework.web.servlet.ModelAndView;

import qspr.OCHEMConfiguration;
import qspr.dao.ChemInfEngine;
import qspr.entities.Announcement;
import qspr.entities.Log;
import qspr.entities.Session;
import qspr.metaserver.configurations.StandartizationOptions;

@XmlRootElement(name = "model")
@SuppressWarnings({"rawtypes","unchecked"})
public class WebModel 
{
	@XmlAnyElement
	public Object element;

	@XmlAttribute(name = "web-root")
	public String webRoot = OCHEMConfiguration.getRootHost() + OCHEMConfiguration.rootDir + "/";

	@XmlElement
	public String innerUrl;

	@XmlElement
	public String templateName;

	@XmlElement
	public Session session;

	@XmlElement(name = "version-info")
	public String versionInfo;

	@XmlElement(name = "announcement")
	public Announcement announcement;

	@XmlElement
	public Log log;

	@XmlElementWrapper(name = "others")
	@XmlAnyElement(lax = true)
	public List<Object> listedObjects;
	
	@XmlElementWrapper(name = "standardizerOptions")
	@XmlAnyElement(lax = true)
	public List<Object> standardizers;
	
	@XmlElement(name = "allowRegistration")
	public boolean allowRegistration = OCHEMConfiguration.registerLogin;

	@XmlAttribute(name = "render-mode")
	public String renderMode;

	@XmlAttribute
	public String context;

	/**
	 * Some UI is different for inhouse installations. 
	 * Thus, pass the inhouse flag to UI.
	 */
	@XmlAttribute
	private boolean inhouse = OCHEMConfiguration.inhouseInstallation;

	@XmlAttribute
	private boolean noanonymous = OCHEMConfiguration.disableAnonymousUsers;

	@XmlAttribute
	static private Boolean nomoreusers = null;


	@XmlElement(name = "param")
	List<Param> params;

	private String redirect;

	/**
	 * Allows to override the default "outer" template
	 */
	public transient String outerTemplate;

	public WebModel(Object object)
	{
		this();
		element = object;
		nomoreusers = false;
	}

	public WebModel setObject(Object object)
	{
		element = object;
		return this;
	}
	
	public WebModel setStandardizers(List object)
	{
		standardizers = object;
		return this;
	}

	public WebModel setList(List object)
	{
		listedObjects = object;
		return this;
	}

	public WebModel setList(Set set)
	{
		if (set != null)
		{
			listedObjects = new ArrayList();
			listedObjects.addAll(set);
		}
		//listedObjects = object;
		return this;
	}

	public WebModel setLog(Log log)
	{
		this.log = log;
		return this;
	}

	public WebModel setRenderMode(String object)
	{
		renderMode = object;
		return this;
	}

	public WebModel setTemplate(String tplName)
	{
		templateName = tplName;
		return this;
	}

	public WebModel addParam(String key, String value)
	{
		if (params == null)
			params = new ArrayList<Param>();

		Param p = new Param();
		p.key = key;
		p.value = value;

		params.add(p);
		return this;
	}

	public WebModel addObject(Object obj)
	{
		if (listedObjects == null)
			listedObjects = new ArrayList<Object>();
		listedObjects.add(obj);

		return this;
	}

	public WebModel addObjects(List objects)
	{
		if (listedObjects == null)
			listedObjects = new ArrayList<Object>();
		listedObjects.addAll(objects);

		return this;
	}

	static public WebModel fromList(List object)
	{
		WebModel wm = new WebModel(null);
		wm.listedObjects = object;
		return wm;
	}


	public ModelAndView getModelAndView()
	{
		if (redirect != null)
			if (redirect.startsWith("http"))
				return new ModelAndView("redirect:" + redirect);
			else
				return new ModelAndView("redirect:" + OCHEMConfiguration.getRootHost() + OCHEMConfiguration.rootDir + "/" + redirect);


		ModelAndView mav = new ModelAndView("xslt");
		mav.addObject("object", this);
		return mav;
	}


	public WebModel setContext(String _context)
	{
		context = _context;
		return this;
	}

	public WebModel()
	{
		setStandardizers(orderStandardizers(StandartizationOptions.getStandardizers()));
	}
	
	private List<ChemInfEngine> orderStandardizers(List<ChemInfEngine> stds) {
		String defaultStd = OCHEMConfiguration.getCheminfEngine().toString();
		int index = 0;
		for (ChemInfEngine std : stds) {
			if (std.toString() == defaultStd) {
				break;
			}
			index += 1;
		}
		ChemInfEngine theDefault = stds.get(index);
		List<ChemInfEngine> stds_new = new ArrayList<>();
		stds_new.add(theDefault);
		for (ChemInfEngine std : stds) {
			if (std != theDefault) {
				stds_new.add(std);
			}
		}
		return stds_new;
	}

	public WebModel(String innerUrl)
	{
		this();
		if(innerUrl != null && OCHEMConfiguration.rootHost.contains("https:") && !innerUrl.contains("https:") && !innerUrl.contains("metaserver")) // Fix for HTTPS
			innerUrl = innerUrl.replace("http:", "https:");
		this.innerUrl = innerUrl;
	}

	public static WebModel redirect(String url)
	{
		WebModel wm = new WebModel();
		wm.redirect = url;
		return wm;
	}
}

class Param
{
	@XmlAttribute
	String key;

	@XmlValue
	String value;
}
