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

package qspr.metaserver.serv;

import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlElementWrapper;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;
import javax.xml.bind.annotation.XmlValue;

import qspr.dao.ChemInfEngine;
import qspr.metaserver.ServerPool;
import qspr.metaserver.protocol.Command;

@XmlRootElement(name = "server-runner-config")
public class ServerRunnerConfiguration 
{
	@XmlTransient
	public static Applications allApplications;

	/**
	 * The configuration singleton available to all calculation servers
	 */
	public static transient ServerRunnerConfiguration instance = new ServerRunnerConfiguration();

	//entries with @XmlAttribute and not parsed from metaserver to servers; i.e. they are local 
	@XmlAttribute(name = "current-version")
	public String currentVersion;

	public String sid;

	public String metaserverURL;

	public String ochemURL;

	public String mongoDbURL;

	/**
	 * If true, will keep local versions of  metaserverURL, ochemURL, mongoDbURL
	 * Required for an advanced configuration through tunneling, re-mapping, etc.
	 */
	@XmlElement
	public Boolean localURL;

	@XmlElement(name = "java-home")
	public String javaHome;

	@XmlElement
	public Integer memoryLimit = 64;

	@XmlElement
	public Integer maximumRuns;

	/**
	 * Indicate directory in which external tools are located
	 */
	@XmlElement
	public String externalToolsDirectory;

	@XmlElement
	// Exit if this lifetime is exceeded, but wait for a running task to finish
	public Integer minimumLifetime; // in seconds

	@XmlElement
	// After running idle more than this value, the server will terminate (e.g., not to run idle LSF servers without necessity)
	public Integer maximumIdleTime; 

	@XmlElement
	public Integer minimumPriority;

	@XmlElement
	public Boolean generateTestReport;

	@XmlElement(name = "scheduled-priority")
	public List<ScheduledPriority> scheduledPriorities;

	@XmlElement
	public Integer sleepTime = 1; // in seconds

	@XmlElement
	public Integer sleepTimeWhileCalculating = 30; // in seconds

	@XmlElement
	public Integer processId;

	@XmlElement
	public Integer parentProcessId;

	@XmlElement
	public Integer gpuCard;

	// default standartizer, required to run tests
	@XmlElement
	public ChemInfEngine chemengine;

	@XmlElementWrapper(name = "supported-applications")
	@XmlElement(name = "application")
	public List<ConfiguredApplication> applications = new ArrayList<ConfiguredApplication>();

	@XmlElement(name = "config")
	public List<LocalAppConfig> applicationConfigurations;

	public boolean isApplicationEnabled(String name)
	{
		ConfiguredApplication application = getApplication(name);
		if (application == null)
			return false;
		if (application.isDisabled())
			return false;
		return true;
	}

	public ConfiguredApplication getApplication(String name)
	{
		for (ConfiguredApplication app : applications) 
		{
			if (app.name.equalsIgnoreCase(name))
				return app;
		}
		return null;
	}

	public boolean enableApplication(String name)
	{
		ConfiguredApplication application = getApplication(name);
		if (application != null)
		{
			if (application.isDisabled())
			{
				application.enable();
				return true;
			}
			return false;
		}

		if (!allApplications.exists(name))
			return false;

		applications.add(new ConfiguredApplication(name));
		return true;
	}

	public boolean disableApplication(String name)
	{
		ConfiguredApplication application = getApplication(name);
		if (application == null || application.isDisabled())
			return false;
		application.disable();
		return true;

	}

	public List<String> getSupportedTaskTypes()
	{
		List<String> taskTypes = new ArrayList<String>();
		for (ConfiguredApplication application : applications)
		{
			if (application.server != null && ServerPool.getInstance().servers.contains(application.server))
				if (application.isAllowedToRun())
					taskTypes.add(application.server.supportedTaskType);
		}

		return taskTypes;
	}

	public Integer getMinimumPriority()
	{
		Integer priority = null;

		if (scheduledPriorities != null)
			for (ScheduledPriority schPriority : scheduledPriorities) 
			{
				int currentHours = Calendar.getInstance().get(Calendar.HOUR_OF_DAY);
				if (schPriority.startTime != null && schPriority.stopTime != null){
					boolean applicable = false;

					if (schPriority.startTime < schPriority.stopTime)
						applicable = currentHours >= schPriority.startTime && currentHours < schPriority.stopTime;
						else
							applicable = currentHours >= schPriority.startTime || currentHours < schPriority.stopTime;

							if(applicable && Calendar.getInstance().get(Calendar.DAY_OF_WEEK) % 6 == 1)  // we have week-end
								applicable = schPriority.weekend != null && schPriority.weekend; // applicable if only true; normally not applicable

							if (applicable)
								priority = schPriority.priority;
				}
			}

		if (priority == null)
			return minimumPriority;
		else
			return priority;
	}

	/**
	 * Merge some parameters from another configuration.
	 * It is useful, when some configuration parameters (e.g., MongoDB URL) are common for all the servers.
	 */
	public void mergeWith(ServerRunnerConfiguration mergeConf)
	{
		if(!Command.FIXED.equals(currentVersion) && (localURL == null || !localURL) ){
			mongoDbURL = mergeConf.mongoDbURL;
			ochemURL = mergeConf.ochemURL;
			if(mergeConf.chemengine != null) chemengine = mergeConf.chemengine;
			chemengine = chemengine == null? ChemInfEngine.CHEMAXON:chemengine;
		}
		if (mergeConf.applicationConfigurations != null)
			for (LocalAppConfig appConfig : mergeConf.applicationConfigurations)
				addApplicationConfig(appConfig);
	}

	/**
	 * Local parameter overrides for an application
	 */
	public List<ApplicationParam> getApplicationParams(String appName)
	{
		LocalAppConfig config = getApplicationConfig(appName);
		if (config != null)
			return config.params;
		else
			return null;
	}

	private LocalAppConfig getApplicationConfig(String appName)
	{
		if (applicationConfigurations != null)
			for (LocalAppConfig localAppConfig : applicationConfigurations) {
				if (localAppConfig.forApplication.equalsIgnoreCase(appName))
					return localAppConfig;
			}

		return null;
	}

	private void addApplicationConfig(LocalAppConfig appConfig)
	{
		if (applicationConfigurations == null)
			applicationConfigurations = new ArrayList<ServerRunnerConfiguration.LocalAppConfig>();
		if (getApplicationConfig(appConfig.forApplication) != null)
			getApplicationConfig(appConfig.forApplication).params = appConfig.params;
		else
			applicationConfigurations.add(appConfig);
	}

	public void setApplicationParam(String app, String param, String value)
	{
		LocalAppConfig appConfig = getApplicationConfig(app);
		if (appConfig == null)
		{
			appConfig = new LocalAppConfig();
			appConfig.forApplication = app;
			addApplicationConfig(appConfig);
		}
		appConfig.setParam(param, value, null);
	}

	public static class LocalAppConfig
	{
		@XmlAttribute(name="for")
		public String forApplication;

		@XmlElement(name = "param")
		public List<ApplicationParam> params = new ArrayList<ApplicationParam>();

		public void setParam(String name, String value, String platform)
		{
			for (ApplicationParam param : params)
			{
				if (param.name.equalsIgnoreCase(name))
				{
					param.set(name, value, platform);
					return;
				}
			}

			params.add(new ApplicationParam(name, value, platform));
		}
	}

	public static class ScheduledPriority
	{
		@XmlAttribute(name = "start")
		public Integer startTime;

		@XmlAttribute(name = "stop")
		public Integer stopTime;

		@XmlAttribute
		public Boolean weekend;

		@XmlValue
		public int priority;
	}

	public void makeApplicationsRedundant() {
		List<ConfiguredApplication> appls = applications;
		applications = new ArrayList<ConfiguredApplication>();
		for(ConfiguredApplication app:appls)if(getApplication(app.name) == null)applications.add(app);
	}
}


