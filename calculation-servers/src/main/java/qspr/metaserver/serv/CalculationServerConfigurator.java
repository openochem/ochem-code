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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.net.URLConnection;
import java.util.HashSet;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Unmarshaller;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.dao.ChemInfEngine;
import qspr.dao.Various;
import qspr.metaserver.CalculationServer;
import qspr.metaserver.ServerPool;
import qspr.metaserver.protocol.Command;
import qspr.metaserver.protocol.ServerInfo;
import qspr.metaserver.transport.CSTransport;

import com.eadmet.utils.CryptUtils;
import com.eadmet.utils.OCHEMUtils;
import com.eadmet.utils.OSType;

public class CalculationServerConfigurator
{
	private static final Logger logger = LogManager.getLogger(CalculationServerConfigurator.class);
	private ServerInfo serverInfo;
	private PrintWriter out;
	public boolean silent = false;

	public final String LINUX ="linux",MAC="mac",WINDOWS="windows",AARCH64 ="aarch64";

	public CalculationServerConfigurator(ServerInfo serverInfo, PrintWriter out)
	{
		this.serverInfo = serverInfo;
		this.out = out;
	}

	public void configureServers(Applications applications, ServerRunnerConfiguration configuration, 
			String workingDirectory) throws InstantiationException,
	IllegalAccessException, ClassNotFoundException, IOException
	{
		boolean useClasspath = workingDirectory.startsWith("classpath:");
		// Load the license key

		String osName = OSType.getOsName();
		String platformType = LINUX;
		if(OSType.isMac())
			platformType = MAC;
		if (OSType.isWindows())
			platformType = WINDOWS;
		if(platformType.equals(LINUX) && OSType.isAarch64())
			platformType = AARCH64;

		serverInfo.platform = osName;
		serverInfo.failures = new HashSet<String>();

		for (ConfiguredApplication configuredApplication : configuration.applications)
		{
			if (configuredApplication.isDisabled())
				continue;
			Application application = applications.getApplication(configuredApplication.name);
			if (application == null)
			{
				out.println("WARNING: Unknown application " + configuredApplication.name + " in " + ServerRunner.VERSIONXML);
				serverInfo.addFailure(configuredApplication.name);
				continue;
			}

			CalculationServer calculationServer;

			try
			{
				calculationServer = (CalculationServer) Class.forName(application.className).newInstance();

				if (calculationServer == null)
					throw new Exception("Calculation server could not be resolved");

			} catch (Throwable e)
			{
				serverInfo.failures.add(application.name);
				out.println("Configuring " + application.name + " failed to resolve " + application.className);
				continue;
			}

			out.println("Configuring " + application.name + "..." + " succeeded");

			// Configure target file for the output from the calculation server
			String logFileId = calculationServer.supportedTaskType != null ? calculationServer.supportedTaskType.toLowerCase() : calculationServer.getClass()
					.getSimpleName();

			if (useClasspath)
				calculationServer.out = out;
			else
				calculationServer.out = new TunnelPrintWriter(new FileOutputStream(workingDirectory + "/output/" + logFileId + ".out", true), logFileId, true,
						silent);

			calculationServer.workingDirectory = workingDirectory;

			calculationServer.availableMemory = configuration.memoryLimit;

			calculationServer.gpuCard = configuration.gpuCard != null ? configuration.gpuCard: OCHEMUtils.findFile("nvidia-smi") != null? 0: CalculationServer.NO_GPU;

			// First global configuration params
			if (application.params != null)
				for (ApplicationParam param : application.params)
				{
					if (param.platform != null && !param.platform.equals(platformType))
						continue;
					param.value = param.value.replace("%TOOLS_DIR%", calculationServer.getToolsDirectory());
					calculationServer.setParam(param.name, param.value);
				}

			// Then local configuration params (provided in config For)
			if (configuration.getApplicationParams(application.name) != null)
			{
				for (ApplicationParam param : configuration.getApplicationParams(application.name))
				{
					if (param.platform != null && !param.platform.equals(platformType))
						continue;
					//					param.value = param.value.replace("%RELEASE_DIR%", workingDirectory);
					param.value = param.value.replace("%TOOLS_DIR%", calculationServer.getToolsDirectory());
					calculationServer.setParam(param.name, param.value);
				}
			}

			if (calculationServer.supportedTaskType == null)
			{
				out.println("Configuring " + application.name + " in" + application.className + " calculationServer.supportedTaskType == null");
				continue;
			}

			calculationServer.setRunsDirectory(workingDirectory + "/runs/");
			calculationServer.setDebugDirectory(workingDirectory + "/debugs/");



			calculationServer.javaHome = (configuration.javaHome == null) ? System.getProperty("java.home") : configuration.javaHome;
			if (new File(serverInfo.workingDirectory+"/lib").exists())
				calculationServer.javaClassPath = OSType.getPath(serverInfo.workingDirectory, "lib", "*");
			else
				calculationServer.javaClassPath = System.getProperty("java.class.path");

			ServerPool.getInstance().servers.add(calculationServer);
			configuredApplication.server = calculationServer;
		}
	}

	/**
	 * Unmarshal the configuration file from version.xml and version-template.xml
	 * @throws Exception 
	 */
	public static ServerRunnerConfiguration readConfigurationFile(JAXBContext context, String serverHomeDirectory) throws Exception
	{
		ServerRunnerConfiguration configuration;
		Unmarshaller unmarshaller = context.createUnmarshaller();
		File localConfigFile = new File(serverHomeDirectory + "/" + ServerRunner.VERSIONXML);
		if (localConfigFile.exists())
			configuration = (ServerRunnerConfiguration) unmarshaller.unmarshal(localConfigFile);
		else
			configuration = new ServerRunnerConfiguration();
		System.out.println(configuration.chemengine);

		configuration.makeApplicationsRedundant();

		downloadVersionTemplate(serverHomeDirectory, configuration.metaserverURL);
		File templateFile = new File(serverHomeDirectory + "/" + Command.VERSIONTEMPLATEXML);
		if (templateFile.exists())
		{
			logger.info("Merging configuration with " + Command.VERSIONTEMPLATEXML);
			ServerRunnerConfiguration templateConf = (ServerRunnerConfiguration) unmarshaller.unmarshal(templateFile);
			configuration.mergeWith(templateConf);
			Various.defaultEngine = configuration.chemengine == null? ChemInfEngine.CHEMAXON:configuration.chemengine;
		}

		return configuration;
	}

	/**
	 * Download an encoded file from metaserver
	 */
	private static void downloadVersionTemplate(String serverHomeDirectory, String metaserverURL)
	{
		try
		{
			logger.info("Downloading the " + Command.VERSIONTEMPLATEXML + " from " + metaserverURL);

			URLConnection openConnection = CSTransport.getConnection(metaserverURL + "/?action=getFile&file=" + Command.VERSIONTEMPLATEXML);

			openConnection.addRequestProperty("User-Agent", "OCHEM server"); // just to have this field not null

			BufferedReader r = new BufferedReader(new InputStreamReader(openConnection.getInputStream()));
			StringWriter writer = new StringWriter();
			String s;
			while ((s = r.readLine()) != null)
			{
				writer.write(s);
			}

			String fileXML = CryptUtils.desDecode(writer.toString());
			FileWriter fWriter = new FileWriter(new File(serverHomeDirectory + "/" + Command.VERSIONTEMPLATEXML));
			fWriter.write(fileXML.toString());
			fWriter.close();
		} catch (Exception e)
		{
			logger.error("Could not download " + Command.VERSIONTEMPLATEXML + " from metaserver");
			e.printStackTrace();
			System.exit(1); // no further calculations is possible if template cannot be downloaded
		}
	}
}
