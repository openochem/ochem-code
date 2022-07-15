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

package qspr;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.servlet.http.HttpServletRequest;

import com.eadmet.utils.config.ConfigurableClass;
import com.eadmet.utils.config.ConfigurableProperty;

import qspr.dao.ChemInfEngine;
import qspr.workflow.utils.QSPRConstants;


/**
 * Some global configuration options loaded from an external file
 * @author midnighter
 */
@ConfigurableClass(name = "ochem", comment = "Global OCHEM configurations")
public class OCHEMConfiguration
{

	@ConfigurableProperty(name = "login", comment = "Configurable properties of the login interface.")
	public static Map<String, String> loginInfo = new HashMap<String,String>();

	@ConfigurableProperty(name = "allow_guest_account", comment = "Allow users to login as guest.")
	public static boolean guestLogin = true;

	@ConfigurableProperty(name = "allow_provider_login", comment = "Allow users to login instantly with LDAP or social accounts.")
	public static boolean providerLogin = false;

	@ConfigurableProperty(name = "allow_registration", comment = "Allow users to register and login with custom credentials.")
	public static boolean registerLogin = true;

	@ConfigurableProperty(name = "root_host")
	static public String rootHost = "http://localhost:8080";

	@ConfigurableProperty(name = "root_dir")
	static public String rootDir = "";

	@ConfigurableProperty(name = "cheminformatics_engine", comment = "Choose between 'CDK' or 'CHEMAXON'.")
	static public String cheminfEngine = null;
	
	@ConfigurableProperty(name = "user_entity", comment = "Choose implementation of the User class and database entity.")
	static public String userEntity = QSPRConstants.DEFAULT_USER;

	@ConfigurableProperty(name = "dynamic_host_detection", comment = "Should we detect the host dynamically based on the requested URL?")
	static public boolean dynamicHostDetection = false;

	@ConfigurableProperty(name = "delete_old_models")
	public static boolean deleteOldModels = true;

	@ConfigurableProperty(name = "allow_external_services")
	public static boolean allowExternalServices = false;

	@ConfigurableProperty(name = "hide_private_molecules")
	static public boolean hidePrivateMolecules = true;

	@ConfigurableProperty(name = "data_download_restrictions")
	static public boolean dataDownloadRestrictions = true;

	@ConfigurableProperty(name = "mirror")
	static public boolean mirror = false;

	@ConfigurableProperty(name = "inhouse_installation", comment = "Inhouse installation means no license agreement, no bonus points accounting and stuff")
	static public boolean inhouseInstallation = true;

	@ConfigurableProperty(name = "trusted_ips")
	static public String trustedIPAddressesString;

	@ConfigurableProperty(name = "testing", comment = "When enabled, the server will run all tests on startup")
	static public boolean testing = false;

	@ConfigurableProperty(name = "xemistry_indexing", comment = "When enabled, the server will index molecules using xemistry")
	static public boolean xemistryIndexing = true;

	@ConfigurableProperty(name = "mmpsimilarity_indexing", comment = "When enabled and MMPFrag is supported, the server will index molecules using MMP")
	static public boolean mmpSimilarityIndexingTask = false;

	@ConfigurableProperty(name = "atommap", comment = "Indicates IP of server for atommaping")
	static public String atommap = "localhost";

	@ConfigurableProperty(name = "solvent", comment = "Indicates ID of article with solvent")
	static public Long solvent = 132547l;

	@ConfigurableProperty(name = "auto_login_user", comment = "User that is automaticlaly logged. Normally, it is used only on the developer machines for testing purposes")
	static public String autoLoginUser = null;

	@ConfigurableProperty(name = "disable_anonymous_users", comment = "Guest users will not be allowed")
	static public boolean disableAnonymousUsers = false;

	@ConfigurableProperty(name = "disable_db", comment = "Disable database (will make most installation types not functional)")
	static public boolean disableDB = false;

	@ConfigurableProperty(name = "preferred_calculation_server", comment = "The default preferred calculation server used for all calculations")
	static public String defaultPreferredServer = null;

	@ConfigurableProperty(name = "verbose_mode", comment = "Control details of OCHEM logs: 0 - detailed output is disabled.")
	static public int verboseMode = 1;

	@ConfigurableProperty(name = "downloadable_logs", comment = "Absolute paths to downloadable logs, required to download logs from OCHEM interface")
	public static Map<String, String> downloadableLogs = new HashMap<String,String>()
	{
		private static final long serialVersionUID = 1L;
		{
		}
	};

	@ConfigurableProperty(name = "maindb")
	public static Map<String, String> mainDbOverrides = new HashMap<String,String>();

	@ConfigurableProperty(name = "fragmentdb")
	public static Map<String, String> fragmentDbOverrides = new HashMap<String,String>();

	@ConfigurableProperty(name = "custom_logs_path", comment = "A place whre to put custom user-generated logs")
	static public String customLogsPath = null;

	@ConfigurableProperty(name = "log_load", comment = "Log load on transaction start")
	static public boolean logLoad = false;

	public static String getRootHost()
	{
		if (dynamicHostDetection)
		{
			// Determine the host based on the request URL
			HttpServletRequest request = null;
			try {
				Object scope = Class.forName("qspr.ThreadScope").getMethod("get").invoke(null);
				request = (HttpServletRequest) scope.getClass().getMethod("getHttpServletRequest").invoke(scope);
			} catch (Exception exp) {
				System.err.println("Failed to get servlet request to detect host. Returning null...");
				exp.printStackTrace();
			}
			if (request != null)
				return getRootHost(request);
		}
		return rootHost;
	}

	/**
	 * Parse the host URL out of HTTP request
	 * @param req
	 * @return root host URL
	 */
	private static String getRootHost(HttpServletRequest req)
	{
		String scheme = req.getScheme();
		String serverName = req.getServerName();
		int serverPort = req.getServerPort(); 
		StringBuffer url =  new StringBuffer();
		url.append(scheme).append("://").append(serverName);

		if ((serverPort != 80) && (serverPort != 443)) {
			url.append(":").append(serverPort);
		}

		return url.toString();
	}

	public static void setRootHost(String newRootHost)
	{
		rootHost = newRootHost;
	}

	private static List<String> trustedIPs = new ArrayList<String>();

	public static boolean isTrustedAddress(String ipAddress) {
		if (!trustedIPs.isEmpty())
			for (String trustedIP : trustedIPs)
				if (ipAddress.contains(trustedIP.trim()))
					return true;

		if (trustedIPAddressesString != null) {
			String[] parts = trustedIPAddressesString.split(",");
			for (String part : parts)
			{
				part = part.trim();
				if (!part.isEmpty() && ipAddress.contains(part))
					return true;
			}
		}

		return false;
	}

	public static void addTrustedAddress(String ipAddress) {
		trustedIPs.add(ipAddress);
	}

	public static ChemInfEngine getCheminfEngine() {
		if(cheminfEngine != null)return ChemInfEngine.valueOf(cheminfEngine.toUpperCase()); // first value in cfg file 
		if(System.getenv(QSPRConstants.ENV) != null) // second is from environment
			cheminfEngine = System.getenv(QSPRConstants.ENV); // second value in environment
		else
			cheminfEngine = ""+ChemInfEngine.NONE; // required only for testing and in servers; will cause crash of OCHEM browsers

		return ChemInfEngine.valueOf(cheminfEngine);
	}


	public static void setCheminfEngine(ChemInfEngine choice) {
		cheminfEngine = ""+choice;
	}

	public static boolean chemaxonLicenceAvailable() {
		return System.getenv("OCHEM_CHEMAXON_LICENSE") != null && System.getenv("OCHEM_CHEMAXON_LICENSE").length() > 0;
	}

	public static Class<?> getUserClass() throws ClassNotFoundException {
		return Class.forName(userEntity);
	}

}
