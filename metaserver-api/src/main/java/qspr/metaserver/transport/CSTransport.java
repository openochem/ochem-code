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

package qspr.metaserver.transport;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;
import java.security.SecureRandom;
import java.security.cert.X509Certificate;
import java.util.Calendar;

import javax.net.ssl.HostnameVerifier;
import javax.net.ssl.HttpsURLConnection;
import javax.net.ssl.SSLContext;
import javax.net.ssl.SSLSession;
import javax.net.ssl.TrustManager;
import javax.net.ssl.X509TrustManager;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.metaserver.protocol.Command;
import qspr.metaserver.protocol.Task;

import com.eadmet.utils.config.ConfigurableClass;
import com.eadmet.utils.config.ConfigurableProperty;

@ConfigurableClass(name = "metaserver", comment = "Metaserver connectivity options")
public class CSTransport implements Transport 
{
	private static transient final Logger logger = LogManager.getLogger(CSTransport.class);

	@ConfigurableProperty(name = "transport_timeout")
	public static int TRANSPORT_TIMEOUT = 1000 * 60 * 1; // 10 seconds
	public static int MAX_DEEP_SLEEP = TRANSPORT_TIMEOUT*50; //500 seconds
	public static final String UNDEFINED ="undefined";

	@ConfigurableProperty(name = "default_url")
	public static String defaultServerURL = UNDEFINED;
	public String serverURL = addSeparator(defaultServerURL);
	public static LocalTransport localTransport = new LocalTransport().setEmbedded(true);
	public static boolean metaserverDown = false;
	public long lastSuccessfullConnect = Calendar.getInstance().getTimeInMillis();
	public int deepSleepTime = TRANSPORT_TIMEOUT;

	public Command executeCommand(Command command) throws MalformedURLException, IOException, ClassNotFoundException{
		return executeCommand(command, true);
	}

	public static String addSeparator(String s) {
		return s == null ? null :
			s.lastIndexOf('/') == s.length() - 1 ? s : s + "/" ;		
	}

	public static HttpURLConnection getConnection(String serverURL) throws IOException {

		URL url = new URL(serverURL);

		//FIX for SSL validation of self-signed certificates
		if (serverURL.startsWith("https:"))
		{
			CSTransport.disableCertificateValidation();
			return (HttpsURLConnection) url.openConnection();
		}
		else
		{
			return (HttpURLConnection) url.openConnection();
		}

	}

	private static void disableCertificateValidation()
	{
		//System.out.println("[WebSiteAvailabiltyTest] disabling certificates");

		// Create a trust manager that does not validate certificate chains
		TrustManager[] trustAllCerts = new TrustManager[] { 
				new X509TrustManager()
				{
					public X509Certificate[] getAcceptedIssuers()
					{
						return new X509Certificate[0];
					}

					public void checkClientTrusted(X509Certificate[] certs, String authType)
					{
					}

					public void checkServerTrusted(X509Certificate[] certs, String authType)
					{
					}
				} 
		};

		// Ignore differences between given hostname and certificate hostname
		HostnameVerifier hv = new HostnameVerifier()
		{
			public boolean verify(String hostname, SSLSession session)
			{
				return true;
			}
		};

		// Install the all-trusting trust manager
		try
		{
			SSLContext sc = SSLContext.getInstance("SSL");
			sc.init(null, trustAllCerts, new SecureRandom());
			HttpsURLConnection.setDefaultSSLSocketFactory(sc.getSocketFactory());
			HttpsURLConnection.setDefaultHostnameVerifier(hv);
		} catch (Exception e)
		{
		}
	}


	public Command executeCommand(Command command, boolean firstCalculateLocally) throws MalformedURLException, IOException, ClassNotFoundException
	{
		try
		{
			if(!firstCalculateLocally)throw new UnsupportedOperationException();
			// Try to execute this command locally
			return localTransport.executeCommand(command);
		}
		catch (UnsupportedOperationException e)
		{
			// Local execution was not possible, send a real request to the Metaserver

			HttpURLConnection conn= null;

			try
			{
				if (command.data instanceof Task)
				{
					Task task = (Task) command.data;
					task.createReferences();
				}

				conn = getConnection(serverURL);		

				conn.addRequestProperty("User-Agent", "Mozilla/4.76"); 		        // We need to set cookies as below.
				conn.setRequestProperty("Cookie", "foo=bar");

				conn.setReadTimeout(TRANSPORT_TIMEOUT); // 4 minutes
				conn.setConnectTimeout(TRANSPORT_TIMEOUT); // connection timeout
				conn.setRequestProperty("Connection", "close");
				conn.setDoOutput(true);

				ObjectOutputStream oos = new ObjectOutputStream(conn.getOutputStream());
				oos.writeObject(command);
				oos.flush();

				// Get the response
				ObjectInputStream ois = new ObjectInputStream(conn.getInputStream());
				Command response = (Command) ois.readObject();
				ois.close();
				oos.close();

				conn.disconnect();

				lastSuccessfullConnect = Calendar.getInstance().getTimeInMillis();
				metaserverDown = false;


				//if (response != null && response.data instanceof Task)
				//{
				//	((Task) response.data).resolveReferences();
				//}
				deepSleepTime = TRANSPORT_TIMEOUT; // we restore sleep to normal 10sec /// ????? Damn it...
				return response;
			}
			catch (IOException e2)
			{
				// we need to sleep until connection is closed automatically
				try
				{

					e2.printStackTrace();
					if (deepSleepTime > 0)
					{
						// FIXME Igor: Please, explain what this means
						long sec=(long)(2*deepSleepTime*Math.random());
						logger.info("The metaserver is unavailable; performing a DEEP SLEEP for "+sec/1000+" seconds to avoid" +
								" problem with metaserver overload and allow automatic closing of connections");
						Thread.sleep(sec);
					}
				}catch(Exception ee){ee.printStackTrace();}
				metaserverDown = true;
				throw e2;
			}finally
			{
				if (conn != null)
					conn.disconnect();
			}
		}
	}

	public CSTransport()
	{
		//Init URL
	}
}
