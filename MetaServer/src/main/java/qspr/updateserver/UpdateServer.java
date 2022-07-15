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

package qspr.updateserver;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.util.Calendar;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.cfg.Configuration;

import qspr.metaserver.AbstractServer;
import qspr.metaserver.ServerServlet;
import qspr.metaserver.protocol.Command;
import qspr.metaserver.protocol.Release;

import com.eadmet.utils.config.ConfigurableClass;
import com.eadmet.utils.config.ConfigurableProperty;

@ConfigurableClass(name = "updateserver", comment="Update server configs")
public class UpdateServer extends AbstractServer
{
	private static transient final Logger logger = LogManager.getLogger(UpdateServer.class);

	public static UpdateServer instance;
	public long lastCheck = 0;
	public Release currentRelease;
	public Guard downloadGuard = new Guard(3);

	@ConfigurableProperty(name = "release_dir")
	public static String releasesDir = "/home/common/Workflow/updateserver/releases";

	@ConfigurableProperty(name = "maindb")
	public static Map<String, String> mainDbOverrides = new HashMap<String,String>()
	{
		private static final long serialVersionUID = 1L;
		{
		}
	};

	public long lastReleaseFileSize = 0;

	public UpdateServer(String id) throws IOException
	{
		super(id);
		instance = this;
	}

	@Override
	protected Configuration configureDbConnection(Configuration hibernateConf) throws Exception
	{
		for (String key : mainDbOverrides.keySet())
			hibernateConf.setProperty(key, mainDbOverrides.get(key));
		return hibernateConf;
	}

	@Override
	public Command executeCommand(Command request) throws Exception
	{
		long currentTime = Calendar.getInstance().getTimeInMillis();

		// Experimental: If somebody is checking for new release, let him finish
		// Presumably this is the reason for high loads during updates
		synchronized (this)
		{
			if (currentTime - lastCheck > 2000)
			{
				File releaseFile = getReleaseFile();
				if (releaseFile != null && (currentRelease == null || !releaseFile.getName().equals(currentRelease.version)))
				{
					long curFileSize = getFileSize(releaseFile.getAbsolutePath());
					logger.info("New release! "+releaseFile.getName() + ", file size is " + curFileSize + " bytes");
					if (curFileSize == lastReleaseFileSize) // check whether file size has not changed (to prevent fetching files being copied at the moment)
					{
						logger.info("Fetching the release file..");
						currentRelease = new Release();
						currentRelease.version = releaseFile.getName();
						currentRelease.fileSize = curFileSize;
						//currentRelease.release = getBytesFromFile(releaseFile);
						logger.info("Release is ready.");
					}
					else
						logger.info("File size validity check in 2 seconds..");
					lastReleaseFileSize = curFileSize;
				}
				lastCheck = currentTime;
			}
		}

		Release clientRelease = (Release) request.data;
		if (currentRelease == null || currentRelease.version.equals(clientRelease.version))
			return null;

		logger.info("Notifying " + request.senderId + " that there is a new version " + currentRelease.version);

		return new Command(0, currentRelease);
	}

	public static long getFileSize(String filename) throws Exception 
	{
		FileInputStream fis = null;
		try {
			File me = new File(filename);
			fis = new FileInputStream(me);
			return fis.getChannel().size();
		} finally {
			fis.close();
		}
	}


	private File getReleaseFile() throws IOException
	{
		File releaseDir = new File(releasesDir);
		File[] files = releaseDir.listFiles(new ZIPFilter());
		long lastModified = 0;
		File releaseFile = null;

		for (File file : files) {
			if (file.lastModified() > lastModified)
			{
				lastModified = file.lastModified();
				releaseFile = file;
			}
		}

		if(releaseFile == null)throw new IOException("No release was found at " + releaseDir.getCanonicalPath());

		return releaseFile;
	}

	/* Get last release file
	 Intended for fast initial import of server runner
	   i.e "wget http://qspr.eu/metaserver/update"
	 */

	@Override
	protected void serveGETRequest(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		// Midnighter

		try
		{
			downloadGuard.enter();
			File releaseFile = getReleaseFile();
			res.setContentType("application/zip");
			res.setHeader("Content-Disposition", "attachment; filename="+releaseFile.getName());
			BufferedOutputStream out = new BufferedOutputStream(res.getOutputStream());
			writeBytesFromFile(releaseFile, out);
			out.flush();
			out.close();
			logger.info("New release has been sent to the client: " + ServerServlet.resolveRemoteAddr(req));
		} catch (GuardException e)
		{
			logger.info("Too many connections, please retry later");
			res.sendError(HttpServletResponse.SC_SERVICE_UNAVAILABLE,"Too many connections, please retry later");
			return;
		}
		finally
		{
			downloadGuard.exit();
		}
	}

	public static void writeBytesFromFile(File file, BufferedOutputStream out) throws IOException 
	{
		InputStream is = new FileInputStream(file);

		byte[] bytes = new byte[2048];

		int numRead = 0;
		while ((numRead = is.read(bytes)) >= 0) 
			out.write(bytes, 0, numRead);

		is.close();
	}

	public static void main(String[] args) throws Exception {
		logger.info(getFileSize("/Users/midnighter/ParStream_2.0_JBT_V2.ps"));
		logger.info(new File("/Users/midnighter/ParStream_2.0_JBT_V2.ps").length());
	}
}

class ZIPFilter implements FilenameFilter {
	public boolean accept(File dir, String name)
	{
		return (name.endsWith(".zip"));
	}
}

/**
 * A simple guard to restrict too many entries to a particular functionality
 * @author midnighter
 *
 */

class GuardException extends IOException
{
	private static final long serialVersionUID = 1L;
	public GuardException(String message)
	{
		super(message);
	}
}

class Guard
{
	public int maximumEntries = 5;
	public AtomicInteger entries = new AtomicInteger(0);

	public void enter() throws GuardException
	{
		entries.addAndGet(1);
		if (entries.get() >= maximumEntries)
			throw new GuardException("Guarded amount exceeded ("+maximumEntries+")");
	}

	public void exit()
	{
		entries.addAndGet(-1);
	}

	public Guard(int maximumEntries)
	{
		this.maximumEntries = maximumEntries;
	}
}
