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

package qspr.metaserver.util;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import com.eadmet.utils.OSType;

import qspr.workflow.utils.QSPRConstants;

/**
 * Creates a new temporary directory in the 'runs' directory with aliases to all files in the passed directory.
 * 
 * The created directory and all its contents can be deleted with {@link destroy}.
 * The created directory can be emptied (all files except the aliases will be removed) using {@link purge}. 
 * This can be used to easily control temporary files during calculation.
 *
 * @author itetko
 * @author Midnighter 
 */
public class AliasedDirectory
{
	/**
	 * Creates a temporary directory and aliases (symbolic links) to files in a source directory.
	 * 
	 * @param sourceDir Source directory. Files in this directory will be linked to. 
	 * @param runDir Parent directory in which the temporary directory will be created.
	 * @param out For logging. Can be null if no logging is desired.
	 * @throws IOException
	 * @throws InterruptedException
	 */

	private PrintWriter out; // Logging.
	private File dirSource = null; // Name of directory with original files.
	private File dirTarget = null; // Name of created alias directory.

	/**
	 *  Current implementation of AliasedDirectory creation was very unclear (Matthias :)?)
	 * @param exefile -- executable file, if available
	 * @param server -- server for which this directory is created
	 * @throws IOException
	 * @throws InterruptedException
	 */

	public AliasedDirectory(String exefile, PrintWriter out, String targetDirectory) throws IOException, InterruptedException
	{
		this.out = out;
		dirTarget = new File(targetDirectory);

		if (exefile != null) {
			File f = new File(exefile);
			if (f.isDirectory())
				setSourceDir(f);
			else
				setSourceDir(f.getParentFile());
		}

		create();
	}

	/**
	 *  It is OS dependent 
	 *  On windows file is copies, on linux we create symbolic links
	 * @throws IOException
	 * @throws InterruptedException
	 */

	private void create() throws IOException, InterruptedException
	{
		if (out != null)
			out.println(String.format("Creating alias directory %s", dirTarget.getAbsolutePath()));
		if (dirTarget.exists())
			recursiveDelete(new File(dirTarget.getPath()));
		dirTarget.mkdir();

		if (dirSource == null)
			return;

		if (out != null)
			out.println(String.format("Creating %d aliases", dirSource.list().length));
		File sourceDir = getSource(dirSource);
		if(sourceDir != null) 
			alias(sourceDir,sourceDir.list());

		alias(dirSource,dirSource.list());
		// Add developer code directory
		File developer = developerDir(dirSource);
		if(developer.exists())
			alias(developer,developer.list());
	}

	File developerDir(File dir) {
		return  new File(dir.getAbsolutePath().replace(QSPRConstants.TOOLS, "/ext-tools/"));
	}
	
	private File getSource(File source) {
		String path = source.getAbsolutePath();
		if(path.contains(QSPRConstants.TOOLS)) {
			File f = new File(QSPRConstants.SOURCE,source.getName());
			if(f.exists()) return f;
		}
		return null;
	}

	public  String alias(File dirSource, String [] filesToLink) throws IOException {
		return alias(dirSource, filesToLink, dirTarget);
	}

	static public  String alias(File dirSource, String [] filesToLink, File dirTarget) throws IOException {
		try {
			final StringBuilder sb = new StringBuilder();
			File targetFile = null;
			for (final String fileName : filesToLink)
			{
				targetFile = new File(dirTarget, fileName);
				File source = new File(dirSource, fileName);
				sb.append(String.format("ln -sf %s %s\n", source.getCanonicalPath(),
						targetFile.getAbsolutePath()));
			}
			java.lang.Runtime.getRuntime().exec(new String[] { "bash", "-c", sb.toString() }).waitFor();
			return targetFile == null ? null: targetFile.toString();
		}catch(Exception e) {
			throw new IOException(e.getMessage());
		}
	}

	private void setSourceDir(File f) throws IOException
	{
		File ff = f;
		if (!f.isDirectory())
			ff=developerDir(f);
		if(!ff.isDirectory())
			ff = getSource(ff);
		if(!ff.isDirectory())
			throw new IOException("Internal error: Non found tools directory while creating the aliased directory: " + f.getAbsolutePath());
		dirSource = ff;
	}

	public String getAliasedPath()
	{
		return this.dirTarget.getAbsolutePath() + OSType.fileSeparator();
	}

	public String getAliasedPath(String file)
	{
		return this.dirTarget.getAbsolutePath() + OSType.fileSeparator() + new File(file).getName();
	}

	/** Returns the alias directory, i.e., the newly created one containing the aliases. */
	public File getTargetDirectory()
	{
		return dirTarget;
	}

	/** Deletes the created alias directory and all its content. 
	 * @throws InterruptedException */
	public void deleteAlias() throws IOException, InterruptedException
	{
		if (out != null)
			out.println(String.format("Deleting aliased directory %s", dirTarget));

		recursiveDelete(dirTarget);
	}

	// Deletes a file or directory.
	static public void recursiveDelete(File path) throws IOException, InterruptedException
	{
		if (!path.exists())
			throw new IllegalArgumentException("Internal error: Attempt to delete non-existing file");

		for (File f : path.listFiles())
		{
			f.delete();
			if (f.isDirectory())
				recursiveDelete(f);
			f.delete();
		}

		if (path.exists())
			path.delete();
	}
}