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

package com.eadmet.utils;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.StringReader;
import java.net.URI;
import java.util.Deque;
import java.util.LinkedList;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.commons.compress.archivers.tar.TarArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import org.apache.commons.compress.archivers.tar.TarArchiveOutputStream;
import org.apache.commons.compress.utils.IOUtils;

/**
 * A set of helper methods to work with files
 * 
 * @author midnighter
 * 
 */
public class FileUtils {
	private final static int TEMP_DIR_ATTEMPTS = 10;

	/*
	 * @data -- data to be saved
	 * 
	 * @file -- file to which the data will be saved
	 */
	public static void saveStringToFile(String data, String file) throws IOException {
		File f = new File(file);
		if(f.exists())f.delete();
		FileWriter pw = new FileWriter(f);
		BufferedReader r = new BufferedReader(new StringReader(data));
		String strLine;
		while ((strLine = r.readLine()) != null)
			pw.write(strLine + OSType.endLine());
		pw.flush();
		pw.close();
	}

	/**
	 * Create a temporary directory
	 */
	public static File createTempDir() {
		File baseDir = new File(System.getProperty("java.io.tmpdir"));
		String baseName = System.currentTimeMillis() + "-";

		for (int counter = 0; counter < TEMP_DIR_ATTEMPTS; counter++) {
			File tempDir = new File(baseDir, baseName + counter);
			if (tempDir.mkdir()) {
				return tempDir;
			}
		}
		throw new IllegalStateException("Failed to create directory within " + TEMP_DIR_ATTEMPTS + " attempts (tried " + baseName + "0 to " + baseName + (TEMP_DIR_ATTEMPTS - 1) + ')');
	}

	public static byte[] getFileAsBytes(String file) throws IOException {
		File myFile = new File(file);

		byte[] result = new byte[(int) myFile.length()];

		int totalBytesRead = 0;
		BufferedInputStream input = new BufferedInputStream(new FileInputStream(file));
		while (totalBytesRead < result.length) {
			int bytesRemaining = result.length - totalBytesRead;
			// input.read() returns -1, 0, or more :
			int bytesRead = input.read(result, totalBytesRead, bytesRemaining);
			if (bytesRead > 0) {
				totalBytesRead = totalBytesRead + bytesRead;
			}
		}
		input.close();

		return result;
	}

	static public String createTarFile(String sourceDir){
		TarArchiveOutputStream tarOs = null;
		try {
			// Using input name to create output name
			FileOutputStream fos = new FileOutputStream(sourceDir +".tar");
			OutputStream gos = new BufferedOutputStream(fos);
			tarOs = new TarArchiveOutputStream(gos);
			addFilesToTar(sourceDir, "", tarOs);     
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}finally{
			try {
				tarOs.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		return sourceDir +".tar";
	}

	static public void unTarFile(File tarFile) throws IOException{
		FileInputStream fis = new FileInputStream(tarFile);
		TarArchiveInputStream tis = new TarArchiveInputStream(fis);
		TarArchiveEntry tarEntry = null;

		// tarIn is a TarArchiveInputStream
		while ((tarEntry = tis.getNextTarEntry()) != null) {
			File outputFile = new File(tarFile.getParent() + File.separator + tarEntry.getName());

			if(outputFile.equals(tarFile)) {
				tarFile.renameTo(new File(tarFile.getAbsolutePath()+".tar"));
			}

			if(tarEntry.isDirectory()){

				System.out.println("outputFile Directory ---- " 
						+ outputFile.getAbsolutePath());
				if(!outputFile.exists()){
					outputFile.mkdirs();
				}
			}else{
				//File outputFile = new File(destFile + File.separator + tarEntry.getName());
				System.out.println("outputFile File ---- " + outputFile.getAbsolutePath());
				outputFile.getParentFile().mkdirs();
				//outputFile.createNewFile();
				FileOutputStream fos = new FileOutputStream(outputFile); 
				IOUtils.copy(tis, fos);
				fos.close();
			}
		}
		tis.close();
	}


	static public void addFilesToTar(String filePath, String parent, TarArchiveOutputStream tarArchive) throws IOException {
		File file = new File(filePath);
		// Create entry name relative to parent file path 
		String entryName = parent + file.getName();
		// add tar ArchiveEntry
		tarArchive.putArchiveEntry(new TarArchiveEntry(file, entryName));
		if(file.isFile()){
			FileInputStream fis = new FileInputStream(file);
			BufferedInputStream bis = new BufferedInputStream(fis);
			IOUtils.copy(bis, tarArchive);
			tarArchive.closeArchiveEntry();
			bis.close();
		}else if(file.isDirectory()){
			tarArchive.closeArchiveEntry();
			for(File f : file.listFiles()){        
				addFilesToTar(f.getAbsolutePath(), entryName+File.separator, tarArchive);
			}
		}          
	}


	public static void saveBytesToFile(byte[] bytes, String file) throws IOException {
		File f = new File(file);
		if(f.exists() && f.isDirectory())org.apache.commons.io.FileUtils.deleteDirectory(f);
		OutputStream output = new BufferedOutputStream(new FileOutputStream(f));
		output.write(bytes, 0, bytes.length);
		output.close();

		try {
			unTarFile(f);
		}catch(Exception e) {
		}

	}

	public static String getStreamAsString(InputStream stream) throws IOException 
	{
		StringBuilder sb = new StringBuilder();
		BufferedReader br = new BufferedReader(new InputStreamReader(stream));
		String line = null;
		while ((line = br.readLine()) != null)
			sb.append(line+"\n");
		br.close();
		return sb.toString();
	}

	public static String getFileAsString(String file) throws IOException {
		BufferedReader input = new BufferedReader(new FileReader(file));
		StringBuffer sb = new StringBuffer();
		String line;
		while ((line = input.readLine()) != null)
			sb.append(line + '\n');
		input.close();
		return sb.toString();
	}

	@SuppressWarnings("resource")
	public static void zip(File directory, File zipfile) throws IOException {
		URI base = directory.toURI();
		Deque<File> queue = new LinkedList<File>();
		queue.push(directory);
		OutputStream out = new FileOutputStream(zipfile);
		Closeable res = out;
		try {
			ZipOutputStream zout = new ZipOutputStream(out);
			res = zout;
			while (!queue.isEmpty()) 
			{
				directory = queue.pop();
				for (File kid : directory.listFiles()) 
				{
					String name = base.relativize(kid.toURI()).getPath();
					if (kid.isDirectory()) 
					{
						queue.push(kid);
						name = name.endsWith("/") ? name : name + "/";
						zout.putNextEntry(new ZipEntry(name));
					} else 
					{
						zout.putNextEntry(new ZipEntry(name));
						copy(kid, zout);
						zout.closeEntry();
					}
				}
			}
		} finally {
			res.close();
		}
	}

	private static void copy(InputStream in, OutputStream out) throws IOException 
	{
		byte[] buffer = new byte[1024];
		while (true) {
			int readCount = in.read(buffer);
			if (readCount < 0) {
				break;
			}
			out.write(buffer, 0, readCount);
		}
	}

	private static void copy(File file, OutputStream out) throws IOException 
	{
		InputStream in = new FileInputStream(file);
		try {
			copy(in, out);
		} finally {
			in.close();
		}
	}

	public static String convertToASCII(String fileName) {
		fileName = fileName.replaceAll("[\\s;\\:]+", "_").toLowerCase();
		fileName = fileName.replaceAll("\u00f6", "o");
		fileName = fileName.replaceAll("\u00fc", "u");
		fileName = fileName.replaceAll("\u00e4", "a");
		fileName = fileName.replaceAll("[^\\x00-\\x7F]", ""); // and eventually deleting all non-ASCII
		return fileName;
	}


	public static void touch(File file)
	{
		long timestamp = System.currentTimeMillis();
		try
		{
			if (!file.exists())
				new FileOutputStream(file).close();
			file.setLastModified(timestamp);
		}
		catch (IOException e)
		{
		}
	}

	//	not in use? (Rob 22.07.13	
	//	private static void copy(InputStream in, File file) throws IOException 
	//	{
	//		OutputStream out = new FileOutputStream(file);
	//		try {
	//			copy(in, out);
	//		} finally {
	//			out.close();
	//		}
	//	}
}
