
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
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.zip.CRC32;
import java.util.zip.DeflaterOutputStream;
import java.util.zip.InflaterInputStream;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import org.apache.commons.codec.binary.Hex;

/**
 * A set of miscellaneous utilities use by the whole the OCHEM infrastructure
 * @author different_people
 *
 */
public class OCHEMUtils 
{
	static public CRC32 crc32; 
	static public MessageDigest md5new;

	private static transient final Logger logger = LogManager.getLogger(OCHEMUtils.class);

	static {
		try {
			crc32 = new CRC32();
			md5new = MessageDigest.getInstance("MD5");
		} catch (Exception e)
		{
			e.printStackTrace();
		}
	}

	public static  String[] loadJars(String libDir) throws Exception{
		String[] jars=null;
		File dependencyDirectory = new File(libDir);
		if(!dependencyDirectory.isDirectory())dependencyDirectory = dependencyDirectory.getParentFile();
		if(!dependencyDirectory.exists() || !dependencyDirectory.isDirectory())return null;
		logger.info("Uploading jars from: " + dependencyDirectory.getAbsolutePath());
		File[] files = dependencyDirectory.listFiles();
		if(files != null){
			int i=0;
			for (File file : files) 
				if (file.getName().endsWith(".jar"))i++;
			jars = new String[i]; i=0;
			for (File file : files) 
				if (file.getName().endsWith(".jar")){ 
					OCHEMUtils.loadLibrary(file);
					jars[i++]=file.getCanonicalPath();
				}
		}
		return jars;
	}

	public static synchronized void loadLibrary(java.io.File jar) throws Exception
	{
		try {
			/*We are using reflection here to circumvent encapsulation; addURL is not public*/
			java.net.URLClassLoader loader = (java.net.URLClassLoader)ClassLoader.getSystemClassLoader();
			java.net.URL url = jar.toURI().toURL();
			/*Disallow if already loaded*/
			for (java.net.URL it : java.util.Arrays.asList(loader.getURLs())){
				if (it.equals(url)){
					return;
				}
			}
			java.lang.reflect.Method method = java.net.URLClassLoader.class.getDeclaredMethod("addURL", new Class[]{java.net.URL.class});
			method.setAccessible(true); /*promote the method to public access*/
			method.invoke(loader, new Object[]{url});
		} catch (final java.lang.NoSuchMethodException | 
				java.lang.IllegalAccessException | 
				java.net.MalformedURLException | 
				java.lang.reflect.InvocationTargetException e){
			throw new Exception(e);
		}
	}


	public static int fibdNullPosition(String commands[]) {
		int found = -1;
		for(int i=0 ; i < commands.length && found == -1; i++)
			if(commands[i] == null) found =i;
		return found;
	}

	/**
	 * Wrap an exception in a RuntimeException (only if necessary) and rethrow it.
	 */
	public static void rethrowSafely(Exception e) {
		if (e instanceof RuntimeException)
			throw (RuntimeException) e;
		else
			throw new RuntimeException(e);
	}

	public static byte[] MySqlCompatibleCompress(String data) throws IOException{
		return data == null? null : MySqlCompatibleCompress(data.getBytes());
	}

	private static byte[] MySqlCompatibleCompress(byte[] data) throws IOException
	{
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BufferedOutputStream bos = new BufferedOutputStream(baos);
		DeflaterOutputStream dos = new DeflaterOutputStream(bos);
		dos.write(data);
		dos.close();
		byte[] cdata = baos.toByteArray();  //compressed file 
		byte[] rdata = new byte[cdata.length + 4]; //For MySQL compatible "compress/uncompress" functions to work we need a 4 byte "length" header
		byte[] datal = ByteBuffer.allocate(Integer.SIZE / Byte.SIZE).order(ByteOrder.LITTLE_ENDIAN).putInt(data.length).array();
		System.arraycopy(datal, 0, rdata, 0, 4); // store size of data
		System.arraycopy(cdata, 0, rdata, 4, cdata.length); // store all data
		return rdata;
	}

	public static byte[] MySqlCompatibleUncompress(byte[] data) throws IOException
	{

		if(data == null) return null;

		ByteArrayInputStream bais = new ByteArrayInputStream(data, 4, data.length - 4); //Skip the 4-byte header
		BufferedInputStream bis = new BufferedInputStream(bais);
		InflaterInputStream iis = new InflaterInputStream(bis);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		byte[] buffer = new byte[1024];
		int count = -1;
		while ((count = iis.read(buffer)) > -1)
			baos.write(buffer, 0, count);
		iis.close();
		return baos.toByteArray();
	}


	static public long getCrc32(String str)
	{
		crc32.reset();
		for (int i = 0; i < str.length(); i++)
			crc32.update(str.codePointAt(i));
		return crc32.getValue();
	}

	static public String getMD5(String data)
	{
		//Below - the correct way to do it... but hey, we need reverse-compatibility
		//return getMD5(data.getBytes());

		byte[] bytedata = new byte[data.length()];
		for (int i=0; i<data.length(); i++)
			bytedata[i] = (byte)data.codePointAt(i);
		return getMD5(bytedata);
	}

	static public String getMD5(byte[]... data)
	{
		byte[] result;
		synchronized(md5new)
		{
			md5new.reset();

			for (int i = 0; i < data.length; i++)
			{
				if (data[i] == null)
					md5new.update("null".getBytes());
				else
					for (int j = 0; j < data[i].length; j++)
						md5new.update(data[i][j]);
			}

			result = md5new.digest();

			return String.valueOf(Hex.encodeHex(result));
		}
	}

	public static String[] splitSDF(final String sdf) {
		String[] sdfs;
		sdfs = sdf.contains(",")?sdf.split(","): sdf.contains(";")?sdf.split(";"):
			sdf.contains("$$$$\n") ? sdf.split("\\$\\$\\$\\$\n") : sdf.split("\\$\\$\\$\\$");
		return sdfs;
	}

	public static byte[] getStreamAsBytes(InputStream is) throws IOException
	{
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		byte[] buffer = new byte[1024];
		int count = -1;
		while ((count = is.read(buffer)) > 0)
			baos.write(buffer, 0, count);
		baos.close();
		return baos.toByteArray();
	}

	static public String getMD5(InputStream is)
	{
		byte[] result;
		try 
		{
			synchronized(md5new)
			{
				md5new.reset();

				byte[] buffer = new byte[1024];
				int count = -1;

				while ((count = is.read(buffer)) > 0)
					md5new.update(buffer, 0, count);

				result = md5new.digest();

				return String.valueOf(Hex.encodeHex(result));
			}
		}
		catch (Exception e)
		{
			return null;
		}
	}



	public static String cleanString(String value)
	{
		if (value == null)
			return null;

		value = value.trim().replaceAll("^[ \\t\\x0B\\f\\r]+", "").replaceAll("[ \\t\\x0B\\f\\r]+$", "");//.replaceAll("\\p{C}", "");//Extended trim

		if (value.equals(""))
			return null;

		return value;
	}

	public static String exceptionToString(Throwable t)
	{
		StringWriter sw = new StringWriter();
		PrintWriter pw = new PrintWriter(sw);
		t.printStackTrace(pw);
		pw.close();
		return sw.toString();
	}

	/**
	 * Generates a random literal name with a specified num of alphabetic characters
	 * @return
	 */
	public static String getRandomName(int numOfCharacters) {
		String name = "";
		for (int i = 0; i < 5; i++)
			name = name + String.valueOf((char) (Math.round(Math.random() * 22) + 'a'));

		return name;
	}

	/**
	 * Is the object of a wrapper boxed type?
	 * @param clazz
	 * @return
	 */
	public static boolean isWrapperType(Class<?> clazz)
	{
		return WRAPPER_TYPES.contains(clazz);
	}

	private static final Set<Class<?>> WRAPPER_TYPES = getWrapperTypes();
	private static Set<Class<?>> getWrapperTypes()
	{
		Set<Class<?>> ret = new HashSet<Class<?>>();
		ret.add(Boolean.class);
		ret.add(Character.class);
		ret.add(Byte.class);
		ret.add(Short.class);
		ret.add(Integer.class);
		ret.add(Long.class);
		ret.add(Float.class);
		ret.add(Double.class);
		ret.add(Void.class);
		return ret;
	}

	public static String getSizeBytes(double size){

		if(size>0.95*1024*1024*1024) return "" + NumericalValueStandardizer.getSignificantDigits(size/(1024*1024*1024.)) + "Gb";
		if(size>0.95*1024*1024) return "" + NumericalValueStandardizer.getSignificantDigits(size/(1024*1024.)) + "Mb";
		if(size>0.95*1024) return "" + NumericalValueStandardizer.getSignificantDigits(size/1024.) + "Kb";
		return ""+size + " bytes";
	}

	public static String readFile(String fileLocation) throws IOException {
		return new String(Files.readAllBytes(Paths.get(fileLocation)));
	}

	public static String[] appendString(String[] arr, String element) {
		if(arr == null) {
			arr = new String[1];
			arr[0] = element;
		}else 
			arr=append(arr,element);
		return arr;
	}

	public static <T> T[] append(T[] arr, T element) {
		final int N = arr.length;
		arr = Arrays.copyOf(arr, N + 1);
		arr[N] = element;
		return arr;
	}

	public static <T> T[] insert(T element, T[] arr) {
		final int N = arr.length;
		T[] arr1 = Arrays.copyOf(arr, N + 1);
		for(int i=0;i<arr.length;i++)
			arr1[i+1]=arr[i];
		arr1[0] = element;
		return arr1;
	}

	public static List<String> getFileAsStringList(String file) throws IOException {
		BufferedReader input = new BufferedReader(new FileReader(file));
		List<String> lines = new ArrayList<String>();
		String line;
		while ((line = input.readLine()) != null)
			lines.add(line);
		input.close();
		return lines;
	}

	public static String findFile(String desiredProgram) {
		String commands[] = {OSType.isWindows() ? "where" : "which", desiredProgram};
		String path = executeCommands(commands);
		if(path != null)path = path.replaceAll("(\\r|\\n)", "");
		logger.info("Was looking for:  "+desiredProgram+" and found: " + path);
		return path == null ? null:Paths.get(path).toString();
	}

	public static String executeCommands(String... commands) {
		String parsed = "", line = null;

		Process proc = null;
		try {
			ProcessBuilder pb = new ProcessBuilder(commands);
			proc = pb.start();
			int errCode = proc.waitFor();
			if (errCode == 0) {
				try (BufferedReader reader = new BufferedReader(new InputStreamReader(proc.getInputStream()))) {
					while((line = reader.readLine()) != null) {
						parsed +=line+"\n";
					}
				}
			} else {
			}
		} catch (IOException | InterruptedException ex) {
		}
		finally {
			if(proc != null)proc.destroyForcibly();
		}

		return parsed.length() == 0 ? null:parsed;
	}

	public static void shuffleLines(String fileName, boolean keepFirst) throws IOException {
		List<String> list = getFileAsStringList(fileName);
		PrintWriter pw = new PrintWriter(new FileWriter(fileName));
		if(keepFirst) {
			pw.write(list.get(0)+"\n");
			list.remove(0);
		}
		Collections.shuffle(list);
		for(String line:list)
			pw.write(line+"\n");
		pw.close();
	}

	public static String getFilteredBasketName(String name) {
		String id = "Any-Latin; Latin-ASCII; [\\u0080-\\u7fff] remove";
		com.ibm.icu.text.ReplaceableString result = new com.ibm.icu.text.ReplaceableString(name);
		com.ibm.icu.text.Transliterator.getInstance(id).transliterate(result); name = result.toString();
		name=name.replaceAll("[^a-zA-Z0-9() ]","_").replaceAll("_+", "_").replaceAll("_$", "").trim();
		if(name.length()>100)return name.substring(0,100);
		return name;
	}

	public static void main(String[] args) throws Exception{
		System.out.println(getFilteredBasketName("тестт(% a)"));
	}

	public static String trim(String article) {
		return article.replaceAll ("^\\p{IsWhite_Space}+|\\p{IsWhite_Space}+$", "");
	}

	public static void duplicateLines(String fileName) throws IOException {
		List<String> list = getFileAsStringList(fileName);
		PrintWriter pw = new PrintWriter(new FileWriter(fileName));

		Collections.shuffle(list);
		for(String line:list)
			pw.write(line+"\n");

		Collections.shuffle(list);
		for(String line:list)
			pw.write(line+"\n");

		pw.close();
	}
}

