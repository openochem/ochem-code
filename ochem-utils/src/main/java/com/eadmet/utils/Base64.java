/**
 * <p>Encodes and decodes to and from Base64 notation.</p>
 * <p>Homepage: <a href="http://iharder.net/base64">http://iharder.net/base64</a>.</p>
 * 
 * <p>Example:</p>
 * 
 * <code>String encoded = Base64.encode( myByteArray );</code>
 * <br />
 * <code>byte[] myByteArray = Base64.decode( encoded );</code>
 *
 * <p>The <tt>options</tt> parameter, which appears in a few places, is used to pass 
 * several pieces of information to the encoder. In the "higher level" methods such as 
 * encodeBytes( bytes, options ) the options parameter can be used to indicate such 
 * things as first gzipping the bytes before encoding them, not inserting linefeeds,
 * and encoding using the URL-safe and Ordered dialects.</p>
 *
 * <p>Note, according to <a href="http://www.faqs.org/rfcs/rfc3548.html">RFC3548</a>,
 * Section 2.1, implementations should not add line feeds unless explicitly told
 * to do so. I've got Base64 set to this behavior now, although earlier versions
 * broke lines by default.</p>
 *
 * <p>The constants defined in Base64 can be OR-ed together to combine options, so you 
 * might make a call like this:</p>
 *
 * <code>String encoded = Base64.encodeBytes( mybytes, Base64.GZIP | Base64.DO_BREAK_LINES );</code>
 * <p>to compress the data before encoding it and then making the output have newline characters.</p>
 * <p>Also...</p>
 * <code>String encoded = Base64.encodeBytes( crazyString.getBytes() );</code>
 *
 *
 *
 * <p>
 * Change Log:
 * </p>
 * <ul>
 *  <li>v2.3.6 - Fixed bug when breaking lines and the final byte of the encoded
 *   string ended in the last column; the buffer was not properly shrunk and
 *   contained an extra (null) byte that made it into the string.</li>
 *  <li>v2.3.5 - Fixed bug in {@link #encodeFromFile} where estimated buffer size
 *   was wrong for files of size 31, 34, and 37 bytes.</li>
 *  <li>v2.3.4 - Fixed bug when working with gzipped streams whereby flushing
 *   the Base64.OutputStream closed the Base64 encoding (by padding with equals
 *   signs) too soon. Also added an option to suppress the automatic decoding
 *   of gzipped streams. Also added experimental support for specifying a
 *   class loader when using the
 *   {@link #decodeToObject(java.lang.String, int, java.lang.ClassLoader)}
 *   method.</li>
 *  <li>v2.3.3 - Changed default char encoding to US-ASCII which reduces the internal Java
 *   footprint with its CharEncoders and so forth. Fixed some javadocs that were
 *   inconsistent. Removed imports and specified things like java.io.IOException
 *   explicitly inline.</li>
 *  <li>v2.3.2 - Reduced memory footprint! Finally refined the "guessing" of how big the
 *   final encoded data will be so that the code doesn't have to create two output
 *   arrays: an oversized initial one and then a final, exact-sized one. Big win
 *   when using the {@link #encodeBytesToBytes(byte[])} family of methods (and not
 *   using the gzip options which uses a different mechanism with streams and stuff).</li>
 *  <li>v2.3.1 - Added {@link #encodeBytesToBytes(byte[], int, int, int)} and some
 *   similar helper methods to be more efficient with memory by not returning a
 *   String but just a byte array.</li>
 *  <li>v2.3 - <strong>This is not a drop-in replacement!</strong> This is two years of comments
 *   and bug fixes queued up and finally executed. Thanks to everyone who sent
 *   me stuff, and I'm sorry I wasn't able to distribute your fixes to everyone else.
 *   Much bad coding was cleaned up including throwing exceptions where necessary 
 *   instead of returning null values or something similar. Here are some changes
 *   that may affect you:
 *   <ul>
 *    <li><em>Does not break lines, by default.</em> This is to keep in compliance with
 *      <a href="http://www.faqs.org/rfcs/rfc3548.html">RFC3548</a>.</li>
 *    <li><em>Throws exceptions instead of returning null values.</em> Because some operations
 *      (especially those that may permit the GZIP option) use IO streams, there
 *      is a possiblity of an java.io.IOException being thrown. After some discussion and
 *      thought, I've changed the behavior of the methods to throw java.io.IOExceptions
 *      rather than return null if ever there's an error. I think this is more
 *      appropriate, though it will require some changes to your code. Sorry,
 *      it should have been done this way to begin with.</li>
 *    <li><em>Removed all references to System.out, System.err, and the like.</em>
 *      Shame on me. All I can say is sorry they were ever there.</li>
 *    <li><em>Throws NullPointerExceptions and IllegalArgumentExceptions</em> as needed
 *      such as when passed arrays are null or offsets are invalid.</li>
 *    <li>Cleaned up as much javadoc as I could to avoid any javadoc warnings.
 *      This was especially annoying before for people who were thorough in their
 *      own projects and then had gobs of javadoc warnings on this file.</li>
 *   </ul>
 *  <li>v2.2.1 - Fixed bug using URL_SAFE and ORDERED encodings. Fixed bug
 *   when using very small files (~&lt; 40 bytes).</li>
 *  <li>v2.2 - Added some helper methods for encoding/decoding directly from
 *   one file to the next. Also added a main() method to support command line
 *   encoding/decoding from one file to the next. Also added these Base64 dialects:
 *   <ol>
 *   <li>The default is RFC3548 format.</li>
 *   <li>Calling Base64.setFormat(Base64.BASE64_FORMAT.URLSAFE_FORMAT) generates
 *   URL and file name friendly format as described in Section 4 of RFC3548.
 *   http://www.faqs.org/rfcs/rfc3548.html</li>
 *   <li>Calling Base64.setFormat(Base64.BASE64_FORMAT.ORDERED_FORMAT) generates
 *   URL and file name friendly format that preserves lexical ordering as described
 *   in http://www.faqs.org/qa/rfcc-1940.html</li>
 *   </ol>
 *   Special thanks to Jim Kellerman at <a href="http://www.powerset.com/">http://www.powerset.com/</a>
 *   for contributing the new Base64 dialects.
 *  </li>
 * 
 *  <li>v2.1 - Cleaned up javadoc comments and unused variables and methods. Added
 *   some convenience methods for reading and writing to and from files.</li>
 *  <li>v2.0.2 - Now specifies UTF-8 encoding in places where the code fails on systems
 *   with other encodings (like EBCDIC).</li>
 *  <li>v2.0.1 - Fixed an error when decoding a single byte, that is, when the
 *   encoded data was a single byte.</li>
 *  <li>v2.0 - I got rid of methods that used booleans to set options. 
 *   Now everything is more consolidated and cleaner. The code now detects
 *   when data that's being decoded is gzip-compressed and will decompress it
 *   automatically. Generally things are cleaner. You'll probably have to
 *   change some method calls that you were making to support the new
 *   options format (<tt>int</tt>s that you "OR" together).</li>
 *  <li>v1.5.1 - Fixed bug when decompressing and decoding to a             
 *   byte[] using <tt>decode( String s, boolean gzipCompressed )</tt>.      
 *   Added the ability to "suspend" encoding in the Output Stream so        
 *   you can turn on and off the encoding if you need to embed base64       
 *   data in an otherwise "normal" stream (like an XML file).</li>  
 *  <li>v1.5 - Output stream pases on flush() command but doesn't do anything itself.
 *      This helps when using GZIP streams.
 *      Added the ability to GZip-compress objects before encoding them.</li>
 *  <li>v1.4 - Added helper methods to read/write files.</li>
 *  <li>v1.3.6 - Fixed OutputStream.flush() so that 'position' is reset.</li>
 *  <li>v1.3.5 - Added flag to turn on and off line breaks. Fixed bug in input stream
 *      where last buffer being read, if not completely full, was not returned.</li>
 *  <li>v1.3.4 - Fixed when "improperly padded stream" error was thrown at the wrong time.</li>
 *  <li>v1.3.3 - Fixed I/O streams which were totally messed up.</li>
 * </ul>
 *
 * <p>
 * I am placing this code in the Public Domain. Do with it as you will.
 * This software comes with no guarantees or warranties but with
 * plenty of well-wishing instead!
 * Please visit <a href="http://iharder.net/base64">http://iharder.net/base64</a>
 * periodically to check for updates or to contribute improvements.
 * </p>
 *
 * @author Robert Harder
 * @author rob@iharder.net
 * @version 2.3.6
 */
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
public class Base64
{
	public static byte[] decode(String encoded)
	{
		return org.apache.commons.codec.binary.Base64.decodeBase64(encoded);
	}
	
	public static byte[] decode(byte[] encoded)
	{
		return org.apache.commons.codec.binary.Base64.decodeBase64(encoded);
	}
	
	public static String encodeBytes(byte[] encoded)
	{
		return org.apache.commons.codec.binary.Base64.encodeBase64String(encoded);
	}
} 