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

package com.eadmet.descriptorcache;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

class CustomSerializer
	{
		static byte[] floatArrayToBytes(float[] values, boolean[] skip)
		{
			if(values == null)return null;
			
			try 
			{
				ByteArrayOutputStream baos = new ByteArrayOutputStream();
				GZIPOutputStream gos = new GZIPOutputStream(baos);
				DataOutputStream dos = new DataOutputStream(gos);
				
				int length = 0;
				for (boolean b : skip)
					length += (b ? 0 : 1);
				
				dos.writeInt(length);
				
				for (int i=0; i<values.length; i++)
					if (!skip[i])
						dos.writeFloat(values[i]);
				
				dos.close();
				return baos.toByteArray();
			} catch (IOException e)
			{
				e.printStackTrace();
				return new byte[0];
			}
		}
		
		static float[] bytesToFloatArray(byte[] data)
		{
			if(data == null)return null;
			try 
			{
				ByteArrayInputStream bais = new ByteArrayInputStream(data);
				GZIPInputStream gis = new GZIPInputStream(bais);
				DataInputStream dis = new DataInputStream(gis);
				int length = dis.readInt();
				float[] result = new float[length];
				for (int i=0; i<length; i++)
					result[i] = dis.readFloat();
				dis.close();
				return result;
			} catch (IOException e)
			{
				e.printStackTrace();
				return new float[0];
			}
		}
		
		static byte[] stringArrayToBytes(String[] values, boolean[] skip) 
		{
			if(values == null)return null;
			try
			{
				ByteArrayOutputStream baos = new ByteArrayOutputStream();
				GZIPOutputStream gos = new GZIPOutputStream(baos);
				DataOutputStream dos = new DataOutputStream(gos);
				int length = 0;
				for (boolean b : skip)
					length += (b ? 0 : 1);
				dos.writeInt(length);
				for (int i=0; i<values.length; i++)
					if (!skip[i])
						dos.writeUTF(values[i]);
				dos.close();
				return baos.toByteArray();
			} catch (IOException e)
			{
				e.printStackTrace();
				return new byte[0];
			}
		}
		
		static String[] bytesToStringArray(byte[] data)
		{
			if(data == null)return null;
			try
			{
				ByteArrayInputStream bais = new ByteArrayInputStream(data);
				GZIPInputStream gis = new GZIPInputStream(bais);
				DataInputStream dis = new DataInputStream(gis);
				int length = dis.readInt();
				String[] result = new String[length];
				for (int i=0; i<length; i++)
					result[i] = dis.readUTF();
				dis.close();
				return result;
			} catch (IOException e)
			{
				e.printStackTrace();
				return new String[0];
			}
		}
	}