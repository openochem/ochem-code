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

package qspr.util;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.ObjectInputStream;
import java.io.Serializable;
import java.sql.Timestamp;
import java.util.BitSet;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.Serializer;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.esotericsoftware.kryo.serializers.CompatibleFieldSerializer;
import com.esotericsoftware.shaded.org.objenesis.strategy.InstantiatorStrategy;
import com.esotericsoftware.shaded.org.objenesis.strategy.SerializingInstantiatorStrategy;

import de.javakaffee.kryoserializers.BitSetSerializer;

public class ClassCompressor 
{
	private static transient final Logger logger = LogManager.getLogger(ClassCompressor.class);
	public static ClassLoader classLoader = null;
	
	private static Kryo getKryo()
	{
		Kryo k = new Kryo();
		k.setDefaultSerializer(CompatibleFieldSerializer.class);
		k.addDefaultSerializer(Timestamp.class, TimestampSerializer.class);
		InstantiatorStrategy is = new SerializingInstantiatorStrategy();
		k.setInstantiatorStrategy(is);
		k.register(BitSet.class, new BitSetSerializer());
		k.register(Timestamp.class, new TimestampSerializer());
		return k;
	}
	
	public static byte[] objectToByte(Serializable obj)
	{
		try 
		{
			if (obj == null)
				return null;
			try {
				obj.getClass().getDeclaredMethod("beforeSerialize").invoke(obj);
			} catch (NoSuchMethodException e) {
				// it ok.
			}
			
			ByteArrayOutputStream os = new ByteArrayOutputStream();
			GZIPOutputStream gos = new GZIPOutputStream(os);
			Output o = new Output(gos);
			Kryo k  = getKryo();
			k.writeClassAndObject(o, obj);
			o.flush();
			o.close();
			return os.toByteArray();
		}
		catch (Exception e)
		{
			throw new RuntimeException(e);
		}
	}
	
	public static Serializable byteToObject(byte[] bytes)
	{
		try 
		{
			return byteToObject(bytes, classLoader);
		}
		catch (Exception e)
		{
			throw new RuntimeException(e);
		}
	}
	
	private static Serializable byteToObject(byte[] bytes, ClassLoader loader)
	{
		try 
		{
			Serializable obj = byteToObjectKryo(bytes, loader);
			return obj;
		}
		catch (Exception e)
		{
			logger.debug(e);
			logger.info("Got exception, possibly old style serialized object is used: "+e.getMessage());
		}
		
		try 
		{
			Serializable obj = byteToObjectJava(bytes, loader);
			logger.info("Successfully deserialized old style java native serialized object");
			return obj;
		}
		catch (Exception e)
		{
			e.printStackTrace();
			throw new RuntimeException(e);
		}
	}
	
	private static Serializable byteToObjectJava(byte[] bytes, ClassLoader loader) throws Exception
	{
		if (bytes == null || bytes.length == 0)
			return null;
		ByteArrayInputStream bais = new ByteArrayInputStream(bytes);
		GZIPInputStream gis = new GZIPInputStream(bais);
		ObjectInputStream ois = (loader == null) ? new ObjectInputStream(gis) : new ObjectInputStreamWithClassLoader(gis, loader);
		Serializable obj = (Serializable) ois.readObject();
		ois.close();
		return obj;
	}
	
	private static Serializable byteToObjectKryo(byte[] bytes, ClassLoader loader) throws Exception
	{
		if (bytes == null || bytes.length == 0)
			return null;
		ByteArrayInputStream bais = new ByteArrayInputStream(bytes);
		GZIPInputStream gis = new GZIPInputStream(bais);
		Input i = new Input(gis);
		Kryo k = getKryo();
		if (loader != null)
			k.setClassLoader(loader);
		Serializable obj = (Serializable)k.readClassAndObject(i);
		i.close();
		return obj;
	}
	
	public static Serializable cloneObject(Serializable object)
	{
		// A naive but universal way to deep-copy an object // Midnighter on Nov 29, 2011 
		return byteToObject(objectToByte(object));
	}
	
	void main() 
	{
		
	}
}


class TimestampSerializer extends Serializer<Timestamp> 
{
    public void write(Kryo kryo, Output output, Timestamp object) {
      output.writeLong(object.getTime(), true);
    }

    public Timestamp read(Kryo kryo, Input input, Class<Timestamp> type) {
      return new Timestamp(input.readLong(true));
    }

    public Timestamp copy(Kryo kryo, Timestamp original) {
      return new Timestamp(original.getTime());
    }
}

