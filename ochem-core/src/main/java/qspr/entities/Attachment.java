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

package qspr.entities;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.Serializable;
import java.io.StringReader;
import java.util.Calendar;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import javax.persistence.Column;
import javax.persistence.EnumType;
import javax.persistence.Enumerated;
import javax.persistence.Transient;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.metaserver.transport.DataReference;
import qspr.metaserver.transport.DataReferenceException;
import qspr.metaserver.transport.DataReferenceFactory;
import qspr.metaserver.transport.DataReferencer;
import qspr.util.ClassCompressor;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.OCHEMUtils;

@XmlRootElement(name = "pdf")
@SuppressWarnings({"rawtypes","unchecked"})
public class Attachment<T> 
{
	String ATTACHMENT ="Attachment";

	public enum AttachmentType
	{
		PDF(1), SERIALIZABLE(2), MARSHALABLE(3);

		private int flag;

		private AttachmentType(int flag)
		{
			this.flag = flag;
		}

		int getInt()
		{
			return flag;
		}

		static AttachmentType fromInt(int arg)
		{
			switch (arg)
			{
			case 1:
				return PDF;
			case 2:
				return SERIALIZABLE;
			case 3:
				return MARSHALABLE;
			default:
				throw new RuntimeException("Invalid argument for AttachmentType " + arg);
			}
		}
	}

	@XmlAttribute
	public Long id;

	@XmlTransient
	public Long getId() {
		return id;
	}

	public void setId(Long id) {
		this.id = id;
	}

	private String md5;
	
	/**
	 * Integer representation of AttachmentType
	 */
	public Integer flag;
	
	@Enumerated(EnumType.STRING)
	@Column
	public AttachmentSource source;


	@XmlTransient
	private byte[] data;

	@Transient
	static DataReferencer referencer;

	@XmlTransient
	private T cachedObject;

	@XmlTransient
	public byte[] getData()
	{
		if(md5 == null && data == null)return null;
		if(data != null)return data;

		try{
			if(referencer == null) referencer = DataReferenceFactory.createReferencer();
			DataReference ref = DataReferenceFactory.createReference(md5, ATTACHMENT);
			return referencer.getDataBytes(ref);
		}catch(DataReferenceException e){
			logger.info("Could not retrieve the data");
			throw new UserFriendlyException("could not retrieve data using DataReference " + md5);
		}

	}

	private boolean isLocalAttachment(){
		switch(source){
		case Article:
		case ModelConfTemplate:
		case ModelTemplate:
		case ReadyPendingTask:
		case UserSettings:
			return true;
		default:
			return false;
		}
	}

	private void setData(byte[] arg)
	{
		if(arg == null){ md5 = null; data = null; return;}
		try{
			if(referencer == null) referencer = DataReferenceFactory.createReferencer();
			DataReference ref = referencer.putDataBytes(arg, ATTACHMENT);
			md5 = ref.getReference();
			if(md5.length() > 32)throw new UserFriendlyException("Too long reference > 32 " + md5);
			if(data != null)
				data = arg;  // this is part of published model; we should save data directly in the database
		}catch(DataReferenceException e){
			data = arg;
		}

		if(isLocalAttachment()) // also storing a local copy
			data = arg;
	}

	/**
	 * "Publishes", i.e. stores data also locally
	 */

	public void publish(){
		data = getData();
	}

	public boolean isLocallyStored(){
		return data != null;
	}
	
	public long getDataLength(){
		if(data != null) return  data.length;

		if(md5 != null) 
			try{
				if(referencer == null) referencer = DataReferenceFactory.createReferencer();
				DataReference ref = DataReferenceFactory.resolveReference(md5, ATTACHMENT);
				Long size =  referencer.getDataSize(ref);
				return size != null ? size : data == null? 0: data.length;
			}catch(RuntimeException e){
			}
		return 0;
	}

	public static Attachment getAttachment(byte[] data, AttachmentSource source, AttachmentType type)
	{
		Attachment pdf = new Attachment(data, type, source);
		List files = Globals.session().createCriteria(Attachment.class).add(Restrictions.eq("md5", pdf.md5)).list();

		if (files.size() != 0)
			return (Attachment)files.get(0);

		Globals.session().saveOrUpdate(pdf);
		return pdf;
	}


	/**
	 * It is not recommended constructor
	 * Use the one with explicit constructor source
	 */

	public Attachment()
	{

	}

	/**
	 * Using default SERIALIZABLE type 
	 * 
	 * @param object
	 * @param source
	 */

	public Attachment(T object, AttachmentSource source)
	{
		this.source = source;
		setObject(object, AttachmentType.SERIALIZABLE);
	}

	public Attachment(T object, AttachmentType type, AttachmentSource source)
	{
		this.source = source;
		setObject(object, type);
	}

	@XmlTransient
	private T getSerializableObject()
	{	
		return cachedObject = (T)ClassCompressor.byteToObject(getData());
	}

	private void setSerializableObject(T object)
	{
		cachedObject = object;
		setData(ClassCompressor.objectToByte((Serializable)object));
	}

	@XmlTransient
	private T getXmlObject() throws JAXBException, IOException
	{
		Unmarshaller unmarshaller = Globals.jaxbContext.createUnmarshaller();
		try{
			cachedObject = (T)unmarshaller.unmarshal(new GZIPInputStream(new ByteArrayInputStream(getData()))); // new type of attachment
		}catch(IOException e){
			cachedObject = (T)unmarshaller.unmarshal(new ByteArrayInputStream(getData()));
		}
		//Unmarshall small object to "clean" unmarshaller before finalization, i.e. this is a FIX
		unmarshaller.unmarshal(new StringReader("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?><message title=\"\"><message></message><time>01-04-2014 15:35:33</time></message>"));
		return cachedObject;
	}

	private void setXmlObject(T object) throws Exception
	{
		cachedObject = object;
		Marshaller marshaller = Globals.jaxbContext.createMarshaller();
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		OutputStream os = new GZIPOutputStream(baos);
		marshaller.marshal(cachedObject, os);
		os.flush();
		os.close();
		setData(baos.toByteArray());
	}

	@XmlTransient
	public T getObject()
	{
		long time = Calendar.getInstance().getTimeInMillis();
		try
		{
			if (cachedObject == null)
			{
				switch (AttachmentType.fromInt(flag))
				{
				case SERIALIZABLE:
					getSerializableObject();
					break;
				case MARSHALABLE:
					getXmlObject();
					break;
				default:
					throw new Exception("Unsupported attachment type");
				}
				logger.info("Getting object " + cachedObject.getClass().getSimpleName() + 
						" from database: "+(Calendar.getInstance().getTimeInMillis() - time) + "ms., size " + OCHEMUtils.getSizeBytes(getDataLength()));
			}

			return cachedObject;
		} catch (Exception e)
		{
			e.printStackTrace();
			throw new RuntimeException(e);
		}
	}

	public void setObject(T object)
	{
		setObject(object, AttachmentType.SERIALIZABLE);
	}

	public void setObject(T object, AttachmentType type)
	{
		try
		{
			cachedObject = object;
			this.flag = type.getInt();
			if (!Globals.isStandalone)
			{
				switch (type)
				{
				case PDF: // already byte[] type
					setData((byte[])object);
					break;
				case SERIALIZABLE:
					setSerializableObject(object);
					break;
				case MARSHALABLE:
					setXmlObject(object);
					break;
				default:
					throw new Exception("Unsupported attachment type");
				}	
				if (cachedObject != null)
					logger.info("Setting object " + cachedObject.getClass().getSimpleName() + ", size " + OCHEMUtils.getSizeBytes(getDataLength()));
			}
		} catch (Exception e)
		{
			throw new RuntimeException(e);
		}
	}

	public void updateObject()
	{
		if (cachedObject != null)
			setObject(cachedObject, AttachmentType.fromInt(flag));
	}

	public String getAttachmentReference()
	{
		logger.debug("Reference is " + md5);
		return md5; 
	}

	public boolean equals(Attachment obj)
	{
		if (this == obj)
			return true;
		return md5 != null && md5.equals(obj.md5);
	}

	public Attachment<T> getCopy()
	{
		Attachment<T> copy = new Attachment<T>();
		copy.data = data;
		copy.flag = flag;
		copy.source = source;
		copy.md5 = md5;
		return copy;
	}

	private static Logger logger = LogManager.getLogger(Attachment.class);

}
