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

import java.lang.reflect.Field;
import java.sql.Timestamp;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.List;

import javax.persistence.CascadeType;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.EnumType;
import javax.persistence.Enumerated;
import javax.persistence.FetchType;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.OneToMany;
import javax.persistence.Transient;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.hibernate.annotations.Cascade;
import org.hibernate.event.spi.PreUpdateEvent;

import qspr.Globals;
import qspr.annotations.Loggable;

@Entity
@XmlRootElement(name = "action")
public class Action 
{
	@Id
	@GeneratedValue
	@Column(name = "action_id")
	@XmlAttribute
	protected Integer id;

	@Column(name = "table_name")
	@XmlAttribute
	protected String entity;

	@Column(name = "primary_key")
	@XmlAttribute(name = "key")
	public Long primaryKey;

	@ManyToOne
	@JoinColumn(name = "session_id")
	public Session session;

	@Column(name = "comm")
	public String comment;

	@Column(name="action_time")
	@XmlTransient
	public Timestamp time;

	@XmlTransient@Transient
	private Object originalEntity;

	@XmlTransient
	@Column
	@Enumerated(EnumType.STRING)
	public ActionType type;


	@OneToMany(fetch = FetchType.LAZY, cascade = CascadeType.ALL, mappedBy = "action")
	@Cascade(org.hibernate.annotations.CascadeType.DELETE_ORPHAN)
	@XmlTransient
	public List<UserAction> useraction;


	@XmlAttribute
	public String getTime()
	{
		SimpleDateFormat dateFormat = new SimpleDateFormat("dd.MM.yyyy HH:mm");
		return dateFormat.format(time.getTime());
		//return time.toGMTString();
	}

	public Action(Object entity) throws Exception
	{
		this.originalEntity = entity;
		this.comment = "";
		@SuppressWarnings("rawtypes")
		Class eClass = entity.getClass();
		@SuppressWarnings("unchecked")
		Loggable classLoggable = (Loggable) eClass.getAnnotation(Loggable.class);
		if (classLoggable.name().equals(""))
		{
			String[] entNames = eClass.getName().split("\\."); 
			this.entity = entNames[entNames.length-1];
		}
		else
			this.entity = classLoggable.name();
		this.session = Globals.userSession();
		this.primaryKey = (Long)eClass.getField("id").get(entity);
		this.time = new Timestamp(Calendar.getInstance().getTimeInMillis()); 
	}

	public String getPropertyNameInLog(String property) throws Exception
	{
		Loggable propertyLoggable = (Loggable) this.originalEntity.getClass().getDeclaredField(property).getAnnotation(Loggable.class);
		Loggable classLoggable = (Loggable) this.originalEntity.getClass().getAnnotation(Loggable.class);

		if ((propertyLoggable != null && !propertyLoggable.exclude()) || (classLoggable.greedy() && propertyLoggable == null))
			return (propertyLoggable == null || propertyLoggable.name().equals("")) ? property : propertyLoggable.name();

		return null;
	}


	/** For the actions on ExperimentalProperty create UserAction entries to notify record introducer
	 * and current owner of the fact, that the record has changed                  
	 **/
	public void publish(PreUpdateEvent event) 
	{
		Object entity = event.getEntity();
		Object[] oldValues = event.getOldState();
		String[] properties = event.getPersister().getPropertyNames();
		User owner = null, introducer = null;

		if (entity instanceof ExperimentalProperty || entity instanceof Article || entity instanceof Property)
		{
			for (int i = 0; i < properties.length; i++)
			{
				if (properties[i].equals("owner"))
					owner = (User) oldValues[i];

				if (properties[i].equals("introducer"))
					introducer = (User) oldValues[i];
			}

			if (owner != null)
				UserAction.save(owner, this);

			if ((introducer != null) && (!introducer.equals(owner)))
				UserAction.save(introducer, this);

		}

		updateTime(entity);
		updateOwner(entity);
	}

	public void publish(Object entity) 
	{
		try
		{
			entity.getClass(); //TODO may not required

			if (entity instanceof ExperimentalProperty || entity instanceof Article || entity instanceof Property)
			{
				User owner = (User)entity.getClass().getField("owner").get(entity);
				User introducer = (User)entity.getClass().getField("introducer").get(entity);

				if (owner != null)
					UserAction.save(owner, this);

				if ((introducer != null) && (!introducer.equals(owner)))
					UserAction.save(introducer, this);
			}

			updateTime(entity);
			updateOwner(entity);

		} catch (Exception e)
		{
			e.printStackTrace();
		}

	}	

	/**
	 * Universal place to track creation and modification time of entity
	 * @param entity
	 */
	public static void updateTime(Object entity)
	{
		try
		{
			Field time = entity.getClass().getField("time");
			Field timeCreated = entity.getClass().getField("timeCreated");
			Timestamp currentTime = new Timestamp(Calendar.getInstance().getTimeInMillis());
			time.set(entity, currentTime);
			if (timeCreated.get(entity) == null)
				timeCreated.set(entity, currentTime);
		}
		catch (NoSuchFieldException e)
		{
			// Its ok.
		}
		catch (IllegalAccessException e)
		{

		}
	}

	public static void updateOwner(Object entity)
	{
		try
		{
			// Change owner of record
			Field owner = entity.getClass().getField("owner");
			if (Globals.userSession() != null && Globals.userSession().user != null && !Globals.userSession().user.isSuperUser()) // a superuser cannot become owner of the record
				owner.set(entity, Globals.userSession().user);

			// If introduced by anonymous, change introducer to this user / Midnighter
			Field introducer = entity.getClass().getField("introducer");
			if (introducer.get(entity) == null)
				if (Globals.userSession() != null && Globals.userSession().user != null && !Globals.userSession().user.isSuperUser())
					introducer.set(entity, Globals.userSession().user);
		}
		catch (NoSuchFieldException e)
		{
			// Its ok.
		}
		catch (IllegalAccessException e)
		{

		}
	}


	public Action()
	{

	}

	public static enum ActionType
	{
		CREATE, DELETE, UPDATE, APPLY
	}
}
