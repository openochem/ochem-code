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

import java.sql.Timestamp;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.UUID;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.EntityListeners;
import javax.persistence.FetchType;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.OneToMany;
import javax.persistence.Transient;
import javax.servlet.http.HttpSessionListener;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.hibernate.annotations.FilterDef;
import org.hibernate.annotations.FilterDefs;
import org.hibernate.annotations.ParamDef;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.metaserver.transport.Transport;
import qspr.metaserver.transport.TransportFactory;
import qspr.workflow.utils.QSPRConstants;
import qspr.util.AccessChecker;

@FilterDefs({
	@FilterDef(
			name="sessionFilter",
			parameters= {
					@ParamDef(
							name = "sessionId",
							type = "long"
							)	
			}
			),
	@FilterDef(
			name="userFilter",
			parameters= {
					@ParamDef(
							name = "userId",
							type = "long"
							)	
			}
			),
	@FilterDef(
			name="basketFilter",
			parameters= {
					@ParamDef(
							name = "basketId",
							type = "java.lang.Long"
							)	
			}
			)
})

@Entity	
@EntityListeners(HttpSessionListener.class)
@XmlRootElement(name="session")
public class Session 
{
	@Id
	@Column(name = "session_id")
	@GeneratedValue
	@XmlAttribute
	public Long id;

	@ManyToOne
	@JoinColumn(name="user_id")
	@XmlElement(name="user")
	public User user;

	@Column(name="ip_address")
	@XmlElement(name="ip-address")
	public String ipAddress;

	@Column(name="activity_time")
	@XmlTransient
	public Timestamp time;

	@Column
	public String session_string_id;

	@Column
	public Integer session_expired;

	@Column
	@XmlTransient
	public Timestamp mod_time;

	@Column
	//@XmlTransient
	public String guid;

	@Transient
	public int userId = 1;

	@Transient
	public int groupId = 1;

	@Transient
	public int groupMask = 1;

	@Column
	public boolean license;

	@Column(name="user_agent")
	public String userAgent;

	//@Transient
	//@XmlTransient
	//public Set<Long> selectionList;

	@Transient@XmlTransient
	public Set<Long> selectionMoleculeList;

	@Transient@XmlTransient
	public Set<Long> visitedModelsIds = new HashSet<Long>();

	@Transient@XmlTransient
	public Set<Long> selectedAlerts = new HashSet<Long>();

	@Transient
	static public Map<String,Map.Entry<String, byte[]>> table = new HashMap<String,Map.Entry<String, byte[]>>();

	/**
	 * A session-scoped CS transport
	 */
	@Transient@XmlTransient
	public Transport defaultTransport;

	@XmlElement(name = "developer")
	@Transient
	public Boolean isDeveloperSession; // Enables the features visible to developers only

	@OneToMany(mappedBy = "session", fetch = FetchType.LAZY)
	@XmlTransient
	public List<UserEvent> events = new ArrayList<UserEvent>();

	/**
	 * The flag to indicate unlimited quota for this session. Useful for physprop guest users.
	 */
	@Transient@XmlTransient
	public boolean disableQuota;

	public Session()
	{
		//		selectionList = new HashSet<Long>();
		//		selectionList.add(-1L);
		selectionMoleculeList = new HashSet<Long>();
		selectionMoleculeList.add(-1L);		
	}

	public String getLogin(){
		return user == null ? "guest" : user.login;
	}

	public static Session createNewSession(User user)
	{
		Session session = new Session();
		session.user = user;
		session.guid = UUID.randomUUID().toString();
		session.time = new Timestamp(Calendar.getInstance().getTimeInMillis());
		session.license = (user != null && user.licenseVersion >= Globals.CURRENT_LICENSE_VERSION);
		return session;
	}

	public void copySessionInfoFrom(Session session)
	{
		visitedModelsIds = session.visitedModelsIds;
		selectedAlerts = session.selectedAlerts;
		selectionMoleculeList = session.selectionMoleculeList;
		isModerator = session.isModerator;
		isDeveloperSession = session.isDeveloperSession;
		disableQuota = session.disableQuota;

		// Check if we have a specific mocked transport for this user
		if (user != null && user.isTestUser())
		{
			Transport sessionSpecificTransport = (Transport) UserAttachment.get(user, "cs-transport");
			if (sessionSpecificTransport != null)
				TransportFactory.setThreadTransport(sessionSpecificTransport);
		}
	}

	@XmlElement(name = "limit")
	public int getLimit()
	{
		if (user == null)
			return 2500;
		else
			return user.isOCHEMDeveloper() || user.isSuperUser() ? 
					QSPRConstants.SUPERUSER_RECORDS_LIMIT : user.isValidated() ?
							QSPRConstants.VALIDATEDUSER_RECORDS_LIMIT : QSPRConstants.USER_RECORDS_LIMIT;
	}

	public transient Boolean isModerator;

	@XmlAttribute(name = "moderator")
	private Boolean getModerator()
	{
		if (isModerator == null)
			isModerator =  AccessChecker.isModerator(user);

		return isModerator ? true : null;
	}

	@XmlAttribute(name = "max-priority")
	public int getMaxPriority()
	{
		return AccessChecker.getMaximumTaskPriority(user);
	}

	@XmlElement(name="activity-time")
	public String getActivityTime()
	{
		SimpleDateFormat sdf = new SimpleDateFormat("dd-MM-yyyy HH:mm:ss");
		return sdf.format(time);
	}

	@XmlElement
	public Boolean isMobile()
	{
		if (userAgent != null && userAgent.toLowerCase().contains("mobile"))
			return true;
		return null;
	}

	public static Session getById(long id)
	{
		return (Session) Globals.session().get(Session.class, id);
	}

	public static Session getByGUID(String guid)
	{
		return (Session) Globals.session().createCriteria(Session.class).add(Restrictions.eq("guid", guid)).uniqueResult();
	}

	public static Session getSuperuserSession()
	{
		return (Session) Globals.session().createCriteria(Session.class).createAlias("user", "u").add(Restrictions.eq("u.rank", 10)).setMaxResults(1).uniqueResult();
	}

	public static Session getFirstSession(String login)
	{
		return (Session) Globals.session().createCriteria(Session.class).createAlias("user", "u").add(Restrictions.eq("u.login", login)).setMaxResults(1).uniqueResult();
	}

	public static Session getLastSession(String login)
	{
		return (Session) Globals.session().createCriteria(Session.class).createAlias("user", "u").add(Restrictions.eq("u.login", login)).addOrder(Order.desc("id")).setMaxResults(1).uniqueResult();
	}

	public void setIpAddress(String ip) {
		if(ip != null && ip.length() > 20)
			ip = ip.substring(0,20);
		ipAddress = ip;
	}	
}
