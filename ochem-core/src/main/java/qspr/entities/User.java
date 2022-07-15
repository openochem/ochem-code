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

import java.util.ArrayList;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.FetchType;
import javax.persistence.GeneratedValue;
import javax.persistence.GenerationType;
import javax.persistence.Id;
import javax.persistence.Inheritance;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.OneToMany;
import javax.persistence.Table;
import javax.persistence.Transient;
import javax.persistence.InheritanceType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlElementWrapper;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;
import javax.xml.bind.annotation.adapters.XmlJavaTypeAdapter;

import org.hibernate.Criteria;
import org.hibernate.Hibernate;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.MarshallingOption;
import qspr.ThreadScope;
import qspr.annotations.Loggable;
import qspr.entities.Attachment.AttachmentType;
import qspr.util.DateAdapter;
import qspr.util.DynaWrap;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.messaging.DialogueService;

@Entity(name="User")
@Table(name="User")
@Inheritance(strategy = InheritanceType.JOINED)
//@MappedSuperclass
//@Inheritance(strategy = InheritanceType.TABLE_PER_CLASS)
@XmlRootElement(name = "user")
@Loggable
public class User {
	
	@Id
	@Column(name = "user_id")
//	@TableGenerator(
//		name = "table-generator", 
//	    table = "baseuser_id_generator", 
//	    pkColumnName = "seq_id", 
//	    valueColumnName = "seq_value",
//	    allocationSize = 1
//    )
//	@GeneratedValue(strategy = GenerationType.TABLE, generator = "table-generator")
	@GeneratedValue(strategy = GenerationType.IDENTITY)
	@XmlAttribute
	public Long id;

	@Column
	@XmlAttribute
	public String login;

	@ManyToOne
	@JoinColumn(name = "group_id")
	public Group group;

	@ManyToOne(fetch = FetchType.LAZY)
	@JoinColumn(name = "settings_id")
	@XmlTransient
	protected Attachment<UserSettings> userSettingsAttachment;

	@Column
	@XmlElement
	public Boolean suspended;

	@Column(name = "default_rights")
	public Integer defaultRights;

	@Column
	@XmlElement
	public Integer rank;

	@Transient
	@XmlElement
	public String message;

	@Column(name = "license_version")
	@XmlTransient
	public int licenseVersion;

	@Column(name = "activities_count_total")
	@XmlAttribute
	public Integer activitiesCountTotal;

	@Column(name = "activities_count_month")
	@XmlAttribute
	public Integer activitiesCountMonth;

	@Column(name = "latest_activity_time")
	@XmlJavaTypeAdapter(DateAdapter.class)
	public Date latestActivityTime;

	@Transient
	public String latestActivity; // XML only

	@OneToMany(mappedBy = "user", fetch = FetchType.LAZY)
	@XmlTransient
	public List<Session> sessions = new ArrayList<Session>();

	@Column(name = "facebook_id")
	public String facebookId;

	@Column(name = "referral_id")
	public Long referralId;

	final public static Integer SUPER = 10, VALIDATED = 1, NOTVALIDATED = 0, ANONYMOUS = -1;

	@XmlTransient
	public UserSettings getUserSettings()
	{
		if (userSettingsAttachment != null && userSettingsAttachment.getObject() != null)
			return userSettingsAttachment.getObject();
		else
			return new UserSettings();
	}

	public void setUserSettings(UserSettings settings)
	{
		if (userSettingsAttachment == null)
			userSettingsAttachment = new Attachment<UserSettings>(settings, AttachmentType.MARSHALABLE, AttachmentSource.UserSettings);
		else
			userSettingsAttachment.setObject(settings);
		userSettingsAttachment.updateObject();
		Globals.session().saveOrUpdate(userSettingsAttachment);
	}

	@XmlElement
	public User getReferral()
	{
		String uri = ThreadScope.get().requestURI;
		if (referralId != null && uri != null && (uri.contains("user/newuser") || uri.contains("user/show")))
			return User.getById(referralId);
		else
			return null;
	}


	@XmlAttribute(name="artcount")
	public Long getArticleCount()
	{
		if (Globals.getMarshallingOption(MarshallingOption.USER_RECORDS))
		{
			if((this.id !=null) && (this.id > 0)){
				Criteria criteria = Globals.session().createCriteria(Article.class)
						.add(Restrictions.or(Restrictions.eq("owner.id", this.id), Restrictions.eq("introducer.id", this.id)));
				criteria.setProjection(Projections.rowCount());
				return ((Long)criteria.list().get(0));
			}
		}
		return null;
	}

	/**
	 * Is this an automatically created user for the tests?
	 * @return
	 */
	@XmlTransient
	public boolean isTestUser()
	{
		return login.startsWith(QSPRConstants.TEST_USER_PREFIX);
	}

	@XmlAttribute(name="reccount")
	public Long getExpCount()
	{
		if (Globals.getMarshallingOption(MarshallingOption.USER_RECORDS))
		{
			if((this.id != null) && (this.id > 0)){
				Criteria criteria = Globals.session().createCriteria(ExperimentalProperty.class)
						.add(Restrictions.isNull("deleted"))
						.add(Restrictions.eq("introducer.id", this.id));

				ExperimentalProperty.addAccessRestrictions(criteria, Globals.RIGHTS_FREELY_AVAILABLE, this, false, true);

				criteria.setProjection(Projections.rowCount());
				return ((Long)criteria.list().get(0));
			}
		}
		return null;
	}

	@XmlElement(name="public-records-count")
	protected Long getPublicRecordsCount()
	{
		if (Globals.getMarshallingOption(MarshallingOption.USER_RECORDS))
		{
			if((this.id != null) && (this.id > 0)){
				Criteria criteria = Globals.session().createCriteria(ExperimentalProperty.class)
						.add(Restrictions.isNull("deleted"))
						.add(Restrictions.eq("rights", Globals.RIGHTS_FREELY_AVAILABLE))
						.add(Restrictions.eq("introducer.id", this.id));

				ExperimentalProperty.addAccessRestrictions(criteria, Globals.RIGHTS_FREELY_AVAILABLE, this, false, true);

				criteria.setProjection(Projections.rowCount());
				return ((Long)criteria.list().get(0));
			}
		}
		return null;
	}

	@XmlAttribute(name="basketcount")
	public Long getBasketCount()
	{
		if (Globals.getMarshallingOption(MarshallingOption.USER_RECORDS))
		{
			if((this.id !=null) && (this.id > 0)){
				Criteria criteria = Globals.session().createCriteria(Basket.class)
						.add(Restrictions.eq("user.id", this.id));
				criteria.setProjection(Projections.rowCount());
				return ((Long)criteria.list().get(0));
			}
		}
		return null;
	}


	private transient Long unreadMessages;

	@XmlAttribute(name="message")
	public Long getUnreadMessage()
	{
		if (Globals.userSession() == null || Globals.userSession().user != this)
			return null;
		if (unreadMessages != null)
			return unreadMessages;
		if((this.id !=null) && (this.id > 0)){
			return unreadMessages = new DialogueService().getUnreadMessagesCount(this);
		}
		return null;
	}

	public boolean isPublisher()
	{
		return this.id == QSPRConstants.PUBLISHER_ID;
	}

	@XmlElement(name = "superuser")
	public boolean isSuperUser()
	{
		return (this.rank != null && this.rank >= SUPER);
	}

	public boolean isValidated()
	{
		return (this.rank != null && this.rank >= VALIDATED);
	}

	public static User getByLogin(String login)
	{
		Criteria c = Globals.session().createCriteria(User.getCurrentClass());
		c.add(Restrictions.eq("login", login));
		return (User) c.uniqueResult();
	}

	public static User getById(long id)
	{
		return (User) Globals.session().get(User.getCurrentClass(), id);
	}

	public static User getByString(String id)
	{
		if (id.matches("[0-9]+"))
			return getById(Long.valueOf(id));
		else
			return User.getByLogin(id);
	}

	@SuppressWarnings("unchecked")
	public static User getByFacebookId(String facebookId)
	{
		Criteria c = Globals.session().createCriteria(User.getCurrentClass());
		c.add(Restrictions.eq("facebookId", facebookId));
		List<User> users = c.list();
		if (users.size() > 0)
			return users.get(0);
		else
			return null;

	}

	@SuppressWarnings("unchecked")
	public static User getByEmail(String email)
	{
		Criteria c = Globals.session().createCriteria(User.getCurrentClass());
		c.add(Restrictions.eq("email", email));
		List<User> users = c.list();
		if (users.size() > 0)
			return users.get(0);
		else
			return null;
	}

	@XmlElement(name="public-models-count")
	protected Long getPublicModelsCount()
	{
		if (!Globals.getMarshallingOption(MarshallingOption.USER_PUBLIC_MODELS))
			return null;

		if((this.id != null) && (this.id > 0))
		{
			Criteria criteria = Globals.session().createCriteria(Model.class)
					.add(Restrictions.eq("published", Boolean.TRUE))
					.createAlias("session", "s")
					.add(Restrictions.eq("s.user", this));

			criteria.setProjection(Projections.rowCount());
			return (Long)criteria.uniqueResult();
		}
		return null;
	}

	@SuppressWarnings("unchecked")
	@XmlTransient
	public List<Session> getAllSessions()
	{
		return Globals.session().createCriteria(Session.class).add(Restrictions.eq("user", this)).list();
	}

	@SuppressWarnings("unchecked")
	@XmlElementWrapper(name="public-models")
	@XmlElement(name = "model")
	protected List<Model> getPublicModels()
	{
		if (!Globals.getMarshallingOption(MarshallingOption.USER_PUBLIC_MODELS))
			return null;

		if((this.id != null) && (this.id > 0))
		{
			Criteria criteria = Globals.session().createCriteria(Model.class)
					.add(Restrictions.eq("published", Boolean.TRUE))
					.createAlias("session", "s")
					.add(Restrictions.eq("s.user", this));

			List<Model> models = criteria.list();
			for (Model model : models) {
				Hibernate.initialize(model.modelMappings);
				Globals.session().evict(model);
				model.session = null;
			}
			return models;
		}
		return null;
	}

	public static Long getActiveUsersCount()
	{
		return (Long) Globals.session().createCriteria(User.getCurrentClass())
				.add(Restrictions.or(Restrictions.eq("suspended", false), Restrictions.isNull("suspended")))
				.add(Restrictions.not(Restrictions.like("login", "test%")))
				.add(Restrictions.not(Restrictions.eq("id", 1l)))
				.setProjection(Projections.count("id"))
				.uniqueResult();
	}

	public String toString()
	{
		return this.login;
	}

	public boolean isSuspended()
	{
		return suspended != null && suspended;
	}

	public boolean equals(Object obj)
	{
		User user = (User) obj;
		return (user == null) ? false : ((id == null) ? (user.id == null) : id.equals(user.id));
	}

	public static String generatePassword(int length) {
		StringBuilder password = new StringBuilder();
		for (int i = 0; i < length; i++)
			password.append(new Character((char)('A' + Math.floor(Math.random() * ('z'- 'A')))));

		return password.toString();
	}

	public static void delete(User user) {
		Globals.session().delete(user);
	}
	
	public static Set<String> getDevelopers() {
		return new HashSet<String>();
	}
	
	private static final Class<?> currentClass = Globals.getCurrentUserClass();
	public static final String defaultExtendedClass = QSPRConstants.EXTENDED_USER;
	
	public static boolean extendedExists() {
		try {
			Class.forName(defaultExtendedClass);
			return true;
		} catch (Exception exp) {
			return false;
		}
	}
	
	public static Class<?> getExtended() {
		if (extendedExists()) {
			try {
				return Class.forName(defaultExtendedClass);
			} catch (ClassNotFoundException e) {
				throw new UserFriendlyException("Could not find extended User entity: " + defaultExtendedClass);
			}
		} else {
			throw new UserFriendlyException("Could not find extended User entity: " + defaultExtendedClass);
		}
	}
	
	public static Class<?> getCurrentClass() {
		return currentClass;
	}
	
	public static User getNewInstance() {
		try {
			return (User) currentClass.newInstance();
		} catch (InstantiationException | IllegalAccessException e) {
			throw new UserFriendlyException(e);
		}
	}
	
	public static User getNewInstance(String classPath) {
		try {
			return (User) (Class.forName(classPath).newInstance());
		} catch (InstantiationException | IllegalAccessException | ClassNotFoundException e) {
			throw new UserFriendlyException(e);
		}
	}
	
	public boolean isOCHEMDeveloper() {
		return false;
	}

	public String getFullName() {
		return login;
	}

	public boolean authorizeWithSecret(String secret) throws Exception {
		throw new UserFriendlyException("This user does not support login with a secret.");
	}

	public String getSecret() {
		// TODO Auto-generated method stub
		return null;
	}
	
	public DynaWrap dynaWrapped() {
		return new DynaWrap(this);
	}
	
	public boolean isExtended() {
		if (extendedExists()) {
			return getExtended().isInstance(this);
		} else {
			return false;
		}
	}

//	public abstract boolean isOCHEMDeveloper();
	
//	public abstract String getFullName();
	
//	public abstract String getSecret();
	
//	public abstract boolean authorizeWithSecret(String secret) throws Exception;
	
//	public abstract void setPassword(String passwd) throws Exception;

}
