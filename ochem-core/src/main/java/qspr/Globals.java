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

package qspr;

import java.io.BufferedWriter;
import java.io.File;
import java.io.StringWriter;
import java.sql.Timestamp;
import java.text.Format;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpSession;
import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;

import org.apache.commons.lang3.mutable.Mutable;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.flywaydb.core.Flyway;
import org.hibernate.Criteria;
import org.hibernate.Hibernate;
import org.hibernate.NonUniqueObjectException;
import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.hibernate.Transaction;
import org.hibernate.boot.registry.StandardServiceRegistryBuilder;
import org.hibernate.cfg.Configuration;
import org.hibernate.criterion.DetachedCriteria;
import org.hibernate.criterion.Disjunction;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;
import org.hibernate.criterion.Subqueries;
import org.hibernate.event.service.spi.EventListenerRegistry;
import org.hibernate.event.spi.EventType;
import org.hibernate.internal.SessionFactoryImpl;
import org.hibernate.proxy.HibernateProxy;
import org.hibernate.service.ServiceRegistry;
import org.hibernate.stat.Statistics;
import org.springframework.web.multipart.MultipartFile;
import org.springframework.web.multipart.MultipartHttpServletRequest;

import qspr.configuration.JAXBContextFactory;
import qspr.dao.ChemInfEngine;
import qspr.dao.Repository;
import qspr.dao.Various;
import qspr.entities.Basket;
import qspr.entities.Journal;
import qspr.entities.Mapping1;
import qspr.entities.Property;
import qspr.entities.Tag;
import qspr.entities.User;
import qspr.util.AnonymousUpdateListener;
import qspr.util.InsertListener;
import qspr.util.PreCollectionUpdateListener;
import qspr.util.UpdateListener;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.exceptions.CriticalException;
import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.OCHEMUtils;
import com.eadmet.utils.config.GlobalConfigurator;
import com.eadmet.utils.mailer.Mailer;

public class Globals
{
	static public double ERROR_RESULT = -999;
	static public boolean ALLOW_ESTIMATED_VALUES = false;
	static public boolean uniqueLoginRequired = false;

	// Primary lets pass clathpass to called jars like this
	// Maybe sometimes we can think of a better way
	// Since *-wildcards are available only in j2se6
	// and we are forced by some unknow reason to use j2se5

	static public String classPath = "";


	static public String defaultSID = null;
	//	static public boolean hidePrivateMolecules = true;
	//	static public boolean dataDownloadRestrictions = true;
	static public String commonUploadDirectory;
	static public String tempWorkingDirectory;
	static public String commonDownloadDirectory;
	static protected String executableDirectory;

	// rt now its a temporary solution for schedule task.
	public static boolean check = true;
	static public boolean debug = true;
	static public boolean debugJaxb = false;
	static public int CURRENT_LICENSE_VERSION = 2;

	public static ClassLoader classLoader;

	public static final SessionFactory sessionFactory;
	public static final SessionFactory alternateSessionFactory;

	public static final ServiceRegistry serviceRegistry;
	public static final ServiceRegistry alternateServiceRegistry;

	public static Statistics stat;

	public static Format formatter = new SimpleDateFormat("HH:mm:ss");

	public static int RIGHTS_NONE = 0;  // private
	//public static int RIGHTS_READ = 1;  // never used
	//public static int RIGHTS_WRITE = 2; // never used
	public static int RIGHTS_FREELY_AVAILABLE = 3; // freely available

	public static JAXBContext jaxbContext;
	public static Configuration ochemConf;
	public static Configuration fragmentsConf;

	public Map<Long, String> operationStatus = new HashMap<Long, String>();
	private static Logger logger = LogManager.getLogger(Globals.class);

	public static final String ochemCfgDir = "/etc/ochem/";

	static
	{
		try
		{
			System.setProperty("https.nonProxyHosts", "*");
			System.setProperty("http.nonProxyHosts", "*");
			System.setProperty("socks.nonProxyHosts", "*");
			System.setProperty("socksNonProxyHosts", "*");
			System.setProperty("socksProxyHost", "");
			System.setProperty("java.net.useSystemProxies", "false");

			loadConfiguration();

			ChemInfEngine engine = OCHEMConfiguration.getCheminfEngine();
			Various.molecule = Various.getCheminfImpl(engine);

			if(Various.molecule == null) throw new CriticalException("No Cheminf engine is available");


			if (!OCHEMConfiguration.disableDB)
			{
				// Configure DB
				ochemConf = new Configuration().configure();
				ochemConf.addAnnotatedClass(getCurrentUserClass());
				fragmentsConf = new Configuration().configure("alternate.cfg.xml");
				for (String key : OCHEMConfiguration.mainDbOverrides.keySet())
					ochemConf.setProperty(key, OCHEMConfiguration.mainDbOverrides.get(key));
				for (String key : OCHEMConfiguration.fragmentDbOverrides.keySet())
					fragmentsConf.setProperty(key, OCHEMConfiguration.fragmentDbOverrides.get(key));

				serviceRegistry = new StandardServiceRegistryBuilder().applySettings(ochemConf.getProperties()).build();
				alternateServiceRegistry = new StandardServiceRegistryBuilder().applySettings(fragmentsConf.getProperties()).build();

				sessionFactory = ochemConf.buildSessionFactory(serviceRegistry);
				alternateSessionFactory = fragmentsConf.buildSessionFactory(alternateServiceRegistry);

				EventListenerRegistry registry = ((SessionFactoryImpl) sessionFactory).getServiceRegistry().getService(EventListenerRegistry.class);
				registry.getEventListenerGroup(EventType.PRE_UPDATE).appendListener(new UpdateListener());
				registry.getEventListenerGroup(EventType.PRE_UPDATE).appendListener(new AnonymousUpdateListener());
				registry.getEventListenerGroup(EventType.PRE_COLLECTION_UPDATE).appendListener(new PreCollectionUpdateListener());
				registry.getEventListenerGroup(EventType.POST_INSERT).appendListener(new InsertListener());

				Mailer.serverSignature = "Host: " + OCHEMConfiguration.getRootHost() + "\nDatabase connection: "
						+ ochemConf.getProperty("hibernate.connection.url");

				stat = sessionFactory.getStatistics();
				stat.setStatisticsEnabled(true);
				migrateDatabase();
				checkDatabaseSetting();

			}
			else
			{
				sessionFactory = null;
				alternateSessionFactory = null;
				serviceRegistry = null;
				alternateServiceRegistry = null;
			}

			// Decreasing Logger information for mongodb
			Logger mongoLogger = LogManager.getLogger( "org.mongodb.driver");
			Configurator.setLevel(mongoLogger.getName(), Level.ERROR);

			//mongoLogger.setLevel(Level.ERROR); // e.g. or Log.WARNING, etc. /// N.B. Check whether it works

			logger.info("QSPR core loaded");

		} catch (Throwable ex)
		{
			ex.printStackTrace();
			System.out.println("\n\n"+ ex.getMessage() + "\nExiting.");
			if(ex instanceof CriticalException) System.exit(1);
			throw new ExceptionInInitializerError(ex.getMessage());		}
	}

	public static boolean isStandalone = false;

	static public String getExecutableDir(String s)
	{
		if (executableDirectory == null)
			executableDirectory = Globals.class.getProtectionDomain().getCodeSource().getLocation().getPath().replace("ochem-core/target/classes/", "") + "ochem-webapp/src/main/webapp/WEB-INF/executables";

		return executableDirectory+"/"+s;
	}

	/**
	 * Check that everything is OK with basic settings after dump
	 * 
	 * @throws CriticalException 
	 */
	private static void checkDatabaseSetting() {

		Globals.startAllTransactions();
		User user = Repository.user.getById(QSPRConstants.PUBLISHER_ID);

		String message = "The database is incorrectly configured: execute addExtended.sql to fix this issue";

		if(user == null || !user.login.equals(QSPRConstants.PUBLISHER)) {
			if(User.extendedExists()) throw new CriticalException(message);
			else
				throw new CriticalException("The database is incorrectly configured and does not have user PUBLISHER = " + QSPRConstants.PUBLISHER + " with id: " + QSPRConstants.PUBLISHER_ID);
		}

		if(Repository.molecule.getEmptyMolecule() == null)
			throw new CriticalException(message);

		if(Globals.session().get(Journal.class, QSPRConstants.UNPUBLISHED_JOURNAL) == null)
			throw new CriticalException(message);

		Globals.rollbackAllTransactions();
	}

	static public void setExecutableDir(String path)
	{
		executableDirectory = path;
	}


	/**
	 * Migrate the database schema to the current version using FlyWay
	 */
	private static void migrateDatabase()
	{
		Flyway flyway = Flyway.configure().dataSource(Globals.ochemConf.getProperty("hibernate.connection.url"), 
				Globals.ochemConf.getProperty("hibernate.connection.username"), Globals.ochemConf.getProperty("hibernate.connection.password")).load();
		flyway.baseline();
		flyway.migrate();
	}

	private static void loadConfiguration() throws Exception
	{
		GlobalConfigurator configurator = new GlobalConfigurator();
		configurator.addResource(ochemCfgDir + "ochem.cfg");
		configurator.addResource(Environment.getWebInfDir() + "/conf/ochem.cfg");
		configurator.addResource(Globals.class.getClassLoader().getResource("ochem.cfg"));
		configurator.configure();

		StringWriter sw = new StringWriter();
		GlobalConfigurator.export(new BufferedWriter(sw));
		logger.info("========= Effective configuration (including default values) ===========");
		logger.info(sw.toString());
		logger.info("========================================================================");
	}

	static private void ensureDirectoryExists(String path)
	{
		File f = new File(path);
		if (!f.exists())
			f.mkdirs();
	}

	static public void initialize() throws Exception
	{
		// Initialize constants, using given context path

		// Common directories
		executableDirectory = Environment.getWebInfDir() + "/executables";
		ensureDirectoryExists(executableDirectory);

		tempWorkingDirectory = Environment.getWebInfDir() + "/temp";
		ensureDirectoryExists(tempWorkingDirectory);

		commonUploadDirectory = Environment.getWebInfDir() + "/upload";
		ensureDirectoryExists(commonUploadDirectory);

		commonDownloadDirectory = Environment.getRootDir() + "/exports";
		ensureDirectoryExists(commonDownloadDirectory);
	}

	static public Session session()
	{
		if (OCHEMConfiguration.disableDB)
			return new StubSession();
		return sessionFactory.getCurrentSession();
	}

	static public Session alternateSession()
	{
		if (OCHEMConfiguration.disableDB)
			return new StubSession();
		return alternateSessionFactory.getCurrentSession();
	}

	static public Basket selectionBasket(boolean trash)
	{
		Basket basket = Basket.getBasket(trash ? "Records in trash" : "Selected records");
		basket.basketType = 1L;
		Globals.session().saveOrUpdate(basket);
		return basket;
	}

	static public boolean isGuestUser()
	{
		return Globals.userSession().user == null;
	}

	static public boolean isSuperUser()
	{
		return Globals.userSession().user != null && Globals.userSession().user.isSuperUser();
	}

	static public boolean isValidatedUser()
	{
		return Globals.userSession().user != null && Globals.userSession().user.isValidated();
	}

	static public User myself() {
		return Globals.userSession().user;
	}

	static public String getUsername()
	{
		if (Globals.userSession() != null && Globals.userSession().user != null)
			return Globals.userSession().user.login;
		else
			return null;
	}

	static public boolean isOCHEMDeveloper() {
		User user = getCurrentUser();
		if(user == null) return false;
		return user.isOCHEMDeveloper();
	}

	static public User getCurrentUser()
	{
		qspr.entities.Session session = Globals.userSession();
		return session != null ? session.user : null;
	}

	static public Class<?> getCurrentUserClass() {
		try {
			return OCHEMConfiguration.getUserClass();
		} catch (ClassNotFoundException e) {
			throw new UserFriendlyException(e);
		}
	}

	static public qspr.entities.Session userSession()
	{
		qspr.entities.Session currentThreadSession = ThreadScope.get().userSession;
		HttpServletRequest req = ThreadScope.get().localRequest;

		// Refetch the session if its not in the current transaction scope
		//		Transaction tx = ThreadScope.get().transaction;
		//		if (currentThreadSession != null)
		//			if (tx != null && tx.isActive() && !Globals.session().contains(currentThreadSession))
		//			{
		//				//qspr.entities.Session oldSession = currentThreadSession;
		//				//currentThreadSession = (qspr.entities.Session) Globals.session().get(qspr.entities.Session.class, currentThreadSession.id);
		//				//currentThreadSession.copySessionInfoFrom(oldSession);
		//				//ThreadScope.get().userSession = currentThreadSession;
		//			}

		if (currentThreadSession == null)
		{
			HttpSession httpSession = null;
			if (req != null)
				httpSession = req.getSession();
			if (httpSession == null)
				return null;
			currentThreadSession = (qspr.entities.Session) httpSession.getAttribute("user-session");
			if (currentThreadSession != null)
			{
				qspr.entities.Session oldSession = currentThreadSession;
				currentThreadSession = (qspr.entities.Session) Globals.session().get(qspr.entities.Session.class, currentThreadSession.id);

				if (currentThreadSession == null)
					httpSession.setAttribute("user-session", null);
				else
				{
					currentThreadSession.copySessionInfoFrom(oldSession);
					httpSession.setAttribute("user-session", currentThreadSession);
					ThreadScope.get().userSession = currentThreadSession;	
				}
			}
		}

		if (currentThreadSession != null && !Boolean.FALSE.equals(ThreadScope.get().updateSessionTime))
		{
			currentThreadSession.time = new Timestamp(Calendar.getInstance().getTimeInMillis());

			if (req != null)
				currentThreadSession.setIpAddress(ThreadScope.resolveRemoteAddr(req));
			try
			{
				if (Globals.session().getTransaction() != null
						&& Globals.session().getTransaction().isActive())
					Globals.session().saveOrUpdate(currentThreadSession);
			} catch (NonUniqueObjectException e)
			{
				// Its ok.
			}
		}

		return currentThreadSession;

	}

	@SuppressWarnings({ "unchecked", "rawtypes" })
	static public Set<Tag> getTaginationFilters(Class targetClass)
	{
		if (ThreadScope.get().localRequest == null)
			return null;

		Set<Tag> allTags = (Set<Tag>) ThreadScope.get().localRequest.getSession().getAttribute("global-filters");
		if (allTags == null)
		{
			allTags = new HashSet<Tag>();
			ThreadScope.get().localRequest.getSession().setAttribute("global-filters", allTags);
		}
		if (targetClass == null)
			return allTags;

		Set<Tag> result = new HashSet<Tag>();
		for (Tag tag : allTags)
		{
			if (targetClass == Property.class ^ "molecule".equals(tag.type))
				result.add(tag);
		}
		return result;
	}

	@SuppressWarnings("rawtypes")
	static public boolean applyTaginationFilters(Criteria criteria, Class targetClass, String targetAlias)
	{
		Set<Tag> globalFilters = Globals.getTaginationFilters(targetClass);
		if (globalFilters != null && globalFilters.size() > 0)
		{
			// The AND condition with tags is a problematic query
			// It is solved by replacing conjunction with disjunction
			// And adding additional restriction on count.
			// Done with a subquery, so performance is under question, to be benchmarked
			// Currently this is the only imaginable way to implement this
			// Midnighter
			if (targetClass == Mapping1.class)
			{
				// Detached query with Molecules works too long (> 80 sec for typical query)
				// So currenty only ONE tag for molecules is supported / Midnighter
				criteria.createAlias(targetAlias + ".tags", targetAlias + "_tg");
				criteria.add(Restrictions.eq(targetAlias + "_tg.id", globalFilters.iterator().next().id));
			}
			else
			{
				Disjunction disjunction = Restrictions.disjunction();
				for (Tag tag : globalFilters)
				{
					disjunction.add(Restrictions.eq("tg.id", tag.id));
				}
				DetachedCriteria tagCount = DetachedCriteria.forClass(targetClass, targetAlias + "_inner").add(disjunction)
						.add(Restrictions.eqProperty(targetAlias + "_inner.id", targetAlias + ".id")).createAlias(targetAlias + "_inner.tags", "tg")
						.setProjection(Projections.countDistinct("tg.id"));

				criteria.add(Subqueries.leAll(Long.valueOf(globalFilters.size()), tagCount));
			}
			return true;
		}
		return false;
	}

	static public void userSession(qspr.entities.Session sess)
	{
		ThreadScope.get().localRequest.getSession().setAttribute("user-session", sess);
	}

	private static boolean isTransactionRunning(Mutable<Transaction> transaction)
	{
		Transaction tx = transaction.getValue();
		return (tx != null && tx.isActive());
	}

	private static void startTransaction(Mutable<Transaction> transaction, SessionFactory sf)
	{
		if (OCHEMConfiguration.disableDB)
			return;
		if (isTransactionRunning(transaction))
		{
			logger.error("Trying to start transaction while another one is already active");
			logger.error(new Throwable());
			Exception e = new Exception();
			e.printStackTrace();
		}
		else
			transaction.setValue(sf.getCurrentSession().beginTransaction());
	}

	private static void rollbackTransaction(Mutable<Transaction> transaction, SessionFactory sf)
	{
		if (OCHEMConfiguration.disableDB)
			return;
		Transaction tx = transaction.getValue();
		if (tx == null)
			logger.error("Trying to rollback the null transaction");
		else if (!tx.isActive())
			logger.error("Trying to rollback the inactive transaction");
		else
		{
			transaction.setValue(null);
			tx.rollback();
		}
	}

	private static void commitTransaction(Mutable<Transaction> transaction, SessionFactory sf)
	{
		if (OCHEMConfiguration.disableDB)
			return;
		Transaction tx = transaction.getValue();
		if (tx == null)
			logger.error("Trying to commit the null transaction");
		else if (!tx.isActive())
			logger.error("Trying to commit the inactive transaction");
		else
		{
			transaction.setValue(null);
			try
			{
				sf.getCurrentSession().flush();
				tx.commit();
			} catch (Exception e)
			{
				tx.rollback();
				throw new RuntimeException(e);
			}
		}
	}

	private static void restartTransaction(Mutable<Transaction> transaction, SessionFactory sf, boolean commit)
	{
		if (commit)
			commitTransaction(transaction, sf);
		else
			rollbackTransaction(transaction, sf);
		startTransaction(transaction, sf);
	}

	static public boolean isMainTransactionRunning()
	{
		return isTransactionRunning(ThreadScope.get().transaction);
	}

	static public void startMainTransaction()
	{
		logger.debug("Starting main transaction");
		startTransaction(ThreadScope.get().transaction, Globals.sessionFactory);
	}

	static public void commitMainTransaction()
	{
		logger.debug("Committing main transaction");
		commitTransaction(ThreadScope.get().transaction, Globals.sessionFactory);
	}

	static public void rollbackMainTransaction()
	{
		logger.debug("Rolling back main transaction");
		rollbackTransaction(ThreadScope.get().transaction, Globals.sessionFactory);
	}

	static public void restartMainTransaction(boolean commit)
	{
		logger.debug("Restarting main transaction, commit="+commit);
		restartTransaction(ThreadScope.get().transaction, Globals.sessionFactory, commit);
	}

	static public boolean isAlternateTransactionRunning()
	{
		return isTransactionRunning(ThreadScope.get().alternateTransaction);
	}

	static public void startAlternateTransaction()
	{
		logger.debug("Starting alternate transaction");
		startTransaction(ThreadScope.get().alternateTransaction, Globals.alternateSessionFactory);
	}

	static public void commitAlternateTransaction()
	{
		logger.debug("Committing alternate transaction");
		commitTransaction(ThreadScope.get().alternateTransaction, Globals.alternateSessionFactory);
	}

	static public void rollbackAlternateTransaction()
	{
		logger.debug("Rolling back alternate transaction");
		rollbackTransaction(ThreadScope.get().alternateTransaction, Globals.alternateSessionFactory);
	}

	static public void restartAlternateTransaction(boolean commit)
	{
		logger.debug("Restarting alternate transaction, commit="+commit);
		restartTransaction(ThreadScope.get().alternateTransaction, Globals.alternateSessionFactory, commit);
	}

	static public boolean areTransactionsRunning()
	{
		return (isTransactionRunning(ThreadScope.get().transaction) && isTransactionRunning(ThreadScope.get().alternateTransaction));
	}

	static public void startAllTransactions()
	{
		logger.debug("Starting all transaction");
		startTransaction(ThreadScope.get().transaction, Globals.sessionFactory);
		startTransaction(ThreadScope.get().alternateTransaction, Globals.alternateSessionFactory);
	}

	static public void commitAllTransactions()
	{
		logger.debug("Committing all transactions");
		Exception exception = null;
		try
		{
			commitTransaction(ThreadScope.get().transaction, Globals.sessionFactory);
		} catch (Exception e)
		{
			logger.debug(e);
			exception = e;
		}

		try
		{
			commitTransaction(ThreadScope.get().alternateTransaction, Globals.alternateSessionFactory);
		} catch (Exception e)
		{
			logger.debug(e);
			exception = e;
		}

		if (exception != null)
			throw new RuntimeException(exception);
	}

	static public void rollbackAllTransactions()
	{
		logger.debug("Rolling back all transactions");
		Exception exception = null;
		try
		{
			rollbackTransaction(ThreadScope.get().transaction, Globals.sessionFactory);
		} catch (Exception e)
		{
			logger.debug(e);
			exception = e;
		}

		try
		{
			rollbackTransaction(ThreadScope.get().alternateTransaction, Globals.alternateSessionFactory);
		} catch (Exception e)
		{
			logger.debug(e);
			exception = e;
		}

		if (exception != null)
			throw new RuntimeException(exception);
	}

	static public void restartAllTransactions(boolean commit)
	{
		logger.debug("Restarting all transactions, commit="+commit);
		Exception exception = null;
		try
		{
			restartTransaction(ThreadScope.get().transaction, Globals.sessionFactory, commit);
		} catch (Exception e)
		{
			logger.debug(e);
			exception = e;
		}

		try
		{
			restartTransaction(ThreadScope.get().alternateTransaction, Globals.alternateSessionFactory, commit);
		} catch (Exception e)
		{
			logger.debug(e);
			exception = e;
		}

		if (exception != null)
			OCHEMUtils.rethrowSafely(exception);
	}


	static public Object getSessionAttribute(SessionVariable var)
	{
		if (ThreadScope.get().localRequest == null || ThreadScope.get().localRequest.getSession() == null)
			return null;
		return ThreadScope.get().localRequest.getSession().getAttribute(var.toString());
	}

	static public void setSessionAttribute(SessionVariable var, Object obj)
	{
		ThreadScope.get().localRequest.getSession().setAttribute(var.toString(), obj);
	}

	public static Statistics getHibernateStatistics()
	{
		return stat;
	}

	public static boolean considerPredicates()
	{
		return Boolean.TRUE.equals(ThreadScope.get().considerPredicatesInStatistics);
	}

	public static String htmlEntityDecode(String s)
	{
		return s == null ? null : s.trim();
	}

	public static boolean getMarshallingOption(MarshallingOption option)
	{
		List<MarshallingOption> options = ThreadScope.get().marshallingOptions;
		return options != null && options.contains(option);
	}

	public static void setMarshallingOption(MarshallingOption option)
	{
		List<MarshallingOption> options = ThreadScope.get().marshallingOptions;
		if (options == null)
			ThreadScope.get().marshallingOptions = options = new ArrayList<MarshallingOption>();
		options.add(option);
	}

	public static void removeMarshallingOption(MarshallingOption option)
	{
		List<MarshallingOption> options = ThreadScope.get().marshallingOptions;
		if (options == null)
			return;
		options.remove(option);
	}

	public static void clearMarshallingOptions()
	{
		ThreadScope.get().marshallingOptions = null;
	}

	public static String getClientSID()
	{
		if (defaultSID != null)
			return defaultSID;

		// Construct the SID on basis of IP address and user name
		String sid = "";
		if (ThreadScope.get().localRequest != null)
			sid += ThreadScope.resolveRemoteAddr(ThreadScope.get().localRequest);
		if (ThreadScope.get().context != null)
			sid += "/" + ThreadScope.get().context;
		if (Globals.userSession() != null && Globals.userSession().user != null)
			sid += "/" + Globals.userSession().user.login;

		if ("".equals(sid))
			sid = QSPRConstants.ANONYMOUS;
		return sid;
	}

	private static File getFileFromMultipart(MultipartFile mf) throws Exception
	{
		String subdir = QSPRConstants.ANONYMOUS;
		if (Globals.userSession().user != null)
			subdir = Globals.userSession().user.id.toString();

		File dir = new File(Globals.commonUploadDirectory +"/" + subdir);

		if (!dir.exists())
			dir.mkdir();

		File f = new File(dir.getAbsoluteFile() + "/" + mf.getOriginalFilename());

		if (f.exists())
			f.delete();

		if (mf.getSize() == 0)
			throw new Exception("The file is zero bytes (" + dir.getAbsoluteFile() + "/" + mf.getOriginalFilename() + ")");

		mf.transferTo(f);
		return f;
	}

	public static File getUploadedFile() throws Exception
	{
		MultipartHttpServletRequest mp = ThreadScope.get().localMpRequest;

		if (mp == null)
			throw new Exception("No files were uploaded");

		Iterator<String> it = mp.getFileNames();

		if (!it.hasNext())
			throw new Exception("No files were uploaded");

		MultipartFile mf = mp.getFile(it.next());

		return getFileFromMultipart(mf);
	}

	public static File getUploadedFile(String name) throws Exception
	{
		MultipartHttpServletRequest mp = ThreadScope.get().localMpRequest;

		if (mp == null)
			return null;

		MultipartFile mf = mp.getFile(name);

		if (mf == null || mf.isEmpty())
			return null;

		return getFileFromMultipart(mf);
	}

	/**
	 * A voodoo solution for unproxying javassist interfaces sometimes created by Hibernate where undesired
	 * @param entity
	 * @return
	 */
	@SuppressWarnings({ "unused", "unchecked" })
	public static <T> T initializeAndUnproxy(T entity) {
		if (entity == null) {
			throw new 
			NullPointerException("Entity passed for initialization is null");
		}

		T oldEntity = entity;

		Hibernate.initialize(entity);
		if (entity instanceof HibernateProxy) {
			entity = (T) ((HibernateProxy) entity).getHibernateLazyInitializer()
					.getImplementation();
		}
		return entity;
	}

	public static Marshaller createMarshaller(boolean formattedOutput) throws JAXBException
	{
		Marshaller m = Globals.jaxbContext.createMarshaller();
		if (formattedOutput)
			m.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, Boolean.TRUE);
		return m;
	}


	public static boolean isLocalhost()
	{
		return OCHEMConfiguration.getRootHost().startsWith("http://localhost:8080");
	}

	/**
	 * Is the request coming form a trusted IP address?
	 * (if so, some privileged operations might be allowed)
	 */
	public static boolean isTrustedAddress()
	{
		if (ThreadScope.get().localRequest == null)
			return false;
		String ipAddress = ThreadScope.resolveRemoteAddr(ThreadScope.get().localRequest);
		return OCHEMConfiguration.isTrustedAddress(ipAddress);
	}

	public static String now()
	{
		Calendar cal = Calendar.getInstance();
		SimpleDateFormat sdf = new SimpleDateFormat("dd-MM-yyyy HH:mm:ss");
		return sdf.format(cal.getTime());
	}

	public static JAXBContext createJAXBContext() throws JAXBException
	{
		return JAXBContextFactory
				.get("qspr.controllers:qspr.frontend:qspr.entities:qspr.fragmententities:qspr.workflow.datatypes:qspr.workflow.structure:qspr.modelling:qspr.modelling.configurations:qspr.metaserver.configurations:qspr.export:qspr.batchupload:qspr.util:qspr.business:com.eadmet:qspr.services");
	}

}
