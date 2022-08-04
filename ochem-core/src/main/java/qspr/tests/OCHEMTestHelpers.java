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

package qspr.tests;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.concurrent.TimeoutException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import com.eadmet.utils.MAILERConstants;

import qspr.Globals;
import qspr.ThreadScope;
import qspr.configuration.JAXBContextFactory;
import qspr.dao.Repository;
import qspr.entities.Article;
import qspr.entities.Author;
import qspr.entities.Basket;
import qspr.entities.BasketEntry;
import qspr.entities.BatchUpload;
import qspr.entities.BatchUploadRow;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Journal;
import qspr.entities.Mapping1;
import qspr.entities.Model;
import qspr.entities.Molecule;
import qspr.entities.Predicate;
import qspr.entities.Property;
import qspr.entities.PropertyOption;
import qspr.entities.Session;
import qspr.entities.Tag;
import qspr.entities.Unit;
import qspr.entities.UnitCategory;
import qspr.entities.User;
import qspr.export.ExportableModel;
import qspr.metaserver.CalculationClient;
import qspr.metaserver.protocol.Task.TaskPriority;
import qspr.metaserver.transport.MockTransport;
import qspr.metaserver.transport.TransportFactory;
import qspr.modelling.ModelFactory;
import qspr.modelling.ModelProcessor;
import qspr.modelling.ModelStatistics;
import qspr.util.DynaWrap;
import qspr.workflow.utils.QSPRConstants;

/**
 * Miscellaneouous utilities for testing:
 * - Random data generation (user, property, basket, etc)
 * - Cleaning up after the test run
 * 
 * @author rkoerner, midnighter
 *
 */
@SuppressWarnings("unchecked")
public class OCHEMTestHelpers
{

	private static transient Logger logger = LogManager.getLogger(OCHEMTestHelpers.class);

	protected static final String TEST_CONDITION = "TestCondition_";
	private static final String TEST_ARTICLE_VOL = "666";
	protected static final String TEST_ARTICLE_TITLE = "TestArticle";
	protected static final String TEST_AUTHOR_NAME = "TestAuthor";
	protected static final String TEST_PROPERTY = "TestProperty_";
	public static final String TEST_USER_PREFIX = QSPRConstants.TEST_USER_PREFIX;

	public enum testModel{ ann_bagging_oestate, ann_cv_class, ann_cv_oestate, ann_cv_alogps };

	public enum loginUserType {guest, user, admin}

	public static Session generateRandomTestUser(String prefix) throws Exception
	{
		return generateRandomTestUser(prefix, loginUserType.user);
	}

	public static Session generateRandomTestUser(String prefix, loginUserType userType) throws Exception
	{
		String randNum = "_" + getRand();

		User user = null;
		if (User.extendedExists() && User.getExtended() == User.getCurrentClass()) {
			user = (User) User.getExtended().newInstance();
			DynaWrap wrapped = user.dynaWrapped();
			wrapped.setField("title", "Mr.");
			wrapped.setField("firstName", "Test");
			wrapped.setField("lastName", "User");
			wrapped.setField("organisation", "Academic");
			wrapped.setField("email", randNum.replaceFirst("_", "") + MAILERConstants.EMAIL_OCHEM);
			User.getExtended().getMethod("setPassword", String.class).invoke(user, randNum);
		} else {
			user = User.getNewInstance();
		}
		user.login = TEST_USER_PREFIX + prefix + randNum;
		user.licenseVersion = Globals.CURRENT_LICENSE_VERSION;
		switch (userType)
		{
		case user:
			user.rank = User.NOTVALIDATED;
			break;
		case admin:
			user.rank = User.SUPER;
			break;
		default:
			throw new RuntimeException("unknow user type selected. Please choose 'user' or 'admin'");
		}

		Globals.session().save(user);

		ThreadScope.get().userSession = Session.createNewSession(user);
		Globals.session().save(ThreadScope.get().userSession);

		logger.info("User with login " + user.login + " created");

		return ThreadScope.get().userSession;
	}

	public static Basket generateBasketFromFile(String basketName, Property property, InputStream is) throws IOException, NumberFormatException, TimeoutException
	{
		MoleculesWithValues mols = MoleculesWithValues.loadFromFile(is);
		return createBasket(mols.molecules, mols.values, property, basketName);
	}

	public static Basket generateBasketFromFile(String basketName, Property property, String fileName) throws IOException, NumberFormatException, TimeoutException
	{
		InputStream is = new FileInputStream(new File(fileName));
		return generateBasketFromFile(basketName, property, is);
	}

	public static Basket generateRandomBasketNumeric(String basketName, int size)
	{
		return generateRandomBasketNumeric(basketName, size, null);
	}

	public static Basket generateRandomBasketNumeric(String basketName, int size, Property property)
	{
		//List<Integer> idss = Globals.session().createSQLQuery("select Molecule_id from Molecule limit 100000," + size).list();
		List<Integer> idss = Globals.session().createSQLQuery("select Molecule_id from Molecule order by RAND() limit " + size).list();
		List<Long> ids = new ArrayList<Long>();
		for (Integer integer : idss)
		{
			ids.add(new Long(integer.longValue()));
		}

		List<Molecule> mols = Globals.session().createCriteria(Molecule.class).add(Restrictions.in("id", ids)).list();
		List<Double> values = new ArrayList<Double>(mols.size());
		for (Molecule molecule : mols)
		{
			values.add(molecule.molWeight);
		}

		if (property == null)
			property = generateRandomPropertyCondition(Property.TYPE_NUMERIC, false);

		return createBasket(mols, values, property, basketName);
	}

	/**
	 * @param mols
	 */
	private static Basket createBasket(List<Molecule> mols, List<Double> values, Property property, String basketName)
	{
		logger.info("Generating basket...");
		String freeBasketName = Basket.getFreeBasketName(basketName);
		Basket randomBasket = Basket.getBasket(freeBasketName, ThreadScope.get().userSession);

		Article article = generateRandomArticle();

		if (property == null)
			property = generateRandomPropertyCondition(Property.TYPE_NUMERIC, false);

		Predicate predicate = Predicate.get("=");

		List<ExperimentalProperty> eps = new ArrayList<ExperimentalProperty>(mols.size());

		for (int i = 0; i < mols.size(); i++)
		{
			ExperimentalProperty ep = new ExperimentalProperty(mols.get(i));
			ep.article = article;
			ep.predicate = predicate;
			ep.property = property;
			ep.introducer = ep.owner = Globals.userSession().user;

			switch (property.type)
			{
			case Property.TYPE_NUMERIC:
				ep.value = values.get(i);
				ep.unit = property.defaultUnit;
				break;
			case Property.TYPE_QUALITATIVE:
				ep.option = ep.property.options.get(i % 2);
				break;
			case Property.TYPE_TEXTUAL:
				//					propCond.unitCategory = UnitCategory.getByName("Dimensionless");
				//					propCond.defaultUnit = propCond.unitCategory.getDefaultUnit();
				break;	
			}

			ep.rights = Globals.RIGHTS_NONE;

			ep.artLineNum = i;
			ep.artMolId = String.valueOf(i);

			ep.canonicalValue = ep.value;

			Globals.session().save(ep);
			eps.add(ep);
		}

		for (int i = 0; i < eps.size(); i++)
		{
			BasketEntry be = new BasketEntry(eps.get(i));
			be.basket = randomBasket;
			randomBasket.entries.add(be);
		}

		Globals.session().saveOrUpdate(randomBasket);

		logger.info("Basket generated: " + randomBasket.name);

		return randomBasket;
	}

	public static Basket generateRandomBasketQualitative(String basketName, int size, Property property)
	{
		//List<Integer> idss = Globals.session().createSQLQuery("select Molecule_id from Molecule limit 100000," + size).list();
		List<Integer> idss = Globals.session().createSQLQuery("select Molecule_id from Molecule order by RAND() limit " + size).list();
		List<Long> ids = new ArrayList<Long>();
		for (Integer integer : idss)
		{
			ids.add(new Long(integer.longValue()));
		}
		List<Molecule> mols = Globals.session().createCriteria(Molecule.class).add(Restrictions.in("id", ids)).list();
		List<Double> values = new ArrayList<Double>(mols.size());
		for (Molecule molecule : mols)
		{
			values.add(molecule.molWeight);
		}

		if (property == null)
			property = generateRandomPropertyCondition(Property.TYPE_QUALITATIVE, false);

		return createBasket(mols, values, property, basketName);
	}


	/**
	 * 
	 * @return A random article
	 */
	public static Article generateRandomArticle()
	{
		User user = ThreadScope.get().userSession.user;

		String randNum = "_" + getRand();

		Article article = new Article(TEST_ARTICLE_TITLE + randNum);
		article.introducer = user;
		article.owner = user;
		article.mediaType = "article";
		article.pmid = 90000000L;
		article.isbn = "978-0-00-000000-2";
		article.volume = TEST_ARTICLE_VOL;
		article.publicationDate = new Date();

		Author author = new Author();
		author.lastName = TEST_AUTHOR_NAME + randNum;
		author.firstName = TEST_AUTHOR_NAME + randNum;
		article.addAuthor(author);
		Globals.session().save(author);

		Journal journal = new Journal();
		journal.setTitle("TestJournal" + randNum);
		journal.introducer = user;
		journal.owner = user;
		article.journal = journal;

		Globals.session().saveOrUpdate(journal);
		Globals.session().saveOrUpdate(article);

		logger.info("Article: " + article.getTitle() + " with author " + article.authors.get(0) + " in journal " + article.journal.getTitle() + " created");

		return article;
	}

	/**
	 * 
	 * @return
	 */
	protected static Property generateRandomPropertyCondition(int type, boolean isCondition)
	{
		String randNum = "_" + getRand();
		Property propCond = new Property((isCondition ? TEST_CONDITION : TEST_PROPERTY) + "type" + type + randNum);

		User user = ThreadScope.get().userSession.user;
		propCond.owner = user;
		propCond.introducer = user;
		propCond.moderator = user;
		propCond.isCondition = isCondition;

		switch (type)
		{
		case Property.TYPE_NUMERIC:
			Unit unit = Unit.getById(1L);
			propCond.unitCategory = unit.category;
			propCond.defaultUnit = unit;
			break;
		case Property.TYPE_QUALITATIVE:
			PropertyOption option = Repository.option.getPropertyOptionByName("TestOption1", propCond.id, true, false);
			propCond.options.add(option);
			option = Repository.option.getPropertyOptionByName("TestOption2", propCond.id, true, false);
			propCond.options.add(option);
			propCond.unitCategory = UnitCategory.getByName("Dimensionless");
			propCond.defaultUnit = propCond.unitCategory.getDefaultUnit();
			break;
		case Property.TYPE_TEXTUAL:
			propCond.unitCategory = UnitCategory.getByName("Dimensionless");
			propCond.defaultUnit = propCond.unitCategory.getDefaultUnit();
			break;	
		}

		propCond.type = type;

		Globals.session().save(propCond);
		logger.info("Random property / condition generated: " + propCond.getName());
		return propCond;
	}

	/**
	 * @return 
	 * 
	 */
	protected static Tag generateRandomTag()
	{
		Tag tag = new Tag();
		tag.name = "TestTag_" + getRand();
		tag.type = "molecule";

		User user = ThreadScope.get().userSession.user;
		tag.owner = user;
		tag.introducer = user;

		Globals.session().save(tag);
		return tag;
	}

	public static Model trainAModel(final String PREFIX, final Basket trainingSet, final Basket validationSet, final String testModel, boolean useMock, boolean storePendingModel) throws Exception
	{
		if (Globals.jaxbContext == null)
			Globals.jaxbContext = JAXBContextFactory.get("qspr.frontend:qspr.entities:qspr.fragmententities:qspr.workflow.datatypes:qspr.workflow.structure:qspr.modelling:qspr.modelling.configurations:qspr.metaserver.configurations:qspr.export:qspr.util:qspr.business");

		ExportableModel eModel = null;
		if (testModel.startsWith("classpath:"))
			eModel = (ExportableModel) Globals.jaxbContext.createUnmarshaller().unmarshal(OCHEMTestHelpers.class.getClassLoader().getResourceAsStream(testModel.substring(10)));
		else
			eModel = (ExportableModel) Globals.jaxbContext.createUnmarshaller().unmarshal(new File(testModel));

		eModel.trainingSetId = trainingSet.id;


		eModel.validationSetId = new ArrayList<Long>();
		if (validationSet != null)
			eModel.validationSetId.add(validationSet.id);
		Model model = eModel.createModel();

		model.name = trainingSet.name + "_modelTest" + "_" + getRand();
		model.session = ThreadScope.get().userSession;

		model.defaultTaskPriority = TaskPriority.EXTRA_HIGH;

		if (useMock)
		{

			MockTransport mockTransport = new ModelCreationMockTransport();
			TransportFactory.setThreadTransport(mockTransport);
		}

		ModelProcessor processor = ModelFactory.getProcessor(model.template);
		processor.model = model;
		Globals.commitAllTransactions();
		processor.start();
		processor.join();

		if (processor.exception != null)
			throw processor.exception;

		if (useMock)
			TransportFactory.clearThreadTransport();

		Globals.startAllTransactions();

		processor.model = Repository.model.getById(processor.model.id);
		processor.saveModel();

		if (storePendingModel)
			processor.model.storePendingModel();

		ModelStatistics mStats = (ModelStatistics) processor.model.modelMappings.get(0).statisticsOriginal.getObject();
		mStats.recalculateStatistics(processor.model.modelMappings.get(0));

		return processor.model;		
	}

	/**
	 * Train a new model based on the configuration files from a specified folder
	 * 
	 * @param testModelFolder the folder with configuration files
	 * @return
	 * @throws Exception
	 */
	public static Model trainAModel(final String PREFIX, final String testBasketFolder, final String testModel, int propertyType, boolean useMock, boolean storePendingModel) throws Exception
	{
		Property property = generateRandomPropertyCondition(propertyType, false);
		Basket trainingSet = OCHEMTestHelpers.generateBasketFromFile(PREFIX, property, testBasketFolder + "/training-set.csv");
		Basket validationSet = OCHEMTestHelpers.generateBasketFromFile(PREFIX, trainingSet.getProperty().get(0), testBasketFolder + "/validation-set.csv");
		return trainAModel(PREFIX, trainingSet, validationSet, testModel, useMock, storePendingModel);
	}

	public static void cleanup(User deletedUser, boolean deleteUser)
	{
		logger.info("Cleaning up user " + deletedUser.login);
		assert deletedUser.login.startsWith(TEST_USER_PREFIX);

		List<Session> sessions = deletedUser.getAllSessions();
		for (int i = 0; i < sessions.size(); i++)
		{
			Session session = sessions.get(i);

			List<BatchUpload> bus = Globals.session().createQuery("from BatchUpload where session=:session").setParameter("session", session).list();
			for (BatchUpload bu : bus)
			{
				List<BatchUploadRow> burs = Globals.session().createQuery("from BatchUploadRow where batchupload_id=:bu_id").setParameter("bu_id", bu.id).list();
				for (BatchUploadRow batchUploadRow : burs)
					Globals.session().delete(batchUploadRow);
				logger.info("Deleting BatchUpload " + bu.id);
				Globals.session().delete(bu);
			}

			final List<Integer> taskIDs = Globals.session().createQuery("select taskId from PendingTask where session=:session and taskId is not null").setParameter("session", session).list();
			Globals.session().createQuery("delete from PendingTask where session=:session").setParameter("session", session).executeUpdate();
			if (!taskIDs.isEmpty())
			{
				new Thread()
				{
					@Override
					public void run()
					{
						CalculationClient client = new CalculationClient("Tests-Cleanup");
						client.setTolerateMetaserverDown();
						for (Integer taskID : taskIDs)
							if (taskID != 0)
								try
						{
									logger.info("Killing pending task with task ID " + taskID);
									client.killTask(taskID);
						} catch (Exception e)
						{
							logger.warn("Task ID " + taskID + " could not be killed");
						}
					}
				}.start();
			}

			// delete modelmapping, pendingTask and model
			Globals.session()
			.createSQLQuery("delete CachedPrediction from CachedPrediction, Model where Model.session_id=:sessionId and CachedPrediction.model_id=Model.model_id")
			.setParameter("sessionId", session.id).executeUpdate();
			Globals.session()
			.createSQLQuery("delete ModelMapping from ModelMapping, Model where Model.session_id=:sessionId and ModelMapping.model_id=Model.model_id")
			.setParameter("sessionId", session.id).executeUpdate();
			logger.info("Deleting model where session id was " + session.id);

			Globals.session().createSQLQuery("delete PendingTask from Model, PendingTask where Model.session_id=:sessionId and PendingTask.model_id=Model.model_id").setParameter("sessionId", session.id).executeUpdate();
			//Globals.session().createSQLQuery("delete from Model where Model.session_id=:sessionId").setParameter("sessionId", session.id).executeUpdate();
			List<Model> models = Globals.session().createCriteria(Model.class).add(Restrictions.eq("session", session)).list();
			for (Model model : models)
			{
				model.delete();
			}

			// delete actions 
			Globals.session().createSQLQuery("delete from ExportAction where ExportAction.session_id=:sessionId").setParameter("sessionId", session.id).executeUpdate();
			Globals.session().createSQLQuery("delete ua from UserAction ua join Action a using (action_id) where a.session_id=:sessionId").setParameter("sessionId", session.id).executeUpdate();
			Globals.session().createSQLQuery("delete from Action where session_id=:sessionId").setParameter("sessionId", session.id).executeUpdate();
		}

		// delete baskets and sessions after all models have been deleted
		for (int i = 0; i < sessions.size(); i++)
		{
			Session session = sessions.get(i);
			List<Basket> baskets = Globals.session().createCriteria(Basket.class).add(Restrictions.or(Restrictions.eq("user", deletedUser), Restrictions.eq("session", session))).list();
			for (Basket basket : baskets)
			{
				logger.info("Deleting basket " + basket.name + " (" + basket.id + ")");
				basket.delete();
			}
			logger.info("Deleting the session itself - "+session.id);
			Globals.session().delete(session);
		}

		// delete Properties and conditions
		List<Long> expPropertyIds = Globals.session().createCriteria(ExperimentalProperty.class).add(Restrictions.eq("introducer", deletedUser)).setProjection(Projections.id()).list();
		if (expPropertyIds.size() > 0)
		{
			Globals.session().createSQLQuery("update ExperimentalProperty set exp_property_id_conn=null where exp_property_id in :list").setParameterList("list", expPropertyIds).executeUpdate();
			Globals.session().createSQLQuery("delete from ExperimentalPropertyName where exp_property_id in :list").setParameterList("list", expPropertyIds).executeUpdate();
			Globals.session().createSQLQuery("delete from ExperimentalProperty where exp_property_id in :list").setParameterList("list", expPropertyIds).executeUpdate();
		}

		List<Property> props = Globals.session().createQuery("from Property where introducer=:user").setParameter("user", deletedUser).list();
		for (Property property : props)
		{
			Globals.session().createSQLQuery("delete from PropertyValue where property_id=:prop_id").setParameter("prop_id", property.id).executeUpdate();
			// This will also handle the options due to cascading
			Globals.session().delete(property);
			logger.info("Deleting property / condition " + property.getName() + " (" + property.id + ")");
		}

		Globals.session().createSQLQuery(
				"delete cs from ConditionSet cs " +
						"left join PropertyValue pv on cs.con_set_id = pv.con_set_id " +
						"left join ExperimentalProperty ep on cs.con_set_id = ep.con_set_id " + 
						"left join SubstructureAlert sa on cs.con_set_id = sa.con_set_id " + 
				"where pv.propertyvalue_id is null and ep.exp_property_id is null and sa.sa_id is null").executeUpdate();

		List<Article> articles = Globals.session().createQuery("from Article where introducer=:user").setParameter("user", deletedUser).list();
		for (Article article : articles)
		{
			for (Author author : article.authors)
			{
				logger.info("Deleting author " + author.lastName + " (" + author.id + ")");
				Globals.session().delete(author);
			}
			article.authors.clear();
			logger.info("Deleting article " + article.getTitle() + " (" + article.id + ")");
			Globals.session().delete(article);
		}

		List<Journal> journals = Globals.session().createQuery("from Journal where introducer=:user").setParameter("user", deletedUser).list();
		for (Journal journal : journals)
		{
			logger.info("Deleting journal " + journal.getTitle() + " (" + journal.id + ")");
			Globals.session().delete(journal);
		}

		List<Tag> tags = Globals.session().createQuery("from Tag where introducer=:user").setParameter("user", deletedUser).list();
		for (Tag tag : tags)
		{
			for (Article a : tag.articles)
				a.tags.remove(tag);

			for (Property p : tag.properties)
				p.tags.remove(tag);

			for (Object m : tag.mapping)
				((Mapping1)m).tags.remove(tag);

			Globals.session().delete(tag);
		}

		Globals.session().createQuery("delete from Unit where introducer=:user").setParameter("user", deletedUser).executeUpdate();
		Globals.session().createQuery("delete from SubstructureAlert where introducer=:user").setParameter("user", deletedUser).executeUpdate();
		Globals.session().createQuery("delete from UserAttachment where user=:user").setParameter("user", deletedUser).executeUpdate();

		if (deleteUser)
			Globals.session().delete(deletedUser);
		logger.info("Cleaning finished");
	}

	public static void cleanupByPrefix(String loginPrefix)
	{ 
		List<User> users = (List<User>) Globals.session().createCriteria(User.getCurrentClass()).add(Restrictions.like("login", TEST_USER_PREFIX + loginPrefix + "%"))
				.list();
		while (users.size() > 0)
		{
			cleanup(users.get(0), true);
			Globals.restartAllTransactions(true);
			users = (List<User>) Globals.session().createCriteria(User.class).add(Restrictions.like("login", TEST_USER_PREFIX + loginPrefix + "%")).list();
		}
	}

	public static String getRand()
	{
		return "" + Math.round(Math.random() * 10000);
	}

	public static void main(String[] args)
	{
		Globals.startAllTransactions();
		OCHEMTestHelpers.cleanupByPrefix("Selenium");
		Globals.commitAllTransactions();
	}
}
