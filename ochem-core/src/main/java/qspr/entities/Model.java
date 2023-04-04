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
import java.io.IOException;
import java.io.StringWriter;
import java.security.NoSuchAlgorithmException;
import java.sql.Timestamp;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import javax.persistence.CascadeType;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.FetchType;
import javax.persistence.GeneratedValue;
import javax.persistence.GenerationType;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.JoinTable;
import javax.persistence.ManyToMany;
import javax.persistence.ManyToOne;
import javax.persistence.OneToMany;
import javax.persistence.Transient;
import javax.servlet.http.HttpServletRequest;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlElementWrapper;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.Hibernate;
import org.hibernate.criterion.Conjunction;
import org.hibernate.criterion.Disjunction;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.SessionVariable;
import qspr.ThreadScope;
import qspr.annotations.Loggable;
import qspr.business.PendingTaskPeer;
import qspr.business.Privileges;
import qspr.dao.Repository;
import qspr.entities.Attachment.AttachmentType;
import qspr.interfaces.ProvidedConditions;
import qspr.metaserver.CalculationClient;
import qspr.metaserver.configurations.BaggingConfiguration;
import qspr.metaserver.configurations.ConsensusModelConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.configurations.NoDescriptors;
import qspr.metaserver.configurations.SupportsInversions;
import qspr.metaserver.configurations.SupportsOneOutputOnly;
import qspr.metaserver.configurations.ConsensusModelConfiguration.IndividualModel;
import qspr.metaserver.cs.WorkflowNodeServer;
import qspr.metaserver.configurations.CrossValidationConfiguration;
import qspr.metaserver.configurations.ValidationConfiguration;
import qspr.metaserver.configurations.ValidationConfiguration.MixtureValidation;
import qspr.metaserver.protocol.Task.TaskPriority;
import qspr.metaserver.util.ShortCondition;
import qspr.modelling.ModelFactory;
import qspr.modelling.ModelProcessor;
import qspr.modelling.ModelStatistics;
import qspr.modelling.PointStatistics;
import qspr.modelling.SetStatistics;
import qspr.modelling.applier.ModelApplier;
import qspr.modelling.configurations.CDSConfiguration;
import qspr.util.AccessChecker;
import qspr.util.BasicRecordMapper;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.OCHEMUtils;
import com.eadmet.utils.TimeUtils;
import com.eadmet.utils.mailer.Mailer;

@Entity
@XmlRootElement(name = "model")
public class Model 
{
	@Id
	@Column(name = "model_id")
	@GeneratedValue(strategy = GenerationType.AUTO)
	@XmlAttribute
	public Long id;

	@Column(name = "published_id")
	@XmlElement
	public Long publicId;

	@Column
	@XmlAttribute
	@Loggable
	public String name;

	/*
	 * The name shown to the user in the OCHEM predictor
	 */
	@Column(name = "featured_name")
	@XmlAttribute
	@Loggable
	public String featuredName;

	@XmlAttribute
	@Transient
	public Boolean selected; 

	@ManyToOne
	@JoinColumn(name = "model_template_id")
	public ModelTemplate template;

	@ManyToOne
	@JoinColumn(name = "training_set_id")
	@XmlElement(name = "training-set")
	public Basket trainingSet;

	@ManyToOne
	@JoinColumn(name = "validation_set_id")
	@XmlElement(name = "validation-set")
	public Basket validationSet;

	@XmlElement
	@Transient
	public Basket selectedValidationSet;

	@Transient
	@XmlTransient
	private Map<String, Basket> filteredSets = new TreeMap<String, Basket>();

	/**
	 * Contains configuration of the model
	 */

	@XmlTransient
	@ManyToOne(cascade={CascadeType.ALL}, fetch = FetchType.LAZY)
	@JoinColumn(name = "attachment")
	public Attachment<ModelAttachment> attachment;

	/**
	 * Contains calculated and recalculated model data
	 */

	@XmlTransient
	@ManyToOne(cascade={CascadeType.ALL}, fetch = FetchType.LAZY)
	@JoinColumn(name = "ready_model_attachment") 
	public Attachment<ReadyModelAttachment> readyModelAttachment;

	@XmlTransient
	@ManyToOne(cascade={CascadeType.ALL}, fetch = FetchType.LAZY)
	@JoinColumn(name = "calculated_descriptors")
	public Attachment<DataTable> calcDescriptors; 

	public DataTable getCalculatedDescriptors(){
		return calcDescriptors.getObject().deCompress();
	}

	/**
	 * Contains excluded basket entries
	 */
	@XmlTransient
	@ManyToOne(cascade={CascadeType.ALL}, fetch = FetchType.EAGER)
	@JoinColumn(name = "microattachment") 
	public Attachment<ModelMicroAttachment> microattachment;

	public Long size;

	@Column(name = "configuration_hash")
	public String configurationHash;

	@Column
	public boolean approved = true;

	@Column(name = "date_created")
	@XmlTransient
	public Timestamp dateCreated;

	/**
	 * An automatically generated description based on the model configuration
	 */
	@Column
	public String description; // A column used for display purposes only 

	/**
	 * A textual description manually entered by the user
	 */
	@Column(name = "user_description")
	public String userDescription;  

	@Column(name = "configuration_xml")
	public String configurationXml; // A column used for display purposes and runALL tasks test (to select specific models)

	@ManyToOne
	@JoinColumn(name = "session_id")
	public Session session;

	@OneToMany (cascade={CascadeType.REMOVE}, fetch = FetchType.LAZY, mappedBy = "model")
	@XmlElement
	public List<ModelMapping> modelMappings = new ArrayList<ModelMapping>();

	@Column
	public boolean published = false;

	@Column(name = "task_id")
	public Integer taskId;

	@Column(name = "ad_task_id")
	public Integer adTaskId;

	@Column(name = "prediction_scenario")
	@XmlElement
	public String predictionScenario = "0";	// see class ModelAbstractConfiguration.PredictionScenario {0 := PREDICTION_ONLY, 1 := TABLE_ONLY, 2 := PREDICTION_AND_TABLE}

	@Column
	public String status;

	@Column(name = "detailed_status")
	public String detailedStatus;

	@Column(name = "time_to_complete")
	public Long timeToComplete;

	@XmlTransient
	@Column(name = "last_modification")
	public Timestamp lastModification;

	@XmlTransient
	@Column(name = "last_access")
	public Timestamp lastAccess;

	@XmlTransient
	@Column(name = "deletion_warning_sent")
	public Timestamp deletionWarningSent; // the model is gonna be deleted. When did we notify a user?

	@XmlTransient
	@Column(name = "delete_after_days")
	public Integer deleteAfterDays; // how many days after a notification a model is kept? (first 14, in one week changed to 7)

	@XmlTransient
	@Column(name = "statistics_calculated")
	public boolean isStatisticsCalculated; // whether the statistics of a currently pending task has been already calculated

	@Column(name = "recalculation")
	public boolean recalculation; // whether there model is being recalculated

	//	@Transient
	//	@XmlTransient
	//	public ModelAbstractConfiguration uploadedModel; //This will store MethodSpecificData whatsoever.... Currently - only MLRAData... NoS 19.10.09
	//	
	@Transient
	@XmlTransient
	public Object preview; //This will store the Excel Sheet preview in a model-specific manner... Currently - only MLRA

	@Transient
	@XmlElement
	public int predictionErrors = 0; //Used in ModelApplierController to keep the per-model error count

	@Column(name = "task_priority")
	@XmlTransient
	public int defaultTaskPriority = TaskPriority.LOW;

	@Column(name = "full_equation")
	@XmlElement
	public String fullEquation;

	@Column(name = "column_list")
	public String columnList; // A helper field to store comma-separated columns names returned by model applier. Required to reconstruct a response from cache.

	@ManyToMany
	(
			targetEntity = Basket.class,
			cascade={CascadeType.PERSIST, CascadeType.MERGE},
			fetch = FetchType.EAGER			
			)
	@JoinTable
	(
			name="ValidationSet",
			joinColumns={@JoinColumn(name="model_id")},
			inverseJoinColumns={@JoinColumn(name="basket_id")}
			)
	@org.hibernate.annotations.IndexColumn(name="number")
	//@XmlElementWrapper(name = "val-sets")
	//@XmlElement(name = "validationSet")
	@XmlTransient
	public List<Basket> validationSets = new ArrayList<Basket>();

	@ManyToMany
	(
			targetEntity = Tag.class,
			cascade={CascadeType.PERSIST, CascadeType.MERGE, CascadeType.REMOVE},
			fetch = FetchType.LAZY
			)
	@JoinTable
	(
			name="ModelTag",
			joinColumns={@JoinColumn(name="model_id")},
			inverseJoinColumns={@JoinColumn(name="tag_id")}
			)
	@XmlTransient
	public Set<Tag> tags = new HashSet<Tag>();

	private static Logger logger = LogManager.getLogger(Model.class);

	private static final String errorModelChanged = "Error is stored predictions; some records were changed or model is very old and should be recalculated before it can be published.";

	private static final Long MAX_PUBLIC_ID = 1000000L;

	public Model()
	{
		assignRandomId();
		if (lastAccess == null)
			lastAccess = new Timestamp(Calendar.getInstance().getTimeInMillis());
		if (lastModification == null)
			lastModification = new Timestamp(Calendar.getInstance().getTimeInMillis());
		if (dateCreated == null)
			dateCreated = new Timestamp(Calendar.getInstance().getTimeInMillis());
	}


	private void checkIfConsensusSubModelsExistAndArePublished(Long articleId) throws Exception{
		ModelAttachment attach = attachment.getObject();
		if(attach.configuration instanceof ConsensusModelConfiguration){
			ConsensusModelConfiguration conf = (ConsensusModelConfiguration) attach.configuration;

			for( IndividualModel model:conf.individualModels) 
				if(Repository.model.getByPublicId(model.id) == null)
					throw new UserFriendlyException("Consensus model can't be used: one of its submodels with Public ID: "+ model.id +" has been deleted.");

			if(articleId == null) return;

			List<HashSet<Long>> perProperty = conf.getModelsForProperties();

			boolean save = false;
			for( IndividualModel model:conf.individualModels) {
				Model mod = Repository.model.getByPublicId(model.id);
				if(!mod.published) {
					mod.publish(articleId);
					if(perProperty != null && model.id != mod.publicId )for(HashSet<Long> s:perProperty) {
						if(s.contains(model.id)) {
							s.remove(model.id);
							s.add(mod.publicId);
						}
					}

					model.id = mod.publicId;
					save = true;

					if(perProperty != null) {
						conf.models = null;
						for(HashSet<Long> set:perProperty)
							conf.addModelsForProperty(new ArrayList<Long>(set));
					}

				}
			}
			if(save) {
				attachment.setObject(attach,AttachmentType.MARSHALABLE);
				attachment.updateObject();
				updateDescription();
				Globals.session().saveOrUpdate(attachment);
			}
		}
	}

	/**
	 * If model can't be applied, check will throw  Exception
	 */

	public void checkIfCanBeApplied() throws Exception{
		switch(template.name){
		case QSPRConstants.UPLOADED_MODEL: 	throw new Exception("Can't apply the uploaded model " + name);
		case QSPRConstants.CONSENSUS: checkIfConsensusSubModelsExistAndArePublished(null);
		}
	}

	/**
	 * If model can't be applied, check will throw  Exception
	 */

	public boolean isStandaloneExportable(){
		if(attachment.getObject().configuration instanceof ConsensusModelConfiguration)return false;
		return true;
	}

	public void publish(Long articleId)
	{
		if(articleId == null)
			throw new UserFriendlyException("You must specify an article to publish a model");

		article = (Article) Globals.session().get(Article.class,articleId);

		List<Property> properties = getNotPublishedProperties();
		if(properties.size()>0){
			String missed = "Your model uses non public or non approved by moderator properties: ";
			for(Property p: properties) {
				missed += "\""+ p.getName() + "\" ; ";
			}
			throw new UserFriendlyException(missed + " Publish these properties or request moderator " + QSPRConstants.INFOEMAIL + " to approve them before publishing the model.");
		}

		// Assign a new unique public ID (incremental for all the published models)
		if (article == null)
			throw new UserFriendlyException("You must specify an article to publish a model");

		if( (columnList == null || columnList.length() == 0) && !template.name.equals(QSPRConstants.UPLOADED_MODEL))
			throw new UserFriendlyException(errorModelChanged);

		try{
			checkIfConsensusSubModelsExistAndArePublished(articleId);
		}catch(Exception e){
			throw new UserFriendlyException(e.getMessage());
		}

		clearCache(); // delete cache each time model is published

		if (publicId == null || publicId > MAX_PUBLIC_ID)
		{
			Long maxPublicId = (Long) Globals.session().createCriteria(Model.class)
					.add(Restrictions.le("publicId", MAX_PUBLIC_ID))
					.setProjection(Projections.max("publicId")).uniqueResult();
			if (maxPublicId == null)
				maxPublicId = 0L;
			publicId = 1L + maxPublicId;
		}
		approved = false; // By default, published models are not approved
		published = true;
		isStatisticsCalculated = true;
		status = null;
		detailedStatus = "Finished";
		taskId = null;

		if(attachment != null)attachment.publish();
		if(microattachment != null)microattachment.publish();
		if(readyModelAttachment != null)readyModelAttachment.publish();
		if(calcDescriptors != null)calcDescriptors.publish();

		for(ModelMapping m: modelMappings){
			if(m.statisticsOriginal != null)m.statisticsOriginal.publish();
			if(m.statisticsRecalculated != null)m.statisticsRecalculated.publish();
		}
	}

	public void markAccessed()
	{
		lastAccess = new Timestamp(Calendar.getInstance().getTimeInMillis());
		deletionWarningSent = null;
		deleteAfterDays = null;
	}

	@XmlElementWrapper(name = "validation-sets")
	@XmlElement(name = "validation-set")
	public List<Basket> getValidationSets()
	{
		List<Basket> sets = new ArrayList<Basket>();
		if (validationSet != null)
			sets.add(validationSet);
		if (!validationSets.isEmpty())
			sets.addAll(validationSets);

		return sets;
	}

	/**
	 * Empirically determine class label by its name. 
	 * For example, "Yes" and "active" are "positive", while "insoluble" is "negative"
	 */
	public static Map<Long, Long> getClassificationRemapping(Property p) {
		Map<Long, Long> mapping = new HashMap<Long, Long>();
		List<String> positives = Arrays.asList(new String[]{"active", "high", "soluble", "inhibitor", "yes"});
		List<String> negatives = Arrays.asList(new String[]{"inactive", "low", "insoluble", "non", "no"});
		p = Property.getById(p.id);
		for (PropertyOption option : p.options)
		{
			for (String name : positives)
				if (option.name.toLowerCase().startsWith(name))
					mapping.put(option.id, 1L);
			for (String name : negatives)
				if (option.name.toLowerCase().startsWith(name))
					mapping.put(option.id, -1L);
			//if (!mapping.containsKey(option.id))
			//	throw new UserFriendlyException("Cannot identify option <" + option.name + "> as positive or negative");
		}

		return mapping;
	}

	public void addValidationSet(Basket basket)
	{
		if (!getValidationSets().contains(basket))
			if (validationSet == null)
				validationSet = basket;
			else
				validationSets.add(basket);
		logger.info("Validation set was added");
	}

	public void clearCache()
	{
		logger.info("Clearing cached predictions for model " + name);
		ModelApplier.clearCachedPredictions(publicId);
	}

	static public void deleteModelsByIds(List <Integer>ids){
		if(ids.size()>0){
			logger.info("Deleting "+ids.size() + " models");

			for(Integer m: ids){
				Model model = Repository.model.getById(m);
				logger.info("deleting anonymous model " + model.name + " id= " + model.id + " published_id=" + model.publicId);
				model.delete();
			}
		}
	}

	public void delete()
	{
		if (this.id != null)
		{

			logger.info("Deleting model "+id);

			if (published)
				throw new UserFriendlyException(
						"You cannot delete a published model! If you want to delete it, please contact the system administrator.");

			clearCache();

			@SuppressWarnings("unchecked")
			final List<Integer> taskIds = Globals.session().createQuery("select taskId from PendingTask where model=:model and taskId > 0 and status in ('assigned', 'init')").setParameter("model", this).list();
			PendingTaskPeer.terminateTaskAsync(taskIds);

			ModelApplier applier = (ModelApplier) Globals.getSessionAttribute(SessionVariable.MODEL_APPLIER);
			if (applier != null)
				applier.removeModel(this);
			//			Globals.session().createQuery("delete from ModelMapping where model=:model").setParameter("model", this).executeUpdate();
			Globals.session().createQuery("delete from PendingTask where model=:model").setParameter("model", this).executeUpdate();
			//			Globals.session().createSQLQuery("delete from ValidationSet where model_id=:id").setParameter("id", id).executeUpdate();
			//			Globals.session().createQuery("delete from Model where id=:id").setParameter("id", id).executeUpdate();
			try{
				Globals.session().delete(this);
			}catch(RuntimeException e){
				Mailer.notifyAdmins("problems to delete model " + publicId + " " +name, "The exception was " + e);
				attachment = null;
				readyModelAttachment = null;
				microattachment = null;
				Globals.session().delete(this);
			}
		}
	}

	public boolean hasValidationSets()
	{
		return validationSet != null || !validationSets.isEmpty();
	}

	public void cleanValidationSets()
	{
		validationSets.clear();
		validationSet = null;
	}

	@XmlAttribute
	public Boolean isToBeDeleted()
	{
		if (deleteAfterDays != null && !published)
			return true;
		return null;
	}

	public void assignRandomId()
	{
		// For every model generate a random publicId
		// everyone who knows it has the access to the model

		if (publicId == null)
			publicId = Math.round(Math.random() * 50000000d + QSPRConstants.MAXIMAL_PUBLIC_ID);
	}

	@Transient
	@XmlAttribute(name="url")
	protected String getURL()
	{
		return "model/action.do?action=discuss&id="+this.id;
	}

	@XmlTransient
	public Long getId()
	{
		return id;
	}

	@Transient
	@XmlAttribute(name="date")
	protected String getDate()
	{
		DateFormat df = new SimpleDateFormat ("yyyy-MM-dd");
		if (dateCreated != null)
			return df.format(dateCreated);
		else
			return null;
	}

	@Transient
	@XmlAttribute(name="last-access")
	protected String getLastAccess()
	{
		DateFormat df = new SimpleDateFormat ("yyyy-MM-dd");
		if (lastAccess != null)
			return df.format(lastAccess);
		else
			return null;
	}

	@Transient
	@XmlAttribute(name="datetime")
	protected String getDateTime()
	{
		DateFormat df = new SimpleDateFormat ("yyyy-MM-dd HH:mm:ss");
		if (dateCreated != null)
			return df.format(dateCreated);
		else
			return null;
	}

	@ManyToOne
	@JoinColumn(name = "article_id")
	@XmlTransient
	public Article article;

	//NB! Both these fields are fix! Should be refactored
	@Transient
	@XmlTransient
	public String implicitValues;

	@Transient
	@XmlTransient
	public String taxonomy;

	@Transient
	@XmlElement(name="descriptors-available")
	protected boolean getDescriptorsAvailable()
	{
		return (calcDescriptors != null);
	}

	public ModelMapping getMappingById(Long id)
	{
		for (ModelMapping mapping : modelMappings) {
			if (mapping.id.equals(id))
				return mapping;
		}

		return null;
	}

	public ModelMapping getMapping(Property property)
	{
		for (ModelMapping mapping : modelMappings) {
			if (mapping.property.equals(property))
				return mapping;
		}

		return null;
	}

	@XmlElement
	protected Article getArticle()
	{
		if (!ThreadScope.get().controller.equals("article"))
			return article;

		return null;
	}



	public void resetFilteredSets()
	{
		filteredSets =  new TreeMap<String, Basket>();
	}

	/**
	 * Implicit values are converted to Double to be used as substitution for real values: only one for all quantitativeProperties
	 * Or we can map implicit values to one of class options, but only one per qualitative Property
	 * If we use implicit mapping, it should be specified for all qualitative properties
	 * @return
	 * @throws IOException
	 */
	public Double[] generateImplicitValues() throws IOException{

		Map<Long, Long> options = attachment.getObject().optionsMapping;

		logger.info("Processing implicit values: "+implicitValues);

		String names[] = implicitValues.split(";");

		List<String> implicitClass = new ArrayList<String>();
		Double implicitFloat = null;

		for(String n : names)try {
			Double val = Double.valueOf(n);
			if(implicitFloat != null)
				throw new IOException("Two float implicit values: " + val + " and " + implicitFloat + " were provided while only one is expected ");
			implicitFloat = val;
		}catch(NumberFormatException e){
			implicitClass.add(n);
		}

		Double opt[] = new Double[modelMappings.size()];
		String optnames[] = new String[modelMappings.size()];

		for (int m = 0;  m < modelMappings.size() ; m++) {
			ModelMapping modelMapping = modelMappings.get(m);
			if (modelMapping.property.isQualitative()) {
				for(String name : implicitClass){
					@SuppressWarnings("unchecked")
					List<Long> implicitOption = (List<Long>) Globals.session()
					.createQuery("select id from PropertyOption where property=:property and name=:name")
					.setParameter("property", modelMapping.property).setParameter("name", name).list();
					if(implicitOption.size() == 0) continue; // no mapping can be done
					if(options.containsKey(implicitOption.get(0))) {
						if(optnames[m] != null)
							throw new IOException("Two implicit class values: " + optnames[m] + " and " + name + " were specified while only one is expected per class");
						opt[m] = options.get(implicitOption.get(0)).doubleValue();
						optnames[m] = name;
						logger.info("Added mapping of implicit values " + implicitValues + " to: " + ((int)(Math.floor(m))) + " using " + name + " for property: " + modelMapping.property);
					}

				}
				if(opt[m] == null) throw new IOException("No mapping of implicit value: " + implicitValues + 
						" is provided for property:" + modelMapping.property.getName());
			}
			else
				opt[m] = implicitFloat; // could be also null
		}



		return opt;
	}



	public String hashOfSets(ProvidedConditions cond){
		String hash = trainingSet.lastModified.toString()+trainingSet.id;
		for(String set: filteredSets.keySet())
			hash += "-"+ set+"_"+filteredSets.get(set).cachedCount;
		if(cond != null && cond.hasConditions())
			for(ShortCondition c:cond.getConditions())
				hash += "__"+c.id;

		return hash;
	}

	public Basket getFilteredSet(String name)
	{
		Basket filteredSet = filteredSets.get(name);
		if (filteredSet == null)
		{
			long timeStarted = Calendar.getInstance().getTimeInMillis();
			filteredSet = new Basket();
			filteredSets.put(name, filteredSet);

			Basket source = trainingSet;
			if (name.startsWith(QSPRConstants.VALIDATION))
			{
				if (getValidationSets().isEmpty())
					return null;
				// Which of the validation sets to choose?
				String numStr = name.substring(10);
				int num = numStr.isEmpty() ? 0 : Integer.valueOf(numStr);
				source = getValidationSets().get(num);
			}

			if (source == null)
				return null;

			filteredSet.id = source.id;

			for (BasketEntry basketEntry : source.entries) 
			{

				boolean excludeBasketEntry = true;
				String psn = basketEntry.ep.predicate.shortName;

				if ("=".equals(psn) || "+-".equals(psn))
					excludeBasketEntry = false;
				else
					if (QSPRConstants.USE.equals(attachment.getObject().datahandling.greaterless) && (">".equals(psn) || "<".equals(psn) || ">=".equals(psn) || "<=".equals(psn) || ">>".equals(psn) || "<<".equals(psn)))
						excludeBasketEntry = false;
					else
						if (QSPRConstants.USE.equals(attachment.getObject().datahandling.intervals) && "-".equals(psn))
							excludeBasketEntry = false;
						else
							if (QSPRConstants.USE.equals(attachment.getObject().datahandling.approximateequals) && ("~".equals(psn) || "~=".equals(psn)))
								excludeBasketEntry = false;					

				if (excludeBasketEntry)
				{
					// Automatically put the records with sophisticated predicates to the excluded set / Midnighter
					microattachment.getObject().excludedBasketEntries.add(basketEntry.id);
				}

				boolean exclude = microattachment.getObject().excludedBasketEntries.contains(basketEntry.id) || (basketEntry.exclude);
				boolean excludedSet = QSPRConstants.EXCLUDED.equals(name) || "excluded-on-the-fly".equals(name);

				if (!name.startsWith(QSPRConstants.VALIDATION))
					if (!(excludedSet ^ !exclude))
						continue;

				boolean foundAppropriateMapping = false;
				if (this.modelMappings.size() == 0)
					foundAppropriateMapping = true; // No mappings, then its a descriptors-only model
				else
					for (ModelMapping modelMapping : this.modelMappings) 
					{
						if (modelMapping.matches(basketEntry.ep))
							foundAppropriateMapping = true;
					}

				if (foundAppropriateMapping)
					filteredSet.entries.add(basketEntry);
				else
					logger.info("Compound "+basketEntry.ep.id+ " was excluded");
			}
			logger.info("Getting filtered set in " + (Calendar.getInstance().getTimeInMillis() - timeStarted) + "ms.");
		}

		filteredSet.cachedCount = (long)filteredSet.entries.size();

		if (filteredSet.entries.size() == 0) {
			if(QSPRConstants.TRAINING.equals(name))
				throw new UserFriendlyException("Training set is empty: check that you still have records in them (have you deleted or excluded all of them?)");
			return null;
		}

		return filteredSet;
	}

	@XmlElement(name = "last-modified")
	protected String getLastModified() // XML only
	{
		if (published)
			return TimeUtils.ago(lastModification != null ? lastModification : dateCreated);
		else
			return "";
	}


	public void storePendingModel() throws IOException, ClassNotFoundException, Exception
	{
		logger.info("Storing the model " + name + " id " + id);

		if (taskId != null && !isStatisticsCalculated)
			fetchCalculatedModel();

		session = Globals.userSession();
		dateCreated = new Timestamp(Calendar.getInstance().getTimeInMillis());
		Integer _taskId = taskId;
		taskId = null;
		Globals.session().saveOrUpdate(this);

		if (_taskId != null)
		{
			PendingTask pTask = PendingTask.getByTaskId(_taskId);
			if (pTask != null)
				Globals.session().delete(pTask);
		}
	}

	/**
	 * Fetch the calculated model from Metaserver
	 */
	public void fetchCalculatedModel() throws IOException, ClassNotFoundException, Exception {
		if (taskId == null)
			return;

		// Check if the model has already been fetched
		if (isStatisticsCalculated)
			return;
		//HERE!

		logger.info("Fetching a finished model from metaserver...");

		// Model is ready, but statistics was not calculated (the case of "fetch later" or "recalculate")
		ModelProcessor processor = ModelFactory.getProcessor(template);
		processor.model = this;
		CalculationClient client = new CalculationClient();
		client.setDeepSleepTime(1);
		client.setTolerateMetaserverDown();
		processor.onTaskReceived(client.getTask(taskId));
		processor.saveModel();

		logger.info("Model has been fetched successfully");
	}

	public void recalculateModelSize()
	{
		if (id != null)
		{
			size =  (attachment == null ?0 : attachment.getDataLength()) + 
					(readyModelAttachment == null ? 0: readyModelAttachment.getDataLength()) +
					(calcDescriptors == null ? 0: calcDescriptors.getDataLength()) + 
					(microattachment == null? 0: microattachment.getDataLength());
			logger.info(String.format("Recalculated size for model (ID=%d) is %d bytes", id, size));
		}
	}

	public void updateDescription() throws JAXBException
	{
		description = "";

		size = 0l;

		if (attachment != null)
		{
			// 1) A textual description of the model
			description += attachment.getObject().toString();
			if (readyModelAttachment != null)
				description += "\n-\n" + readyModelAttachment.getObject().toString();

			// 2) An XML file with model configuration
			// Do not marshall model data

			StringWriter writer = new StringWriter();
			Marshaller marshaller = Globals.jaxbContext.createMarshaller();
			marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT,
					new Boolean(true));
			marshaller.marshal(attachment.getObject(), writer);
			configurationXml = writer.toString()
					// a patch: in case of feature nets, do not store the sub-models in the XML, otherwise it can become huge
					.replaceAll("nodesConfiguration[^\\uFFFF]*nodesConfiguration", "nodesConfiguration>(skipped)</nodesConfiguration") 
					.replace(" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"", "");
		}

		logger.info("Model description: \n" + description);
	}

	public void addModelMapping(ModelMapping modelMapping)
	{
		this.modelMappings.add(modelMapping);
		modelMapping.model = this;
	}

	@XmlElement(name = "owner")
	public boolean getOwnerPrivileges()
	{
		if (Globals.userSession() == null)
			return false;
		if (session == null)
			return true;
		User user = Globals.userSession().user;
		return session == Globals.userSession() || (session.user != null && session.user.equals(user));
	}

	public static void addAuthRestrictions(Criteria c, boolean showGroupModels, boolean showUnapporvedModels)
	{
		c.createAlias("session", "sess");
		Disjunction authCriteria = Restrictions.disjunction();
		authCriteria.add(Restrictions.eq("sess.id", Globals.userSession().id));
		Conjunction publishedAndApproved = Restrictions.conjunction();
		publishedAndApproved.add(Restrictions.eq("published", new Boolean(true)));
		publishedAndApproved.add(Restrictions.isNull("taskId"));
		if (!showUnapporvedModels) 
			publishedAndApproved.add(Restrictions.eq("approved", Boolean.TRUE));
		authCriteria.add(publishedAndApproved);
		if (Globals.userSession().user != null)
			if (Globals.userSession().user.group != null && showGroupModels)
			{
				c.createAlias("sess.user", "u", Criteria.LEFT_JOIN);
				authCriteria.add(Restrictions.eq("u.group", Globals.userSession().user.group));
			}
			else
				authCriteria.add(Restrictions.eq("sess.user", Globals.userSession().user));
		c.add(authCriteria);
	}


	// Interactive DM calculation stuff. Wolfram
	public static Map<Long, Thread> dmThreads = new HashMap<Long, Thread>();
	@XmlElement(name = "is-dm-calculated")
	protected boolean isDMCalculated()
	{
		return id != null && dmThreads.get(id) != null && dmThreads.get(id).isAlive();
	}

	@Override
	public boolean equals(Object obj)
	{
		Model mm = (Model)obj; 
		if (this.id != null)
			return mm != null && this.id.equals(mm.id);
		else
			return this == obj;
	}

	@Override
	public int hashCode()
	{
		if (this.id != null)
			return this.id.hashCode();
		else
			return super.hashCode();
	}

	public Object getModelData(boolean recalculated)
	{
		if (recalculated && (readyModelAttachment != null && readyModelAttachment.getObject().modelDataRecalculated != null)) 
			return readyModelAttachment.getObject().modelDataRecalculated;

		return readyModelAttachment == null ? null: readyModelAttachment.getObject().modelData;
	}

	public void setModelData(Object object, boolean recalculated)
	{
		if (readyModelAttachment == null)
			readyModelAttachment = new Attachment<ReadyModelAttachment>(new ReadyModelAttachment(), AttachmentType.MARSHALABLE, AttachmentSource.Model);
		if (!recalculated)
		{
			readyModelAttachment.getObject().modelData = object;
		}
		else
			if (object != null && !object.equals(readyModelAttachment.getObject().modelData) && recalculation)
				readyModelAttachment.getObject().modelDataRecalculated = object;
	}

	public void createMappings(Map<Property, Unit> units) 
	{
		BasicRecordMapper brm = new BasicRecordMapper(trainingSet);
		List<Property> propList = trainingSet.getProperty();
		Unit prevUnit = modelMappings.size() > 0 ? modelMappings.get(0).unit : null;

		modelMappings.clear();
		for (Property property : propList) 
		{
			ModelMapping modelMapping = new ModelMapping();
			modelMapping.property = property;
			if (units != null && units.get(modelMapping.property) != null)
				modelMapping.unit = units.get(modelMapping.property);
			else if (prevUnit != null)
				modelMapping.unit = prevUnit;
			else
				modelMapping.unit = property.defaultUnit;
			modelMapping._class = brm.getClass(modelMapping.property);
			addModelMapping(modelMapping);
		}
	}

	public Privileges getPrivileges(HttpServletRequest request)
	{
		Privileges privileges = new Privileges("model");
		User user = Globals.userSession().user;

		// Only the owner of the a model can edit it
		privileges.canEdit = session == Globals.userSession()
				|| (session.user != null && session.user.equals(user));

		// Groups members can view the model
		boolean groupPrivileges = AccessChecker.isFromGroup(session.user, user);

		// For published baskets we allow to see all models
		boolean publishedPrivileges = published || (trainingSet.user != null && trainingSet.user.isPublisher());

		// Any user can view this model, if he knows its public identifier
		boolean retrievedByPublicId = request.getParameter("public_id") != null
				&& request.getParameter("public_id").equals("" + publicId);

		retrievedByPublicId |= Globals.userSession().visitedModelsIds.contains(publicId);

		privileges.canView = privileges.canEdit || groupPrivileges || retrievedByPublicId || publishedPrivileges || 
				Globals.isOCHEMDeveloper();

		return privileges;
	}

	public void replaceOriginalWithRecalculated()
	{
		setModelData(getModelData(true), false);
		attachment.updateObject();
		readyModelAttachment.updateObject();
		for (ModelMapping mm : modelMappings)
		{
			mm.statisticsOriginal.setObject(mm.statisticsRecalculated.getObject());
			Globals.session().saveOrUpdate(mm);
		}
		Globals.session().saveOrUpdate(this);
	}

	/*
	 * If setId is null, delete all the validation sets for this model
	 */
	public void deleteValidationSet(Long setId) {
		boolean deleted = false;
		Iterator<Basket> iter = validationSets.iterator();
		while (iter.hasNext())
			if (iter.next().id.equals(setId) || setId == null){
				iter.remove();
				deleted = true; 
			}
		if (validationSet != null && (validationSet.id.equals(setId) || setId == null))
		{validationSet = null;deleted = true;}

		if(!deleted)return; // skipping other steps in case nothing to be updated

		Globals.session().saveOrUpdate(this);
		for (ModelMapping mm : modelMappings)
		{
			((ModelStatistics)mm.statisticsOriginal.getObject()).deleteValidationSet(setId);
			((ModelStatistics)mm.statisticsRecalculated.getObject()).deleteValidationSet(setId);
			mm.statisticsOriginal.updateObject();
			mm.statisticsRecalculated.updateObject();
			Globals.session().saveOrUpdate(mm);
		}

		logger.info("sets were deleted");
	}

	public List<String> getStoredColumnsNames() {
		List<String> columns = new ArrayList<String>();
		if(columnList == null)return columns;
		String[] cols = columnList.split(",");
		for (String col : cols)
			if (!col.startsWith(QSPRConstants.INDIVIDUAL_PREDICTIONS)) //We do		return null;
				columns.add(col);
		return columns;
	}

	public void setColumnsNames(List<String> columns) {
		String columnCandidates =  StringUtils.join(columns, ",");

		if (columnCandidates != null && !"".equals(columnCandidates))
		{
			columnList = columnCandidates;
			Globals.session().saveOrUpdate(this);
		}		
	}

	public void sanitize() throws IOException {
		if(!(attachment.getObject().configuration instanceof ConsensusModelConfiguration))
			calcDescriptors = null;
		for (ModelMapping mm : modelMappings)
		{
			for(int set = 0; set <2;set++) {
				ModelStatistics msAnalysed = set == 0? (ModelStatistics) mm.statisticsOriginal.getObject() : (ModelStatistics) mm.statisticsRecalculated.getObject();
				for (int i = 0 ;  msAnalysed != null && msAnalysed.sets != null && i < msAnalysed.sets.size() ; i++ ){
					for(PointStatistics ps: msAnalysed.sets.get(i).points){
						ps.id = Integer.MAX_VALUE;
						ps.moleculeId = 0l;	
					} 
				}
			}

			mm.statisticsOriginal.updateObject();
			mm.statisticsRecalculated.updateObject();
			Globals.session().saveOrUpdate(mm);
		}
		Globals.session().saveOrUpdate(this);
	}


	public void publishModelPredictions(Long articleId) throws IOException {

		logger.info("Saving data records for model: " + name);
		String messages = "";

		List<Long> ids = retrieveRecordIds();

		List<Long> toUpdate = new ArrayList<Long>();

		for(Long epId: ids)try{
			if(epId == Integer.MAX_VALUE || epId == null) {messages += "\nno molecule id for mol = " + epId; continue;}
			ExperimentalProperty ep = (ExperimentalProperty) Globals.session().get(ExperimentalProperty.class, epId); 
			if(ep == null) {messages += "\nrecord was deleted R" + epId; continue;	}	
			if(ep.error != null) {messages += "record R" + ep.id + " is error. It can't be published."; continue;}
			if(ep.article.id.longValue() == articleId.longValue() && ep.owner.id == QSPRConstants.PUBLISHER_ID)continue;
			toUpdate.add(epId);
		}catch(Exception e){
			messages += "\n" + e.getMessage();
		}

		if(messages.length() > 0)
			throw new IOException(errorModelChanged + messages);

		// publishing
		int SIZE = QSPRConstants.MAXRECORDS_PER_MYSQL_UPDATE;

		if(articleId != null)
			do{
				List<Long>  update = toUpdate.subList(0, SIZE < toUpdate.size() ? SIZE : toUpdate.size());
				ExperimentalProperty.publishRecords(update, articleId);
				update.clear();
			}while(toUpdate.size()>0);

	}

	/**
	 *  Gets DataTable with prediction results for the given set
	 *  return null if the set does not exist
	 *  BUT CURRENTLY does not work for models with 2 or more properties
	 * @return
	 */

	private List<Long>  retrieveRecordIds(){

		List<Long> ids = new ArrayList<Long>();

		for(ModelMapping mapping: modelMappings) {		
			List<SetStatistics> sets = ((ModelStatistics) mapping.statisticsRecalculated.getObject()).sets;
			for(SetStatistics ss: sets) 		
				for(PointStatistics ps: ss.points) // for all datapoints
					ids.add(ps.id);
		}

		return ids;
	}

	public boolean restoreModifiedBasketEntries(boolean recalculated)
	{
		// Put the modified points back on their places
		boolean changesDone = false;

		for (ModelMapping mm : modelMappings)
		{
			ModelStatistics msReference = (ModelStatistics) (recalculated ? mm.statisticsRecalculated.getObject() : mm.statisticsOriginal.getObject());
			msReference.actualizeStatistics(mm);
			SetStatistics ssModified = msReference.getSetStatistics(QSPRConstants.MODIFIED_POINTS_ID);
			if (ssModified != null && !ssModified.points.isEmpty())
			{
				List<Long> modifiedIds = new ArrayList<Long>();
				for (PointStatistics ps : ssModified.points)
					modifiedIds.add(ps.id);
				mm.model.trainingSet.invertExclusion(modifiedIds);
				changesDone = true;
			}
		}
		return changesDone;
	}

	public boolean excludeDuplicates()
	{
		boolean changesDone = false;

		for (ModelMapping mm : modelMappings)
		{
			ModelStatistics msReference = (ModelStatistics) (mm.statisticsOriginal.getObject());
			msReference.actualizeStatistics(mm);

			String sets[] = {QSPRConstants.TRAINING,QSPRConstants.MODIFIED_POINTS_ID};

			Map<String,Double> minimalErrors = new HashMap<String,Double>();

			for(String set:sets) {
				SetStatistics ss = msReference.getSetStatistics(set);
				if(ss == null) continue;

				for (PointStatistics ps : ss.points) {
					if(ps.error != null) continue;

					String hashMol =  ps.getMoleculeHashWithValue();// hash by molecule and predicted value - to account for conditions
					double error = Math.abs(ps.predicted - ps.real);

					if(!minimalErrors.containsKey(hashMol) || minimalErrors.get(hashMol) > error)
						minimalErrors.put(hashMol, error); // first store an error
				}
			}

			List<Long> toExclude = new ArrayList<Long>();

			for(String set:sets) {
				SetStatistics ss = msReference.getSetStatistics(set);
				if(ss==null) continue;

				for (PointStatistics ps : ss.points) {
					if(ps.error != null) continue;

					String hashMol = ps.getMoleculeHashWithValue(); // hash by molecule and predicted value - to account for conditions
					double error = Math.abs(ps.predicted - ps.real);

					if(!minimalErrors.containsKey(hashMol) || error >  minimalErrors.get(hashMol)) { // was already selected; we keep only one record
						toExclude.add(ps.id);
						continue;
					}

					minimalErrors.remove(hashMol); // this molecule will be part of the training set, only one record
				}
			}

			if(!toExclude.isEmpty()) {
				mm.model.trainingSet.excludeOrIncludeEntries(toExclude, true);
				changesDone = true;
			}
		}

		return changesDone;
	}

	//TODO: are both functions below identical?

	public void mergeTrainingAndValidationSet() {
		trainingSet = (Basket) Globals.session().merge(trainingSet);

		if (validationSet != null)
			validationSet = (Basket) Globals.session().merge(validationSet);
	}

	public void initTrainingAndTestSetEntries() {
		if (trainingSet.id != null)
			trainingSet = (Basket) Globals.session().get(Basket.class, trainingSet.id);

		if (validationSet != null && validationSet.id != null)
			validationSet = (Basket) Globals.session().get(Basket.class, validationSet.id);

		Hibernate.initialize(trainingSet.entries);
		if (validationSet != null)
			Hibernate.initialize(validationSet.entries);	
	}

	public boolean isStratifyValidationConfiguration(){
		if(attachment != null && attachment.getObject() != null && attachment.getObject().protocol != null && attachment.getObject().protocol.validationConfiguration !=null){
			ValidationConfiguration config = (ValidationConfiguration) attachment.getObject().protocol.validationConfiguration;
			return config.validationType != null && config.validationType == BaggingConfiguration.STRATIFIED;
		}
		return false;	
	}


	public boolean isStoredConsistently() throws IOException, ClassNotFoundException
	{
		for (ModelMapping mm : modelMappings)
			if (mm.statisticsOriginal != null && mm.statisticsRecalculated != null && 
			((Attachment<ModelStatistics>)mm.statisticsOriginal).id != ((Attachment<ModelStatistics>)mm.statisticsRecalculated).id)
				return false;

		return true;
	}

	/**
	 * Identifies inconsistencies in configuration and conditions
	 * Requires to enumerate conditions to understand whether we have correct correspondence
	 * Correct some small inconsistencies
	 */

	public void isConsistent(){

		try {

			if(!(attachment.getObject().configuration instanceof CDSConfiguration))return;

			CDSConfiguration cds = (CDSConfiguration)attachment.getObject().configuration;

			ModelAbstractConfiguration conf =  cds.modelConfiguration;

			if(cds.hasConditions() && !conf.isSupportConditions())
				throw new UserFriendlyException("Conditions currently cannot be used with: " + conf.getDefaultName());

			if(isStratifyValidationConfiguration() && !conf.areClassificationData())
				throw new UserFriendlyException("Stratified protocols can be used only for classification tasks. Select another validation protocol.");

			if(!cds.hasDescriptors() && !cds.hasConditions() && !(conf instanceof NoDescriptors))
				throw new UserFriendlyException("No descriptors were selected for : " + conf.getDefaultName() + ". Calculations for this method without descriptors are not possible.");

			if(!conf.isSupportRegression() &&  !conf.areClassificationData()) 
				throw new UserFriendlyException("The method \"" + conf.getInformativeName() + "\" only supports classification. It should not be used with regression data. Use other methods.");

		}catch(UserFriendlyException e) {
			detailedStatus = e.getMessage();
			status = QSPRConstants.ERROR_STATUS;
			throw e;			
		}
	}

	static public Model unproxyModel(Model model) {
		if (model.id != null){
			model = initialiseModel(model);
			model = Globals.initializeAndUnproxy(model);
		}
		return model;
	}

	/**
	 * To initialize lazy modelMappings
	 */
	static public Model initialiseModel(Model model){
		if(model.id != null && !Globals.session().contains(model))
			model = (Model) Globals.session().get(Model.class, model.id);
		return model;
	}

	@Override
	public String toString() {
		return " id:" + publicId + " name: \"" + name + "\" ";
	}

	/**
	 * @returns id for mapping of an option
	 */

	public Long getMappedOption(long iDInDatabase) {
		return attachment.getObject().optionsMapping.get(iDInDatabase);
	}

	public boolean hasMappedOptions() {
		return attachment.getObject().optionsMapping != null && !attachment.getObject().optionsMapping.isEmpty();
	}

	public String getTrainingSetHash() throws NoSuchAlgorithmException {

		String hash ="";
		for(ModelMapping map: modelMappings) {
			if(map == null || map.statisticsRecalculated == null || map.statisticsRecalculated.getObject() == null)continue;
			ModelStatistics stat = (ModelStatistics)map.statisticsRecalculated.getObject();
			SetStatistics statis = stat.getSetStatistics(QSPRConstants.TRAINING);
			if(statis == null)continue;
			hash += "_" + 	statis.getHash();
		}
		return OCHEMUtils.getMD5(hash);
	}


	public List<Property> getNotPublishedProperties() {
		List<Property> propList = trainingSet.getProperty();
		List<Property> unpublished = new ArrayList<Property>();
		for(Property p: propList)
			if(!p.isPublished())
				unpublished.add(p);
		return unpublished;
	}

	public String canBeMixtures() {
		ValidationConfiguration valid = attachment.getObject().protocol.validationConfiguration;
		if(valid == null || valid.mixtureValidation != null && valid.mixtureValidation == MixtureValidation.RECORD) return "false";
		return attachment.getObject().standartization.desaltWith == null?"true":"";
	}

	public boolean isCompatibleModelAndDescriptors() {

		if(OCHEMConfiguration.mirror && OCHEMConfiguration.autoLoginUser == null)
			throw new UserFriendlyException("Submission of models using mirror web site is not allowed");

		if(attachment.getObject().configuration instanceof ConsensusModelConfiguration) return true;

		CDSConfiguration cdsconf = (CDSConfiguration)  attachment.getObject().configuration;

		ModelAbstractConfiguration modelcfg = (ModelAbstractConfiguration) cdsconf.modelConfiguration;

		if(modelcfg instanceof SupportsOneOutputOnly && modelcfg.containMultiClassData()) return false;

		int count = trainingSet.entries.size();

		if(count<modelcfg.requireMinimumRecords() && OCHEMConfiguration.autoLoginUser == null)
			throw new UserFriendlyException(modelcfg.getInformativeName() + " requires minumum " +modelcfg.requireMinimumRecords()+ " data points but your training set has only " + count + " data points. Task cannot be started.");

		ValidationConfiguration valid = attachment.getObject().protocol != null?attachment.getObject().protocol.validationConfiguration:null;

		if(valid != null) {

			Map<Property,Integer> properties = new HashMap<Property,Integer>();

			for(BasketEntry e: trainingSet.entries)
				if(properties.containsKey(e.ep.property))
					properties.put(e.ep.property, properties.get(e.ep.property)+1);
				else
					properties.put(e.ep.property,1);

			Property min = null;
			int countOne = count;
			for(Property p:properties.keySet())
				if(properties.get(p)<=countOne) {
					countOne = properties.get(p);
					min = p;
				}

			if(!valid.compatibleWithDesalt(attachment.getObject().standartization.desaltWith != null))
				throw new UserFriendlyException("Validation " + valid.getInformativeName() + " is not compatible with \"Remove salts\" option in \"Preprocessing of molecules\"");

			if(count < valid.requireMinimumRecords() && OCHEMConfiguration.autoLoginUser == null)
				throw new UserFriendlyException(valid.getInformativeName() + " protocol requires minumum " + valid.requireMinimumRecords()+ " data points but your training set has only " + count + " data points. Task cannot be started.");

			if(!Globals.isOCHEMDeveloper() && !WorkflowNodeServer.isRunningTest()){
				if(valid.ensembleSize > QSPRConstants.MAXBAGGING)
					throw new UserFriendlyException("Maximum bag size per model is: " + QSPRConstants.MAXBAGGING);
				if(valid instanceof BaggingConfiguration && modelcfg.isLarge())
					throw new UserFriendlyException("This model " + modelcfg.getInformativeName() + "  can not be used with the bagging since it will produce too large models.");
				if(valid.ensembleSize > QSPRConstants.MAXCV && valid instanceof CrossValidationConfiguration)
					throw new UserFriendlyException("Maximum CV size is " + QSPRConstants.MAXCV);
				if(valid.ensembleSize < QSPRConstants.MINCV)
					throw new UserFriendlyException("Minimum CV size is " + QSPRConstants.MINCV);
				if(valid.ensembleSize < QSPRConstants.MINBAGGING && valid instanceof BaggingConfiguration)
					throw new UserFriendlyException("Minimum Bagging size is " + QSPRConstants.MINBAGGING);
			}

			if(valid.ensembleSize >= countOne && valid instanceof CrossValidationConfiguration && !(modelcfg instanceof SupportsInversions))
				throw new UserFriendlyException("The number of records (" + countOne + ") for property " + min + " in your basket is less than number of CV folds. Task cannot be started.");

		}

		if(modelcfg instanceof NoDescriptors) { // only method supporting augmentation should be used
			NoDescriptors cfg = (NoDescriptors)modelcfg;
			if(!cfg.isSupportAugmentation() && (cfg.getAugementTraining()>1 || cfg.getAugmentApply()>1 || cfg.getBalanceData()))return false;
		}else
			if(valid != null && valid.mixtureValidation != null && 
			valid.mixtureValidation != MixtureValidation.RECORD && 
			cdsconf.descriptors.mixtures == null) //always fraction is added
				throw new UserFriendlyException("Validation " + valid.getInformativeName() + " should be used with mixture descriptors only. Configuration: " + cdsconf.toString());

		return cdsconf.isCompatibleDescriptorsAndMethod(modelcfg);
	}

}


