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
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Date;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import javax.persistence.CascadeType;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.FetchType;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.JoinTable;
import javax.persistence.ManyToMany;
import javax.persistence.ManyToOne;
import javax.persistence.OneToMany;
import javax.persistence.OrderBy;
import javax.persistence.Transient;
import javax.xml.bind.JAXBException;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlElementWrapper;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.FlushMode;
import org.hibernate.annotations.Cascade;
import org.hibernate.annotations.Filter;
import org.hibernate.annotations.FilterDef;
import org.hibernate.annotations.Filters;
import org.hibernate.criterion.Disjunction;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.MarshallingOption;
import qspr.ThreadScope;
import qspr.annotations.Loggable;
import qspr.business.PendingTaskFilter;
import qspr.entities.Attachment.AttachmentType;
import qspr.util.AccessChecker;
import qspr.util.DynaWrap;
import qspr.util.ISBNUtility;
import qspr.util.PubMedUtility;
import qspr.util.UploadContext;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.mmpa.MMPAnnotationService;
import com.eadmet.mmpa.domain.MMPAnnotationSet;
import com.eadmet.utils.OCHEMUtils;

@FilterDef(name="no-experimental-data")

@Entity
@XmlRootElement(name = "article")
@Loggable
@SuppressWarnings("unchecked")
public class Article
{
	private static transient final Logger logger = LogManager.getLogger(Article.class);
	public static final String defaultDateFormat = "yyyy-MM-dd HH:mm:ss";

	@Id
	@Column(name = "article_id")
	@GeneratedValue
	@XmlAttribute
	public Long id;

	private static int PUBMED = 100000;

	private String title;

	/**
	 * Required for web form. Should not be removed!
	 */

	public String short_title;

	public String affiliation;

	@Loggable(name = "pubmed ID")
	public Long pmid;

	@Column(name = "abstract")
	@Loggable(name = "abstract")
	public String articleAbstract;

	@ManyToOne
	@JoinColumn(name="journal_id")
	@XmlElement
	public Journal journal;

	public String volume;

	@Column(length = 20)
	public String issue;

	@Column(name = "page_nums")
	@Loggable(name = "page numbers")
	public String pageNumbers;

	@Column(name = "pub_date")
	@XmlTransient
	@Loggable(name = "publication date")
	public Date publicationDate;

	@XmlElement
	private String url;

	public String doi;

	@ManyToOne
	@JoinColumn(name = "modifier_id")
	@XmlTransient
	@Loggable(name = "modifier")
	public User owner;

	@ManyToOne
	@JoinColumn(name = "introducer_id")
	@XmlTransient
	public User introducer;

	@Column(name = "media_type")
	public String mediaType;

	public String publisher;

	public String isbn;

	public String isbn13;

	@Column(name = "is_chapter")
	public Boolean isChapter; 

	@ManyToOne
	@JoinColumn(name = "parent_id")
	public Article parent;

	@Column(name = "time_modified")
	@XmlTransient
	@Loggable (exclude = true)
	public Timestamp time;

	@Column(name = "time_created")
	@XmlTransient
	@Loggable (exclude = true)
	public Timestamp timeCreated;

	@Column
	public String comment;


	@ManyToMany
	(
			targetEntity = Author.class,
			cascade={CascadeType.PERSIST, CascadeType.MERGE}
			//fetch = FetchType.EAGER			
			)
	@JoinTable
	(
			name="ArticleAuthor",
			joinColumns={@JoinColumn(name="article_id")},
			inverseJoinColumns={@JoinColumn(name="author_id")}
			)
	@org.hibernate.annotations.IndexColumn(name="number")
	@XmlElementWrapper(name = "authors")
	@XmlElement(name = "author")
	public List<Author> authors = new ArrayList<Author>();

	@ManyToMany
	(
			targetEntity = Property.class,
			cascade={CascadeType.PERSIST, CascadeType.MERGE}, 
			fetch = FetchType.LAZY
			)
	@JoinTable
	(
			name="ArticleProperty",
			joinColumns={@JoinColumn(name="article_id")},
			inverseJoinColumns={@JoinColumn(name="property_id")}
			)
	@XmlTransient
	@Loggable(exclude = true)
	// Auxiliary superfluous data, showing what properties are present in given article
	// Needed mostly for fast tagination of the article
	// Midnighter
	public Set<Property> properties = new HashSet<Property>();

	@ManyToMany
	(
			targetEntity = Tag.class,
			cascade={CascadeType.PERSIST, CascadeType.MERGE},
			fetch = FetchType.LAZY
			)
	@JoinTable
	(
			name="ArticleTag",
			joinColumns={@JoinColumn(name="article_id")},
			inverseJoinColumns={@JoinColumn(name="tag_id")}
			)
	@XmlTransient
	public Set<Tag> tags = new HashSet<Tag>();

	@OneToMany(mappedBy = "article", cascade = CascadeType.ALL, fetch = FetchType.LAZY)
	@Filters({
		@Filter(name="userFilter", condition="user_id = :userId")
	})
	@Cascade(org.hibernate.annotations.CascadeType.DELETE_ORPHAN)
	@XmlTransient
	@OrderBy
	@Loggable(name = "files")
	public Set<ArticleUserPdf> pdfs = new HashSet<ArticleUserPdf>();

	@OneToMany (fetch = FetchType.LAZY, mappedBy = "article")
	@XmlTransient
	public List<ExperimentalProperty> experimentalProperties;

	@OneToMany (fetch = FetchType.LAZY, mappedBy = "article")
	@XmlTransient
	public List<Model> models;

	@Transient
	@XmlElement
	public Article dublicate;

	@XmlAttribute(name="basketCount")
	@Transient
	public Long count;

	@XmlElement(name = "publication-date")
	XmlDate getPublicationDate()
	{
		return new XmlDate(publicationDate);
	}

	public Article()
	{
	}

	public Article(String _title)
	{
		setTitle(_title);
	}

	@XmlAttribute@Column(unique = true)
	public Article setTitle(String _title)
	{
		this.title = _title;
		this.short_title = Article.shortTitle(_title);
		return this;
	}

	public void addComment(String message){
		message = message == null ? "" : message;
		comment = comment != null ? comment + message : message;
		if(comment.length()>=255)comment = comment.substring(0, 255);
	}

	public void setComment(String message){
		comment = "";
		addComment(message);
	}

	public String getTitle()
	{
		return this.title;
	}

	public String getUrl()
	{
		return url;
	}

	public void setLink(String link)
	{
		if(link != null && !link.equals("") && !link.startsWith("http://"))
			this.url = "http://"+link;
		else
			this.url = link;
	}

	private Long getRecordSize(Integer dummy){
		if (ThreadScope.get().controller.equals("article"))
		{
			if((this.id != null) && (this.id > 0)){
				Criteria criteria = Globals.session().createCriteria(ExperimentalProperty.class)
						.add(Restrictions.isNull("deleted"));

				ExperimentalProperty.addAccessRestrictions(criteria, Globals.RIGHTS_FREELY_AVAILABLE, null, false, true);
				criteria.add(Restrictions.eq("article.id", this.id));
				criteria.createAlias("property", "p");
				if (dummy != null)
				{
					if (dummy == 1)
						criteria.add(Restrictions.eq("p.name", QSPRConstants.DUMMY));
					else
						criteria.add(Restrictions.ne("p.name", QSPRConstants.DUMMY));
				}
				criteria.setProjection(Projections.rowCount());
				return (Long)criteria.list().get(0);
			}
		}
		return 0L;
	}

	@XmlAttribute(name="property-record")
	public Long getRecordSize()
	{
		return getRecordSize(null);
	}

	@XmlAttribute(name="property-record-nondummy")
	public Long getRecordSizeNondummy()
	{
		return getRecordSize(0);
	}

	@XmlAttribute(name="property-record-dummy")
	public Long getRecordSizeDummy()
	{
		return getRecordSize(1);
	}

	@XmlAttribute(name="structural-alerts")
	public Long getAlertsSize(){
		if (ThreadScope.get().controller.equals("article"))
		{
			if((this.id != null) && (this.id > 0)){
				Criteria criteria = Globals.session().createCriteria(SubstructureAlert.class);
				AccessChecker.addAccessRestrictions(criteria, Globals.RIGHTS_FREELY_AVAILABLE, Globals.userSession().user, null, true);
				criteria.add(Restrictions.eq("article.id", this.id));
				criteria.setProjection(Projections.rowCount());
				return (Long)criteria.list().get(0);
			}
		}
		return 0L;
	}

	@XmlAttribute(name="model-list")
	protected Long getModelList(){
		if (ThreadScope.get().controller.equals("article") && this.id != null && this.id > 0)
		{
			Criteria criteria = Globals.session().createCriteria(Model.class).add(Restrictions.eq("article", this)).add(Restrictions.eq("published", true));
			Model.addAuthRestrictions(criteria, true, true);
			criteria.setProjection(Projections.rowCount());
			return (Long) criteria.list().get(0);
		}
		return null;
	}

	@XmlElement(name="models")
	protected List<Model> getModels() throws JAXBException{
		if (ThreadScope.get().controller.equals("article") && this.id != null && this.id > 0)
		{
			Criteria criteria = Globals.session().createCriteria(Model.class).add(Restrictions.eq("article", this)).add(Restrictions.eq("published", true));
			Model.addAuthRestrictions(criteria, true, true);
			List<Model> models = criteria.list();

			Set<Model> depdupeModelIds = new LinkedHashSet<Model>(models);
			models.clear();
			models.addAll(depdupeModelIds);

			return models;
		}
		return null;
	}

	@XmlElement(name = "mmp-set")
	protected List<MMPAnnotationSet> getMMPSets() {
		if (Globals.getMarshallingOption(MarshallingOption.ARTICLE_PENDING_TASKS))
		{
			List<MMPAnnotationSet> sets = Globals.session().createCriteria(MMPAnnotationSet.class).add(Restrictions.eq("article", this)).add(Restrictions.eq("published", true)).list();
			for (MMPAnnotationSet set : sets) {
				MMPAnnotationService.calculateAnnotationCounts(set);
			}
			return sets;
		}
		else
			return null;
	}

	@XmlAttribute(name="tasks-count")
	public Long getTasksCount() throws Exception{
		if (ThreadScope.get().controller.equals("article") && this.id != null && this.id > 0)
		{
			PendingTaskFilter filter = new PendingTaskFilter();
			filter.ownTasksOnly = false;
			filter.articleID = this.id;
			Criteria criteria = filter.createCriteria();
			criteria.setProjection(Projections.rowCount());
			return (Long) criteria.list().get(0);
		}
		return null;
	}

	@XmlElement(name="pending-task")
	protected List<PendingTask> getPendingTasks() throws Exception{
		if (Globals.getMarshallingOption(MarshallingOption.ARTICLE_PENDING_TASKS))
		{
			PendingTaskFilter filter = new PendingTaskFilter();
			filter.ownTasksOnly = false;
			filter.articleID = this.id;
			Criteria criteria = filter.createCriteria();
			return criteria.list();
		}
		return null;
	}

	@XmlAttribute(name = "pdf-available")
	String getPdfAvailable()
	{
		if (id != null)
			if (getPdf() != null)
				return "true";
		return null;
	}

	private List<ArticleUserPdf> getAttachedFiles(int type)
	{
		if (Globals.userSession() == null || Globals.userSession().user == null)
			return null;

		if (mediaType == null)
			return null;

		Criteria criteria = Globals.session().createCriteria(ArticleUserPdf.class)
				.add(Restrictions.eq("article", this))
				.add(Restrictions.eq("type", type));

		Disjunction disj = Restrictions.disjunction();
		disj.add(Restrictions.eq("user", Globals.userSession().user));
		if (Globals.userSession().user.group != null)
		{
			criteria.createAlias("user", "u");
			disj.add(Restrictions.eq("u.group", Globals.userSession().user.group));
		}

		criteria.add(disj);

		return criteria.list();
	}

	public ArticleUserPdf getPdf()
	{
		if(Globals.userSession() != null && Globals.userSession().user != null)
		{
			List<ArticleUserPdf> pdfs = getAttachedFiles(ArticleUserPdf.PDF);
			if(pdfs != null && pdfs.size() > 0)
				return pdfs.get(0);
		}

		return null;
	}

	@XmlAttribute(name = "batch-file")
	String getBatchFile()
	{
		if(Globals.userSession() != null && Globals.userSession().user != null && pdfs != null)
		{
			for (ArticleUserPdf pdf : pdfs) {
				if(pdf.type == ArticleUserPdf.EXCEL || pdf.type == ArticleUserPdf.SDF)
					return "true";
			}
		}
		return null;
	}

	@XmlAttribute(name = "pdf-exist")
	public String isPdfexist()
	{
		if(this.id != null && this.id > 0 && ThreadScope.get().controller.equals("article"))
		{
			Criteria criteria = Globals.session().createCriteria(ArticleUserPdf.class)
					.add(Restrictions.eq("article", this))
					.add(Restrictions.eq("type", ArticleUserPdf.PDF));
			if(criteria.list().size() > 0)
				return "true";
		}

		return null;
	}

	@XmlElement(name = "owner")
	public String getOwner()
	{
		if (owner != null)
			return owner.login;
		else
			return null;
	}

	@XmlElement(name = "introducer")
	public String getIntroducer()
	{
		if (introducer != null)
			return introducer.login;
		else
			return null;
	}

	public void addAuthor(Author author)
	{
		if (!authors.contains(author))
		{
			authors.add(author);
			author.articles.add(this);			
		}
	}

	public void addExperimentalProperty(ExperimentalProperty expProperty)
	{
		//experimentalProperties.add(expProperty);
		expProperty.article = this;
	}

	public static Date parseDate(String value)
	{
		Date temp = null;
		String[] formats = {"d MMM y", "d MM y", "MMM d y", "MM d y", "y/MMM/d", "y/MM/d", "y MM d", "MM y","MMM y", "y"};
		for (int i=0; i<formats.length; i++)
		{
			try
			{
				temp = (new SimpleDateFormat(formats[i])).parse(value);
				break;
			} catch (Exception e)
			{	
			}
		}
		return temp;
	}

	public void setPublicationDate(String value)
	{
		this.publicationDate = parseDate(value);
	}

	public void setCurrentPublicationDate()
	{
		Calendar pub_date = Calendar.getInstance();
		this.publicationDate = pub_date.getTime();
	}


	public void doCheckRights() throws Exception
	{
		if (!(this.owner == null || this.owner.equals(Globals.userSession().user)))
			if (!(Globals.userSession().user != null && Globals.userSession().user.rank > this.owner.rank))
			{
				if (!(Globals.userSession().user != null && Globals.userSession().user.rank.equals(this.owner.rank)))
				{
					logger.info("Not permitted edit has been ignored");
					throw new UserFriendlyException("Not permitted edit has been ignored");
				}
			}
	}

	public static Article getArticle(String article, UploadContext context) throws Exception
	{
		article = OCHEMUtils.trim(article);

		Article art = null;
		try
		{
			art = context.articleCache.get(article);

			if (art != null)
				return art;

			if (article.equals(""))
				throw new Exception("Empty article string provided");

			if (article.startsWith("Q") || article.startsWith("A"))
			{
				String ids = article.replaceFirst("Q", "").replaceFirst("A", "");
				Long id = Long.valueOf(ids);

				List<Article> articles = Globals.session().createCriteria(Article.class).add(Restrictions.eq("id", id)).list();

				if (articles.size() > 0)
				{
					art = articles.get(0);
					return art;
				}

				throw new Exception("Article not found in OCHEM by internal ID: "+article);
			} else // expecting pubmed id contains only digits
				if (article.matches("(\\d)*"))
				{
					Long id = Long.valueOf(article);
					List<Article> articles = Globals.session().createCriteria(Article.class).add(Restrictions.eq("pmid", id)).list();

					if (articles.size() > 0)
					{
						art = articles.get(0);
						return art;
					}

					if (!context.allowArticlePubmedSearch)
						throw new UserFriendlyException("Article not found in OCHEM, and PubMed search is disabled");

					if (id < PUBMED)
						throw new UserFriendlyException("PubMed ID can not be less than " + PUBMED + ". Use a valid PubMed ID or fetch article manually and  use ARTICLEID.");

					Globals.restartAllTransactions(true);
					art = PubMedUtility.getArticleByPubmedId(id);
					Globals.restartAllTransactions(true);

					if (art != null)
					{
						Globals.session().saveOrUpdate(art);
						Globals.session().flush();  
						return art;
					}

					throw new Exception("Article not found neither in OCHEM nor in PubMed");

				} else
				{
					throw new Exception("Wrong value provided - can only fetch article by Internal ID or PubMed ID");
				}
		} catch (Exception e)
		{
			context.articleCache.put(article,e);
			throw e;
		} finally
		{
			if (art != null)
				context.articleCache.put(article, art);
		}
	}

	public static Article getDefaultArticle(User user, UploadContext context) throws Exception {
		DynaWrap extended = user.dynaWrapped();
		String firstName = user.isExtended() ? extended.getString("firstName") : user.login;
		String lastName = user.isExtended() ? extended.getString("lastName") : "";
		
		if(firstName.length()==0)firstName = user.login;
		
		String name = QSPRConstants.UNPUBLISHED_ARTICLE+firstName+" " + lastName;

		if (context.articleCache.get(name) != null)
			return context.articleCache.get(name);
		Article article = Article.getByTitle(name);

		if (article == null)
		{
			article = new Article();
			article.mediaType = "article";
			article.addAuthor(Author.get(firstName, lastName, ""));

			article.journal = (Journal) Globals.session().get(Journal.class, QSPRConstants.UNPUBLISHED_JOURNAL); // re-fetching the model

			if (article.journal == null)
				throw new UserFriendlyException("Unpublished journal with journal_id = 1 is absent. Add it to the Journal table!");

			article.setTitle(name);
			context.articleCache.put(name, article);
			Calendar calendar = Calendar.getInstance();
			article.time=article.timeCreated=new java.sql.Timestamp(calendar.getTime().getTime());
			article.publicationDate=calendar.getTime();
		}
		return article;
	}

	public static Article getArticle(String internalOrPubMedID) throws Exception
	{
		internalOrPubMedID = OCHEMUtils.trim(internalOrPubMedID);

		if (internalOrPubMedID.equals("")) 
		{
			return null;
		} 

		if (internalOrPubMedID.startsWith("Q") || internalOrPubMedID.startsWith("A"))
		{
			String ids = internalOrPubMedID.replaceFirst("Q", "").replaceFirst("A", "");
			Long id = Long.parseLong(ids);

			List<Article> articles = Globals.session().createCriteria(Article.class).add(Restrictions.eq("id", id)).list();

			if (articles.size() > 0)
			{
				return articles.get(0);
			}
		} else // expecting pubmed id contains only digits
			if (internalOrPubMedID.matches("(\\d)*"))
			{
				Long id = Long.parseLong(internalOrPubMedID);
				List<Article> articles = Globals.session().createCriteria(Article.class).add(Restrictions.eq("pmid", id)).list();

				if (articles.size() > 0)
				{
					return articles.get(0);
				}
				else
				{
					if (id < PUBMED)
						throw new UserFriendlyException("PubMed ID can not be less than " + PUBMED + ". Use a valid PubMed ID or fetch article manually and  use ARTICLEID.");
					Article article = PubMedUtility.getArticleByPubmedId(id);
					Globals.session().saveOrUpdate(article);
					Globals.session().flush();  
					return article;
				}

			} else
			{
				throw new UserFriendlyException("Error uploading an article: Wrong value provieded - can only fetch Article by internal ID or PubMed ID");
			}
		return null;
	}

	public static Article getById(long id)
	{
		Article article = (Article) Globals.session().createCriteria(Article.class)
				.add(Restrictions.eq("id", id)).uniqueResult();
		return article;
	}

	/**
	 * @param pmid pubmed ID
	 * @return article from ochem db with given pmid or null 
	 */
	public static Article getByPmId(long pmid)
	{
		//return getByPmId(pmid, new Article());
		Article article = (Article) Globals.session().createCriteria(Article.class)
				.add(Restrictions.eq("pmid", pmid)).uniqueResult(); 
		return article;
	}

	public static Article getByShortTitle(String shortitle)
	{
		Article article = (Article) Globals.session().createCriteria(Article.class)
				.add(Restrictions.eq("short_title", shortTitle(shortitle)))
				.uniqueResult();
		return article;
	}

	public static Article getByISBN(String isbn13) throws Exception
	{
		List<Article> articles = 
				Globals.session()
				.createCriteria(Article.class).add(Restrictions.or(Restrictions.eq("isbn", isbn13), Restrictions.eq("isbn13", isbn13)))
				.list();

		return articles.size() > 0 ? articles.get(0) : ISBNUtility.getArticleByISBN(isbn13, new Article());
	}

	public static Article getByTitle(String title)
	{
		List<Article> articles = Globals.session().createQuery("from Article where title=:title").setString("title", title.trim()).list();
		return articles.size() > 0 ? articles.get(0) : null;
	}

	public static Article getByTitle(String title, boolean createIfMissing)
	{
		Article a = getByTitle(title);
		if (a != null)
			return a;
		if (!createIfMissing)
			return null;
		a = new Article();
		a.setTitle(title);
		return a;
	}

	public void attachFile(byte[] data, Integer fileTypeFlag)
	{
		if (this.id != null)
		{
			List<ArticleUserPdf> attachments = getAttachedFiles(fileTypeFlag);
			for (ArticleUserPdf aFile : attachments)
			{
				Globals.session().delete(aFile);
				pdfs.remove(aFile);
			}
			Globals.session().flush();
		}

		// Only one file of each type is allowed per user/article
		ArticleUserPdf aup = new ArticleUserPdf();
		aup.article = this;
		aup.user = Globals.userSession().user;
		aup.attachedFile = Attachment.getAttachment(data, AttachmentSource.Article, AttachmentType.PDF);
		aup.type = fileTypeFlag;
		aup.article.pdfs.add(aup);

		if (this.id != null)
			Globals.session().saveOrUpdate(aup);
	}

	public void attachFile(Attachment<ArticleUserPdf> attachment, Integer fileTypeFlag)
	{
		if (this.id != null)
		{
			List<ArticleUserPdf> attachments = getAttachedFiles(fileTypeFlag);
			for (ArticleUserPdf aFile : attachments)
			{
				Globals.session().delete(aFile);
				pdfs.remove(aFile);
			}
		}

		// Only one file of each type is allowed per user/article
		ArticleUserPdf aup = new ArticleUserPdf();
		aup.article = this;
		aup.user = Globals.userSession().user;
		aup.attachedFile = attachment;
		aup.type = fileTypeFlag;
		aup.article.pdfs.add(aup);

		if (this.id != null)
			Globals.session().saveOrUpdate(aup);
	}

	public void removeAttachment(int fileTypeFlag)
	{
		List<ArticleUserPdf> aups = getAttachedFiles(fileTypeFlag);
		for (ArticleUserPdf articleUserPdf : aups)
			Globals.session().delete(articleUserPdf);
	}

	public void updateProperties()
	{
		List<Property> list = Globals.session()
				.createCriteria(Article.class)
				.add(Restrictions.eq("id", this.id))
				.createCriteria("experimentalProperties", "ep")
				.createCriteria("property")
				.setProjection(Projections.distinct(Projections.property("ep.property")))
				.list();

		properties.clear();
		properties.addAll(list);
	}

	public String toString()
	{
		return this.title;
	}

	public Article hasArticle() throws Exception
	{
		// Added to avoid session flush before query to database. 
		// Session flush may result in dublicate constraint violation.
		Globals.session().setFlushMode(FlushMode.MANUAL);

		try {
			DateFormat format = new SimpleDateFormat("yyyy-MM-dd hh:mm:ss");
			Calendar pub_date = Calendar.getInstance();
			pub_date.setTime(this.publicationDate);
			Date startDate = (Date) format.parse("" + pub_date.get(Calendar.YEAR) + "-01-01 00:00:00");
			Date endDate = (Date) format.parse("" + pub_date.get(Calendar.YEAR) + "-12-31 23:59:59");

			if(this.title == null)title="";

			Criteria criteria1 = Globals.session().createCriteria(Article.class)
					.add(Restrictions.eq("short_title", Article.shortTitle(this.title.trim())))
					.add(Restrictions.between("publicationDate", startDate, endDate));
			if (id != null)
			{
				criteria1.add(Restrictions.ne("id", this.id));
			}

			List<Article> article = (List<Article>)criteria1.list();
			if (article.size() > 0)
				return article.get(0);

			if(journal == null || journal.id == QSPRConstants.UNPUBLISHED_JOURNAL)return null; // we do not check Unpublished; can have multiple entries with the same name, etc.

			else if (this.mediaType.equals("article"))
			{
				Criteria criteria2 = Globals.session().createCriteria(Article.class).add(Restrictions.between("publicationDate", startDate, endDate));
				Calendar calendar = Calendar.getInstance();
				calendar.setTime(publicationDate);

				if (calendar.get(Calendar.YEAR) == 0)
					return null;
				if (volume != null)
					criteria2.add(Restrictions.eq("volume", volume.trim()));
				if (id != null)
					criteria2.add(Restrictions.ne("id", id));
				if (pageNumbers != null)
					criteria2.add(Restrictions.like("pageNumbers", this.pageNumbers.split("-")[0].trim() + "-%"));
				if (journal != null)
					criteria2.add(Restrictions.eq("journal", journal));
				article = criteria2.list();
				if (article.size() > 0)
					return article.get(0);
			}
			return null;
		} finally
		{
			Globals.session().setFlushMode(FlushMode.AUTO);
		}
	}

	public boolean hasConflicts() throws Exception
	{
		dublicate = hasArticle();
		Globals.session().setFlushMode(FlushMode.AUTO);
		return (dublicate != null);
	}

	public static String shortTitle(String name)
	{
		String result = name.toLowerCase().replaceAll("\\W", "");
		return result;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((id == null) ? 0 : id.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Article other = (Article) obj;
		if (id == null) {
			if (other.id != null)
				return false;
		} else if (!id.equals(other.id))
			return false;
		return true;
	}

}
