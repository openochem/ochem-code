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

package qspr.controllers;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.sql.Timestamp;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Date;
import java.util.List;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.Hibernate;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.MarshallingOption;
import qspr.ThreadScope;
import qspr.business.WebFilters;
import qspr.entities.Alert;
import qspr.entities.Article;
import qspr.entities.ArticleUserPdf;
import qspr.entities.Author;
import qspr.entities.Basket;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Journal;
import qspr.entities.Property;
import qspr.exception.DublicateArticleException;
import qspr.frontend.BrowserModel;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;
import qspr.util.EndNoteParser;
import qspr.util.EndNoteXmlParser;
import qspr.util.ISIParser;
import qspr.util.JournalUtility;
import qspr.util.PubMedUtility;
import qspr.util.RISParser;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.exceptions.UserFriendlyException;

@Controller
public class ArticleController extends BrowserWrapper 
{
	private static transient final Logger logger = LogManager.getLogger(ArticleController.class);

	Article article;

	public ArticleController()
	{
		sessionRequired = true;
	}

	public ModelAndView show(HttpServletRequest request, HttpServletResponse response) 
			throws Exception 
	{
		return new WebModel().setList(Globals.getTaginationFilters(Property.class)).setTemplate("article-browser").getModelAndView();
	}

	@SuppressWarnings("unchecked")
	public ModelAndView edit(HttpServletRequest request, HttpServletResponse response) 
			throws Exception
	{
		if (Globals.userSession().user == null)
			return redirect("user/pleaseregister.do");

		List<Author> authors = ((List<Author>) request.getSession().getAttribute("AuthorList"));
		if(authors != null)
			request.getSession().removeAttribute("AuthorList");

		Article art;
		if (assertParam("fromsession"))
			art = (Article) ThreadScope.get().localRequest.getSession().getAttribute("article");
		else
			art = (Article) Globals.session().get(Article.class, getLongParam("id"));
		if (art == null)
		{
			art = new Article();
			art.mediaType = request.getParameter("media-type");
			if(art.mediaType.equals("all") || art.mediaType.equals("temporal")){
				if(art.mediaType.equals("temporal")){
					art.journal = (Journal) Globals.session().get(Journal.class, QSPRConstants.UNPUBLISHED_JOURNAL);
					art.setCurrentPublicationDate();
				}
				art.mediaType="article";
			}
			art.id = Long.valueOf(-1);
			art.introducer = Globals.userSession().user;
			art.owner = Globals.userSession().user;

		}
		return new WebModel(art).setTemplate("article-edit").getModelAndView();
	}

	/**
	 * The article profile page
	 */
	public ModelAndView profile(HttpServletRequest request, HttpServletResponse response)
	{
		Globals.setMarshallingOption(MarshallingOption.ARTICLE_PENDING_TASKS);
		Article art = (Article) Globals.session().get(Article.class, getLongParam("id"));
		return new WebModel(art).setTemplate("article-profile").getModelAndView();
	}

	@SuppressWarnings("unchecked")
	public ModelAndView saveauthors(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		Article art = (Article) Globals.session().get(Article.class, getLongParam("id"));
		art.doCheckRights();

		if (assertParam("name"))
		{
			String[] authorName = req.getParameterValues("name");

			art.authors.clear();

			for (int i = 0; i < authorName.length; i++)
			{
				if (authorName[i].contains(","))
				{
					String[] splittedNames = authorName[i].split(";");
					Author author;
					if (splittedNames.length > 1)
						author = Author.get("", splittedNames[0].trim(), splittedNames[1].trim());
					else
						author = Author.get("", splittedNames[0].trim(), "");

					art.authors.add(author);
				}
				else
				{
					String[] splittedNames = authorName[i].split(";");
					for (String name : splittedNames)
					{
						Author author = Author.get(name.trim());
						if (!art.authors.contains(author))
							art.authors.add(author);
					}
				}
			}

			//save author
			art.owner = Globals.userSession().user;
			Globals.session().save(art);
		}
		else if ("swap".equals(req.getParameter("action")))
		{
			List<Author> authors = new ArrayList<Author>();
			List<Author> tempAuthors = (List<Author>) req.getSession().getAttribute("AuthorList") == null ? art.authors : (List<Author>) req.getSession()
					.getAttribute("AuthorList");
			for (Author author : tempAuthors)
			{
				String lastName = author.firstName;
				String firstName = author.lastName;
				author.firstName = firstName.replace(",", "");
				author.lastName = lastName.replace(",", "");
				if (author.initials != null && !author.initials.equals("") && !firstName.equals(""))
					author.initials = firstName.substring(0, 1) + ".";
				authors.add(author);
				Globals.session().evict(author);
			}
			req.getSession().setAttribute("AuthorList", authors);
			return new WebModel(new WebList().loadFromList(authors)).getModelAndView();
		}
		//		}
		return new WebModel(new WebList().loadFromList(art.authors)).getModelAndView();
	}

	public ModelAndView listauthors(HttpServletRequest req, HttpServletResponse res)
	{	
		Article art;
		if (assertParam("fromsession"))
			art = (Article) ThreadScope.get().localRequest.getSession().getAttribute("article");
		else
			art = (Article) Globals.session().get(Article.class, getLongParam("id"));
		if (art != null)
		{
			List<Author> authors = art.authors;
			return new WebModel(new WebList().loadFromList(authors)).getModelAndView();
		}
		return new WebModel().getModelAndView();
	}



	public ModelAndView managePdfs(HttpServletRequest req, HttpServletResponse res) throws Exception {
		Article a = Article.getById(Long.valueOf(req.getParameter("id")));
		String action = req.getParameter("action");
		if (!Globals.isValidatedUser())
			throw new UserFriendlyException("Non-validated users cannot attach files");
		if ("attach".equals(action))
			uploadPdf(a);
		else if ("remove".equals(action)) {
			a.removeAttachment(ArticleUserPdf.PDF);
		}
		return redirect("article/profile.do?id=" + a.id);
	}

	private ModelAndView actionDelete(Article art) throws Exception{
		// Delete
		if (art != null)
		{
			if (assertParam("deletion-confirmation"))
			{
				// Delete my records from trash
				// Dirty.. But no other solution except direct SQL query, since HQL does not delete cascade in HQL qith many to many / Midnighter
				Globals.session().createSQLQuery("delete ExperimentalPropertyName from ExperimentalPropertyName left join ExperimentalProperty using (exp_property_id) " +
						"where article_id=" + art.id + " and introducer_id=" + Globals.userSession().user.id
						).executeUpdate();

				Globals.session().flush();

				Globals.session().createQuery("delete from ExperimentalProperty where deleted is not null and introducer=:user and article=:article")
				.setParameter("user", Globals.userSession().user)
				.setParameter("article", art)
				.executeUpdate();

				// Deattach my files from article
				art.pdfs.clear();

				Globals.session().flush();
				Globals.restartAllTransactions(true);
			}

			art.doCheckRights();
			long myTrushCount = (Long) Globals.session().createCriteria(ExperimentalProperty.class).add(Restrictions.isNotNull("deleted"))
					.add(Restrictions.eq("article", art)).add(Restrictions.eq("introducer", Globals.userSession().user))
					.setProjection(Projections.count("id")).list().get(0);
			int myPdfCount = art.pdfs.size();
			if (myTrushCount != 0 || myPdfCount != 0)
			{
				String message = "Following entities are connected to the article being deleted:\n";
				if (myTrushCount != 0)
					message += "" + myTrushCount + " records in your trash\n";
				if (myPdfCount != 0)
					message += "" + myPdfCount + " file(s) attached to this article\n";

				return new WebModel(new Alert(message)).getModelAndView();
			}

			int pdfCount = Globals.session().createCriteria(ArticleUserPdf.class).add(Restrictions.eq("article", art))
					.add(Restrictions.ne("user", Globals.userSession().user)).list().size();
			long trushCount = (Long) Globals.session().createCriteria(ExperimentalProperty.class).add(Restrictions.isNotNull("deleted"))
					.add(Restrictions.eq("article", art)).setProjection(Projections.count("id")).list().get(0);
			if (pdfCount != 0 || trushCount != 0)
				throw new UserFriendlyException("Article cannot be deleted since there are hidden records or attached files added by other users");

			Globals.session().delete(art);
		}
		return new WebModel().getModelAndView();
	}

	private ModelAndView actionReloadPM(Article art) throws Exception{
		Long pmID = getLongParam("n-pubmed");
		Article pubMedArticle = PubMedUtility.getArticleByPubmedId(pmID);
		if (pubMedArticle == null)
			throw new UserFriendlyException("No article found in PubMed for given PMID: \"" + pmID + "\"");

		if (art == null)
			return new WebModel(pubMedArticle).getModelAndView();
		else
		{
			art = pubMedArticle.hasArticle(); 

			if (art != null)
			{
				ThreadScope.get().localRequest.getSession().setAttribute("article", art);
				Hibernate.initialize(art.pdfs);
				Hibernate.initialize(art.authors);
				Globals.session().evict(art);

				pubMedArticle.setLink(art.getUrl());
				pubMedArticle.setComment(art.comment);
			}
			return new WebModel(pubMedArticle).getModelAndView();
		}		
	}


	private ModelAndView actionLoadISBN(String isbn13 ) throws Exception{
		ThreadScope.get().localRequest.getSession().removeAttribute("article");
		Article art = Article.getByISBN(isbn13.trim());
		if (art.id != null)
			throw new DublicateArticleException(art);
		if (art.publicationDate == null)
			art.setCurrentPublicationDate();
		ThreadScope.get().localRequest.getSession().setAttribute("article", art);
		Hibernate.initialize(art.pdfs);
		Hibernate.initialize(art.authors);
		Globals.session().evict(art);
		return new WebModel(art).getModelAndView();
	}

	private Article actionEditArticle() throws Exception{

		Article art;
		HttpServletRequest mp = ThreadScope.get().localRequest;


		String[] ids = mp.getParameterValues("id");
		if ( ! ids[0].equals("") && ! ids[0].equals("-1"))
		{
			Long id = Long.valueOf(ids[0]);
			art = (Article)Globals.session().get(Article.class, id);
		}
		else
		{
			art = new Article();
			art.mediaType = Globals.htmlEntityDecode(mp.getParameter("media-type"));
		}

		if (assertParam("n-pubmed"))
			art.pmid = Long.valueOf(Globals.htmlEntityDecode(mp.getParameter("n-pubmed")));
		
		if (assertParam("n-journal-key"))
		{	
			Long journalID = Long.valueOf(Globals.htmlEntityDecode(mp.getParameter("n-journal-key")));
			art.journal = (Journal) Globals.session().get(Journal.class, journalID);
		} 
		else if (assertParam("n-journal-title"))
		{
			String title = Globals.htmlEntityDecode(mp.getParameter("n-journal-title"));
			Journal dbJournal = Journal.getByTitle(title);
			if (dbJournal == null) // no journal found in ochem
			{
				dbJournal = JournalUtility.fetchJournal(title);
				if(dbJournal ==null )throw new UserFriendlyException("This journal " + title + " is abdent in OCHEM. Please, create it before uploading the article");
			}
			art.journal = dbJournal;
		}

		if (assertParam("n-title"))
			art.setTitle(Globals.htmlEntityDecode(mp.getParameter("n-title")));
		if (assertParam("n-abstract"))
			art.articleAbstract = Globals.htmlEntityDecode(mp.getParameter("n-abstract"));

		// only book stuff
		if (assertParam("n-isbn"))
			art.isbn = Globals.htmlEntityDecode(mp.getParameter("n-isbn"));
		if (assertParam("n-isbn13"))
			art.isbn13 = Globals.htmlEntityDecode(mp.getParameter("n-isbn13"));
		if (assertParam("n-publisher"))
			art.publisher = Globals.htmlEntityDecode(mp.getParameter("n-publisher"));
		// only book stuff end

		Date date = art.publicationDate;
		if (assertParam("n-date"))
			art.setPublicationDate(Globals.htmlEntityDecode(mp.getParameter("n-date")));				

		if (art.id != null && art.publicationDate != null && !art.publicationDate.equals(date))
		{
			// Invalidate the "primary record" column for the related records
			Globals.session().createQuery("update ExperimentalProperty set firstEntry=null where article=:article")
			.setParameter("article", art)
			.executeUpdate();
		}

		art.isChapter = assertParam("is-chapter");

		if (assertParam("n-volume"))
			art.volume = Globals.htmlEntityDecode(mp.getParameter("n-volume"));
		if (assertParam("n-issue"))
			art.issue = Globals.htmlEntityDecode(mp.getParameter("n-issue"));
		if (assertParam("n-pages"))
			art.pageNumbers = Globals.htmlEntityDecode(mp.getParameter("n-pages"));

		if (assertParam("n-affiliation"))
			art.affiliation = Globals.htmlEntityDecode(mp.getParameter("n-affiliation"));

		if (assertParam("n-doi"))
			art.doi = Globals.htmlEntityDecode(mp.getParameter("n-doi"));

		// save whatever is in the form
		art.setLink(Globals.htmlEntityDecode(mp.getParameter("n-url")));
		art.setComment(Globals.htmlEntityDecode(mp.getParameter("n-comment")));


		if (assertParam("parent-book-id"))
		{
			art.parent = (Article) Globals.session().get(Article.class, getLongParam("parent-book-id"));
			assert art.parent.mediaType.equals("book");
		}

		if (Globals.userSession().user != null)
		{
			if (assertParam("n-delpdf"))
				art.removeAttachment(ArticleUserPdf.PDF);
			if (assertParam("n-batch"))
				art.removeAttachment(ArticleUserPdf.EXCEL);

			uploadPdf(art);
		}

		return art;
	}

	private void uploadPdf(Article art) throws IOException{
		FileInputStream ff = null;
		try
		{
			File f = Globals.getUploadedFile();
			byte[] pdf = new byte[(int)f.length()];
			ff=new FileInputStream(f);
			ff.read(pdf);
			art.attachFile(pdf, getFileType(f.getName()));
		} catch (Exception e)
		{
			logger.info("No file uploaded");
		}finally{
			if(ff!=null)ff.close();
		}
	}

	private int getFileType(String fileName) {
		if (fileName.toLowerCase().endsWith(QSPRConstants.PDF))
			return ArticleUserPdf.PDF;
		throw new UserFriendlyException("Unsupported attachment type for file " + fileName+ " Only  PDF file can be stored");
	}


	public ModelAndView action(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		// TODO: Needs some refactoring. Forbid edit on server-side when PubMed present
		Article art = null;

		if (assertParam("id"))
			art = (Article) Globals.session().get(Article.class, getLongParam("id"));

		String action = request.getParameter("action");

		if (action.equals("delete"))
			return actionDelete(art);

		if (action.equals("reloadpm"))
			return actionReloadPM(art);

		if (action.equals("load_isbn"))
			return actionLoadISBN(request.getParameter("load-isbn"));

		art=actionEditArticle();
		
		if(art.pmid != null && art.authors.size() == 0 && art.publicationDate == null) // new entry!
			art = PubMedUtility.getArticleByPubmedId(art.pmid);

		// Ownership and rights	
		art.doCheckRights();
		if (Globals.userSession().user != null)
		{
			if (art.id == null || art.id < 0)
				art.introducer = Globals.userSession().user;
			if(art.introducer == null)
				art.introducer = Globals.userSession().user;

			art.owner = Globals.userSession().user;
		}

		art.time = new Timestamp(Calendar.getInstance().getTimeInMillis());
		if (art.publicationDate == null)
			art.setCurrentPublicationDate();
		if (art.hasConflicts())
			throw new DublicateArticleException(art.dublicate);

		Globals.session().saveOrUpdate(art);
		return new WebModel(art).getModelAndView();	
	}

	public ModelAndView upload(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		HttpServletRequest mp = ThreadScope.get().localRequest;
		Article art;
		File f = Globals.getUploadedFile();

		if (mp.getParameter("format").equals("ris"))
		{
			RISParser rp = new RISParser();
			art = rp.parse(f);
		} 
		else if (mp.getParameter("format").equals("en"))
		{
			EndNoteParser rp = new EndNoteParser();
			art = rp.parse(f);
		} 
		else if (mp.getParameter("format").equals("ex"))
		{
			EndNoteXmlParser rp = new EndNoteXmlParser();
			art = rp.parse(f);					
		} 
		else
		{
			ISIParser isip = new ISIParser();
			art = isip.parse(f);					
		}

		art.mediaType = "article";
		// Ownership and rights	
		art.doCheckRights();
		if (Globals.userSession().user != null)
		{
			if (art.id == null || art.id < 0)
				art.introducer = Globals.userSession().user;
			if(art.introducer == null)
				art.introducer = Globals.userSession().user;

			art.owner = Globals.userSession().user;
		}

		art.time = new Timestamp(Calendar.getInstance().getTimeInMillis());

		if(art.publicationDate == null)
			art.setCurrentPublicationDate();

		if(art.getTitle() != null)
		{
			if(mp.getParameter("overwrite") != null){
				Article existArticle = art.hasArticle();
				if(existArticle != null)
				{
					existArticle.doCheckRights();
					art.id = existArticle.id;
					art.introducer = existArticle.introducer;
					art.pdfs = existArticle.pdfs;
					Globals.session().merge(art);
					return new WebModel(art).getModelAndView();	
				}

			}else{
				if(art.hasConflicts())
					throw new DublicateArticleException(art.dublicate);
			}
		}
		else
			throw new UserFriendlyException("Please choose correct file format before uploading");

		if (art.issue != null)
			art.issue = art.issue.substring(0, Math.min(art.issue.length(), 20));

		if (art.pageNumbers != null)
			art.pageNumbers = art.pageNumbers.substring(0, Math.min(art.pageNumbers.length(), 45));

		Globals.session().saveOrUpdate(art);

		return new WebModel(art).getModelAndView();	
	}

	public ModelAndView list(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		WebList list = new WebList();
		WebFilters filters = formFilters(request);

		Globals.session().enableFilter("no-experimental-data"); //Important!

		Criteria articleCriteria = Globals.session().createCriteria(Article.class, "art");

		if (filters.has("id") && ! "".equals(filters.get("id")))
			filterById(articleCriteria);
		else
			if (!filters.has("query") && assertParam("media-type")){ //In autocomplete we would want both - articles and books
				if(filters.get("media-type").equals("article"))
					articleCriteria.add(Restrictions.and(Restrictions.eq("mediaType", "article"),Restrictions.ne("journal.id",QSPRConstants.UNPUBLISHED_JOURNAL)));
				else
					if(filters.get("media-type").equals("book"))
						articleCriteria.add(Restrictions.eq("mediaType", "book"));
					else
						if(filters.get("media-type").equals("temporal"))
							articleCriteria.add(Restrictions.eq("journal.id",QSPRConstants.UNPUBLISHED_JOURNAL));
			}
		if(filters.has("journal_id"))
			articleCriteria.add(Restrictions.eq("journal.id", filters.getLong("journal_id")));
		if (filters.has("title"))
			articleCriteria.add(Restrictions.like("title", "%"+filters.get("title")+"%"));
		if (filters.has("query"))
			articleCriteria.add(Restrictions.like("title", "%"+filters.get("query")+"%"));
		if (filters.has("year"))
		{
			DateFormat format = new SimpleDateFormat("yyyy-MM-dd hh:mm:ss");
			Date startDate = (Date)format.parse(""+filters.getInteger("year")+"-01-01 00:00:00");
			Date endDate = (Date)format.parse(""+filters.getInteger("year")+"-12-31 23:59:59");
			articleCriteria.add(Restrictions.between("publicationDate", startDate, endDate));
		}
		if (filters.has("name"))
		{
			articleCriteria.createCriteria("authors").
			add(Restrictions.or(
					Restrictions.like("lastName", "%"+filters.get("name")+"%"),
					Restrictions.like("initials", "%"+filters.get("name")+"%"))
					);
		}
		if (filters.has("username"))
		{
			articleCriteria.createAlias("introducer","intr");
			articleCriteria.createAlias("owner","modif");
			articleCriteria.add(Restrictions.or(Restrictions.like("intr.login", "%"+filters.get("username")+"%"), Restrictions.like("modif.login", "%"+filters.get("username")+"%")));
		}

		if (Globals.applyTaginationFilters(articleCriteria, Property.class, "p"))
			articleCriteria.createAlias("properties", "p");

		if (filters.has(QSPRConstants.PDF)) {
			if (!Globals.isValidatedUser())
				throw new UserFriendlyException("You do not have sufficient privileges for this action");

			articleCriteria.createAlias("pdfs", "pdfs")
			.createAlias("pdfs.attachedFile", "attachedFile")
			.add(Restrictions.eq("attachedFile.id", getLongParam(QSPRConstants.PDF)));
		}
		if (assertParam("basket"))
		{
			Basket basket = (Basket) Globals.session().get(Basket.class, getLongParam("basket"));
			articleCriteria.createAlias("experimentalProperties", "ep").createAlias("ep.basketEntries", "be")
			.add(Restrictions.eq("be.basket", basket));
		}

		if (filters.has("article-identifier"))
		{
			String identifier=filters.get("article-identifier").toUpperCase();
			if (identifier.startsWith("Q") || identifier.startsWith("A"))
				articleCriteria.add(Restrictions.eq("id", Long.valueOf(identifier.substring(1))));
			else
				articleCriteria.add(Restrictions.eq("pmid", Long.valueOf(identifier)));
		}

		if ("article".equals(filters.get("media-type")))
		{

			if (filters.has("volume"))
				articleCriteria.add(Restrictions.eq("volume", filters.get("volume")));

			if (filters.has("journal"))
			{
				articleCriteria.createAlias("journal", "jrn");
				articleCriteria.add(Restrictions.like("jrn.title", "%"+filters.get("journal")+"%"));
			}
		}

		if ("book".equals(filters.get("media-type")))
		{
			if (filters.has("identifier"))
			{
				String identifier = filters.get("identifier").replaceAll("-", "").replaceAll(" ", "");
				if (identifier.startsWith("Q") || identifier.startsWith("q") || identifier.startsWith("A") || identifier.startsWith("a"))
					articleCriteria.add(Restrictions.eq("id", Long.valueOf(identifier.substring(1))));
				else
					articleCriteria.add(Restrictions.or(Restrictions.eq("isbn", identifier), Restrictions.eq("isbn13", identifier)));
			}
			if (filters.has("hide-chapters"))
				articleCriteria.add(Restrictions.eq("isChapter", false));
		}
		if(filters.has("models"))
			articleCriteria.add(Restrictions.isNotEmpty("models"));
		Integer pageNum;

		if (filters.has("pagenum"))
			pageNum = filters.getInteger("pagenum");
		else
			pageNum = 1;

		articleCriteria.addOrder(Order.desc("id"));

		list.column = "art.id";
		list.useEntity(Article.class).loadDistinctFromCriteria(articleCriteria, pageNum, getPageSize(5));

		return new BrowserModel().setFilters(filters).setObject(list).getModelAndView();
	}

	public ModelAndView refreshtags(HttpServletRequest request, HttpServletResponse response)
	{
		Globals.session().createSQLQuery("delete from ArticleProperty").executeUpdate();
		Globals.session().createSQLQuery("insert into ArticleProperty(article_id, property_id) (select article_id, property_id from  ExperimentalProperty group by article_id, property_id)").executeUpdate();
		return redirect("article/show.do");
	}

}
