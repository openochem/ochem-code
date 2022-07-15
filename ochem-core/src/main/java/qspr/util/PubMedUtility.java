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

import java.io.InputStream;
import java.sql.Timestamp;
import java.text.SimpleDateFormat;
import java.util.Calendar;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathExpression;
import javax.xml.xpath.XPathFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.w3c.dom.Document;
import org.w3c.dom.NodeList;

import qspr.Globals;
import qspr.entities.Article;
import qspr.entities.Author;
import qspr.entities.Journal;

import com.eadmet.exceptions.UserFriendlyException;

public class PubMedUtility 
{  
	private static transient final Logger logger = LogManager.getLogger(PubMedUtility.class);

	/*	public static void main(String[] args) 
	   throws Exception {
		  Article art = getArticleByPubmedId(456);
	  }*/

	public static Article getArticleByPubmedId(long pmId) throws Exception
	{
		Document doc = fetchById(pmId);
		return fillFromDocument(doc);
	}

	private static Document fetchById(long pmId) throws Exception
	{
		URLS pubmedurl = new URLS(String.format("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=%d&retmode=xml",pmId));
		logger.info(pubmedurl.toString());
		DocumentBuilderFactory domFactory = DocumentBuilderFactory.newInstance();
		domFactory.setNamespaceAware(true);
		DocumentBuilder builder = domFactory.newDocumentBuilder();
		InputStream is = pubmedurl.openStream();
		return builder.parse(is);
	}

	private static Article fillFromDocument(Document doc) throws Exception
	{
		XPathFactory factory = XPathFactory.newInstance();
		XPath xpath = factory.newXPath();

		String pmId = xpath.evaluate("//PMID", doc);
		if (pmId.equalsIgnoreCase(""))
			throw new UserFriendlyException("Invalid PubMedId or Pubmed id not found");

		String abbreviation = xpath.evaluate("//Journal/ISOAbbreviation", doc);
		String jouTitle = xpath.evaluate("//Journal/ArticleTitle", doc);
		String volume = xpath.evaluate("//JournalIssue/Volume", doc);
		String issue = xpath.evaluate("//JournalIssue/Issue", doc);
		String year = xpath.evaluate("//PubDate/Year", doc);
		String month = xpath.evaluate("//PubDate/Month", doc);
		String articleAbstract = xpath.evaluate("//AbstractText", doc);
		String articleTitle = xpath.evaluate("//ArticleTitle", doc);
		String affiliation = xpath.evaluate("//Affiliation", doc);
		String articlepn = xpath.evaluate("//MedlinePgn", doc);
		String issn = xpath.evaluate("//ISSN", doc);
		String doi = xpath.evaluate("//ELocationID", doc);

		if (articleAbstract.equals(""))
			articleAbstract = xpath.evaluate("//Abstract", doc);

		if (jouTitle.equals(""))
			jouTitle = xpath.evaluate("//Journal/Title", doc);

		XPathExpression expr = xpath.compile("//Author");
		Object result = expr.evaluate(doc, XPathConstants.NODESET);
		NodeList nodes = (NodeList) result;

		if (month.equals("") && year.equals(""))
		{
			String meddate = xpath.evaluate("//PubDate/MedlineDate", doc);
			if (!meddate.equals(""))
			{
				String[] date = meddate.split(" ");
				year = date[0];
				if(date.length > 1) {
					if (date[1].contains("-"))
						month = date[1].split("-")[0];
					else
						month = date[1];
				}
				else
					month = "Jan";
			}
		}

		if (month.equals(""))
			month = "Jan";
		if (year.equals(""))
			year = "1900";

		Journal journal = Journal.getByISSN(issn);

		if (journal == null)
			journal = Journal.getByTitle(jouTitle);

		// There is no such journal in our database, create one
		if (journal == null)
		{
			journal = new Journal();
			journal.setISSN(issn);
			journal.setTitle(jouTitle);
			journal.setAbbreviation(abbreviation);
			if (Globals.userSession().user != null)
			{
				journal.owner = Globals.userSession().user;
				if (journal.introducer == null)
					journal.introducer = Globals.userSession().user;
			}
			Globals.session().saveOrUpdate(journal);
		}

		Article article = new Article();
		article.mediaType = "article";
		article.pmid = Long.valueOf(pmId);
		article.doi = doi;
		article.time = new Timestamp(Calendar.getInstance().getTimeInMillis());
		article.journal = journal;
		if (!volume.equals(""))
			article.volume = volume;
		article.issue = issue;
		article.setTitle(articleTitle);
		article.pageNumbers = articlepn;try {
			article.publicationDate = (new SimpleDateFormat("MMM yyyy")).parse(month + " " + year);
		}catch(java.text.ParseException e) {
			article.publicationDate = (new SimpleDateFormat("mm yyyy")).parse(month + " " + year);
		}
		article.articleAbstract = articleAbstract;
		article.affiliation = affiliation;
		if (article.authors != null)
			article.authors.clear();


		for (int i = 0; i < nodes.getLength(); i++)
		{
			String lastName = xpath.evaluate("LastName", nodes.item(i));
			String firstName = xpath.evaluate("FirstName", nodes.item(i));
			String initials = xpath.evaluate("Initials", nodes.item(i));
			if (initials.equals(""))
			{
				initials = xpath.evaluate("ForeName", nodes.item(i));
				if (initials.equals(""))
				{
					if (!firstName.equals(""))
					{
						String[] iniarr = firstName.split(" ");
						for (int j = 0; j < iniarr.length; j++)
						{
							initials += iniarr[j].substring(0, 1) + " ";
						}
					}
				}
			}

			if (firstName.equals(""))
				firstName = xpath.evaluate("ForeName", nodes.item(i));
			logger.info("Name " + lastName + " " + firstName + " " + initials.trim());
			article.addAuthor(Author.get(firstName, lastName, initials.trim()));
		}
		return article;
	}

}
