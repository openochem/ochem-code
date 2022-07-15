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
import java.net.URL;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.w3c.dom.Document;
import org.xml.sax.InputSource;

import qspr.entities.Article;
import qspr.entities.Author;

import com.eadmet.exceptions.UserFriendlyException;

public class ISBNUtility 
{  
	private static transient final Logger logger = LogManager.getLogger(ISBNUtility.class);

	/*public static void main(String[] args) 
	   throws Exception {

		  getArticleByISBN("0948404019", new Article());

	  }*/
	public static Article getArticleByISBN(String isbn, Article article) throws Exception{
		URL isbnurl;
		InputStream is;	  

		DocumentBuilderFactory domFactory = DocumentBuilderFactory.newInstance();
		domFactory.setNamespaceAware(true); // never forget this!
		DocumentBuilder builder = domFactory.newDocumentBuilder();

		isbnurl = new URL("http://isbndb.com/api/books.xml?access_key=D5ZWOTE5&index1=isbn&value1="+isbn);
		is = isbnurl.openStream();			
		InputSource inputSrc = new InputSource(is);			

		Document doc = builder.parse(inputSrc);
		XPathFactory factory = XPathFactory.newInstance();
		XPath xpath = factory.newXPath();
		if(xpath.evaluate("//BookData", doc).equalsIgnoreCase(""))
			throw new UserFriendlyException("Invalid ISBN no or book not found");

		String articleTitle = xpath.evaluate("//TitleLong",doc);
		if(articleTitle == null || articleTitle.equals(""))
			articleTitle = xpath.evaluate("//Title",doc);
		String authorString = xpath.evaluate("//AuthorsText",doc);
		String publisher = xpath.evaluate("//PublisherText",doc);
		String articleAbstract = xpath.evaluate("//Summary",doc);
		String isbn13 = xpath.evaluate("//BookData/@isbn13",doc);
		isbn = xpath.evaluate("//BookData/@isbn",doc);

		//get year
		String string = publisher.replaceAll("\\W", "");
		Pattern p = Pattern.compile("(19|20)\\d\\d");
		Matcher m = p.matcher(string);
		if(m.find())
		{
			String year = m.group();
			logger.info("publication year - "+year);
			article.setPublicationDate(year);
		}
		//compare isbn
		//add author
		authorString = authorString.replaceAll("\\([Ee]ditor\\)|\\[by\\]|edited by|translated by|et al.|\\.\\.\\.", "");
		authorString = authorString.replaceAll("\\[\\]|\\(\\)", "");

		String[] authorList;
		authorList = authorString.split(",|;|and");

		if(article != null)
			article.authors.clear();

		for (String author : authorList) {
			logger.info(author);
			Author auth = getAuthor(author.trim());
			article.addAuthor(auth);				
		}
		article.publisher = publisher;
		article.mediaType = "book";
		article.isbn = isbn.trim();
		article.isbn13 = isbn13.trim();
		article.time = new Timestamp(Calendar.getInstance().getTimeInMillis());
		article.setTitle(articleTitle);
		article.articleAbstract = articleAbstract;

		if(article.publicationDate == null)
			article.setCurrentPublicationDate();
		return article;
	}
	private static Author getAuthor(String author)
	{
		final List<String> nameaffixes = new ArrayList<String>(Arrays.asList(new String[]{"van", "von", "de", "le", "di"}));
		String lastname = "", firstname = "", initials = "";
		String[] tokens = author.split(" ");
		if(tokens.length == 1) {
			lastname = tokens[0];
		} else if(tokens.length == 2) {
			firstname = tokens[0];
			initials = tokens[0].substring(0,1) + ".";
			lastname = tokens[1];
		} else {
			firstname = tokens[0];
			initials = tokens[0].substring(0,1) + ".";
			lastname = tokens[tokens.length-1];
			for(int i = tokens.length-2; i > 0 && nameaffixes.contains(tokens[i]) ; --i) lastname = tokens[i] + " " + lastname;
		}
		return Author.get(firstname, lastname, initials);
	}

}
