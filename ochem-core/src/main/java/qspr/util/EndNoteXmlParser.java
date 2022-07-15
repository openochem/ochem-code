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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.List;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathExpression;
import javax.xml.xpath.XPathFactory;

import org.w3c.dom.Document;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;

import qspr.Globals;
import qspr.entities.Article;
import qspr.entities.Author;
import qspr.entities.Journal;

// TODO: Consider cases of books


public class EndNoteXmlParser
{

	Journal jrn;

	Article art;

	public Article parse(File f) throws Exception
	{
		FileInputStream fis = new FileInputStream(f);
		BufferedReader in = new BufferedReader(new InputStreamReader(fis));

		DocumentBuilderFactory domFactory = DocumentBuilderFactory.newInstance();
		domFactory.setNamespaceAware(true); // never forget this!
		DocumentBuilder builder = domFactory.newDocumentBuilder();
		InputSource inputSrc = new InputSource(fis);	

		Document doc = builder.parse(inputSrc);
		XPathFactory factory = XPathFactory.newInstance();
		XPath xpath = factory.newXPath();

		String title = xpath.evaluate("//xml/records/record/titles/title",doc);

		art = new Article();
		art.setTitle(title);

		String jtitle = xpath.evaluate("//xml/records/record/titles/secondary-title",doc);		

		jrn = Journal.getByTitle(jtitle);

		if (jrn == null)
		{
			jrn = new Journal();
			xpath.evaluate("//xml/records/record/ref-type/@name", doc);
			jrn.setTitle(jtitle);
			if (Globals.userSession().user != null)
			{
				jrn.owner = Globals.userSession().user;
				if(jrn.introducer == null)
					jrn.introducer = Globals.userSession().user;
			}
		}

		XPathExpression expr = xpath.compile("//xml/records/record/contributors/authors/author");
		Object result = expr.evaluate(doc, XPathConstants.NODESET);
		NodeList nodes = (NodeList) result;

		String s;

		for (int i = 0; i < nodes.getLength(); i++) 
		{
			s = xpath.evaluate(".", nodes.item(i));
			List<Author> auth = Author.getFromEN(s);
			for (Author author : auth) {
				art.addAuthor(author);
			}				
		}

		String pages = xpath.evaluate("//xml/records/record/pages", doc);
		art.pageNumbers = pages;

		String vol = xpath.evaluate("//xml/records/record/volume", doc);
		art.volume = vol;

		String iss = xpath.evaluate("//xml/records/record/number", doc);
		art.issue = iss;

		String date = xpath.evaluate("//xml/records/record/dates/year", doc);
		art.setPublicationDate(date);

		String abstr = xpath.evaluate("//xml/records/record/abstract", doc);
		art.articleAbstract = abstr;

		art.journal = jrn;
		Globals.session().save(jrn);
		in.close();

		return art; // return a newly created article
	}
}
