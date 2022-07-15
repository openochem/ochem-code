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
import java.net.URLEncoder;
import java.text.SimpleDateFormat;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathFactory;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.junit.Test;
import org.w3c.dom.Document;

import qspr.entities.Journal;

import com.eadmet.exceptions.UserFriendlyException;

public class JournalUtility
{
	private static transient final Logger logger = LogManager.getLogger(JournalUtility.class);

	public static void main(String[] args) throws Exception
	{
		//		
		//		Globals.startMainTransaction();
		//1060-1325 1474-0346 0198-6325
		
		Configurator.setLevel(logger.getName(), Level.DEBUG);
		//logger.setLevel(Level.DEBUG); # N.B.!
		fetchJournal("0198-6325");

		//		  Globals.session().saveOrUpdate(journal);
		//		  
	}

	@Test(timeout = 180000)
	public static Journal fetchJournal(String searchTerm) throws Exception
	{
		// TODO Rob 11.07.13: how to handle more than one ID
		String journal_id = searchJournal(searchTerm); 
		Document doc = fetchJournalByJournalId(journal_id);
		Journal journal = fillFromDocument(doc);
		return journal;
	}

	/** example esummary result
	 * 	<eSummaryResult>
	 * 		<DocSum>
	 * 			<Id>5606</Id>
	 * 			<Item Name="Title" Type="String">Medicinal research reviews</Item>
	 * 			<Item Name="MedAbbr" Type="String">Med Res Rev</Item>
	 * 			<Item Name="IsoAbbr" Type="String">Med Res Rev</Item>
	 * 			<Item Name="NlmId" Type="String">8103150</Item>
	 * 			<Item Name="pISSN" Type="String">0198-6325</Item>
	 * 			<Item Name="eISSN" Type="String">1098-1128</Item>
	 * 			<Item Name="PublicationStartYear" Type="String">1981</Item>
	 * 			<Item Name="PublicationEndYear" Type="String"/>
	 * 			<Item Name="Publisher" Type="String">Wiley</Item>
	 * 			<Item Name="Language" Type="String">eng</Item>
	 * 			<Item Name="Country" Type="String">United States</Item>
	 * 			<Item Name="BroadHeading" Type="List">
	 * 				<Item Name="string" Type="String">Pharmacology</Item>
	 * 			</Item>
	 * 			<Item Name="ContinuationNotes" Type="String"/>
	 * 		</DocSum>
	 * 	</eSummaryResult>
	 * @param journal2 
	 */

	@Test(timeout = 180000)
	private static String searchJournal(String searchTerm) throws Exception
	{
		DocumentBuilderFactory domFactory = DocumentBuilderFactory.newInstance();
		domFactory.setNamespaceAware(true); // never forget this!
		DocumentBuilder builder = domFactory.newDocumentBuilder();

		URLS journalSearch = new URLS("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nlmcatalog&term=\"" + URLEncoder.encode(searchTerm, "UTF-8")+"\"");
		logger.info(journalSearch);

		InputStream is = journalSearch.openStream();	
		Document doc = builder.parse(is);
		is.close();

		XPathFactory factory = XPathFactory.newInstance();
		XPath xpath = factory.newXPath();		

		return xpath.evaluate("//IdList/Id", doc);
	}

	@Test(timeout = 180000)
	private static Document fetchJournalByJournalId(String journalID) throws Exception
	{
		DocumentBuilderFactory domFactory = DocumentBuilderFactory.newInstance();
		domFactory.setNamespaceAware(true); // never forget this!
		DocumentBuilder builder = domFactory.newDocumentBuilder();

		if (journalID.equalsIgnoreCase(""))
			throw new UserFriendlyException("Invalid journal ID. No journal found for: \"" + journalID + "\"");

		//get the details for journal from pubmed
		URLS journalURL = new URLS("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nlmcatalog&id=" + journalID + "&retmode=xml");
		logger.warn(journalURL);
		InputStream is = journalURL.openStream();
		Document doc = builder.parse(is);
		//id.close();

		return doc;
	}

	@Test(timeout = 180000)
	private static Journal fillFromDocument(Document doc) throws Exception
	{
		XPathFactory factory = XPathFactory.newInstance();
		XPath xpath = factory.newXPath();

		Journal	journal = new Journal();
		journal.setTitle(xpath.evaluate("//TitleMain/Title", doc));

		String abbreviation = xpath.evaluate("//MedlineTA", doc).trim();
		journal.setAbbreviation(abbreviation);

		journal.setISSN(xpath.evaluate("//ISSNLinking", doc));

		journal.setLink(xpath.evaluate("//ELocationID", doc));

		journal.publisher = xpath.evaluate("//ImprintFull", doc).trim();

		String publicationStartYear = xpath.evaluate("//PublicationFirstYear", doc);
		SimpleDateFormat sdf = new SimpleDateFormat("yyyy");
		journal.publish_date = sdf.parse(publicationStartYear);

		logger.info(journal.getTitle() + " - " + journal.getAbbreviation() + " - " + journal.getISSN() + " - " + journal.publish_date + " - " + journal.publisher);

		return journal;
	}


}
