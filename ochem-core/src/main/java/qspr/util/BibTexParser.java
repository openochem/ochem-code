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
//
//import java.io.File;
//import java.io.FileReader;
//
//import qspr.Globals;
//import qspr.entities.Article;
//import qspr.entities.Author;
//import qspr.entities.Journal;
//import bibtex.dom.BibtexEntry;
//import bibtex.dom.BibtexFile;
//import bibtex.expansions.CrossReferenceExpander;
//import bibtex.expansions.ExpansionException;
//import bibtex.expansions.MacroReferenceExpander;
//import bibtex.expansions.PersonListExpander;
//import bibtex.parser.BibtexMultipleFieldValuesPolicy;
//import bibtex.parser.BibtexParser;
//
//public class BibTexParser
//{
//
//	Journal jrn;
//	
//	Article art;
//	
//	public String removeBraces(String s)
//	{
//		if(s.startsWith("{{") && s.endsWith("}}"))
//			return s.substring(2,s.length()-2);
//		if (s.startsWith("{") && s.endsWith("}"))
//			return s.substring(1, s.length()-1);
//		else return s;
//	}
//	
//	public Article parse(File f) throws Exception
//	{
//		BibtexFile bibtexFile = new BibtexFile();
//		BibtexParser parser = new BibtexParser(false);
//		parser.setMultipleFieldValuesPolicy(BibtexMultipleFieldValuesPolicy.KEEP_ALL);
//		String issn = "";
//		try {
//			parser.parse(bibtexFile, new FileReader(f));
//		} catch (Exception e) {
//			System.err.println("Fatal exception: ");
//			e.printStackTrace();
//			return null;
//		} 
//		
//		try {
//				MacroReferenceExpander mrexpander = new MacroReferenceExpander(true, true, false, false);
//				mrexpander.expand(bibtexFile);			
//				CrossReferenceExpander crexpander = new CrossReferenceExpander(false);
//				crexpander.expand(bibtexFile);
//				PersonListExpander plexpander = new PersonListExpander(true, true, true);
//				plexpander.expand(bibtexFile);
//		} catch (ExpansionException e1) {
//			e1.printStackTrace();
//			return null;
//		}
//		
//		for (Object elem : bibtexFile.getEntries()) {
//			if (elem instanceof BibtexEntry) 
//			{
//				BibtexEntry entry = (BibtexEntry)elem;
//				art = new Article();
//				//Article title
//				art.setTitle(removeBraces(entry.getFieldValue("title").toString()));
//				//get the author 
//				art.authors = Author.getFromBT(removeBraces(entry.getFieldValue("author").toString()));
//				if(entry.getFieldValue("abstract") != null)
//					art.articleAbstract = removeBraces(entry.getFieldValue("abstract").toString());
//				if(entry.getFieldValue("volume") != null)
//					art.volume = removeBraces(entry.getFieldValue("volume").toString());
//				if(entry.getFieldValue("number") != null)
//					art.issue = removeBraces(entry.getFieldValue("number").toString());
//				if(entry.getFieldValue("pages") != null)
//					art.pageNumbers = removeBraces(entry.getFieldValue("pages").toString());
//				if(entry.getFieldValue("doi") != null)
//				{
//					String doi = removeBraces(entry.getFieldValue("doi").toString());
//					art.doi = doi.replaceAll("DOI:", "").trim();
//				}
//				if(entry.getFieldValue("issn") != null)
//				{
//					jrn = Journal.getByISSN(removeBraces(entry.getFieldValue("issn").toString()).replaceAll("-", "").trim());
//				}
//				else{
//					jrn = Journal.getByTitle(removeBraces(entry.getFieldValue("journal").toString()).trim());
//				}
//					
//				
//				if (jrn.id == null)
//				{
//					if(entry.getFieldValue("journal") != null)
//						jrn.title = removeBraces(entry.getFieldValue("journal").toString()).trim();
//					if(entry.getFieldValue("publisher") != null)
//						jrn.publisher = removeBraces(entry.getFieldValue("publisher").toString()).trim();
//					if(entry.getFieldValue("journal-iso") != null)
//						jrn.abbreviation = removeBraces(entry.getFieldValue("journal-iso").toString()).trim();
//					else
//						jrn.abbreviation = jrn.title;
//					jrn.issn = removeBraces(entry.getFieldValue("issn").toString()).replaceAll("-", "").trim();
//					if (Globals.userSession().user != null)
//					{
//						jrn.owner = Globals.userSession().user;
//						if(jrn.introducer == null)
//							jrn.introducer = Globals.userSession().user;
//					}
//					Globals.session().save(jrn);
//				}
//				art.journal = jrn;
//				
//				if(entry.getFieldValue("month") != null && entry.getFieldValue("year") != null)
//					art.setDate(removeBraces(entry.getFieldValue("month").toString())+" "+removeBraces(entry.getFieldValue("year").toString()));
//				else if(entry.getFieldValue("month") == null && entry.getFieldValue("year") != null)
//					art.setDate(removeBraces(entry.getFieldValue("year").toString()));
//				else
//					art.setDate("Jan 1900");
//				
//				return art;
//			}
//		}
//		return null;
//	}
//}
