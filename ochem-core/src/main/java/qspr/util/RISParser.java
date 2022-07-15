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
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;

import qspr.Globals;
import qspr.entities.Article;
import qspr.entities.Author;
import qspr.entities.Journal;

enum RISTag
{
	UNKNOWN, CITATION_TYPE, CITATION_TITLE, AUTHOR, YEAR, 
	JOURNAL_TITLE, JOURNAL_ABBREVIATION, VOLUME, ISSUE, PAGE_FROM, PAGE_TO, DOI, 
	ABSTRACT, ISSN, END_OF_REFERENCE, URL, 
	RIS_SOURCE, RIS_COMMENT
}

public class RISParser
{
	static final  Map<String, RISTag> tagMap = new HashMap<String, RISTag>()
	{
		private static final long serialVersionUID = 1L;
		{
			put("ER  -", RISTag.END_OF_REFERENCE);
			put("TY  -", RISTag.CITATION_TYPE);
			put("TI  -", RISTag.CITATION_TITLE);
			put("T1  -", RISTag.CITATION_TITLE);
			put("BT  -", RISTag.CITATION_TITLE);
			put("CT  -", RISTag.CITATION_TITLE);
			put("A1  -", RISTag.AUTHOR);
			put("A2  -", RISTag.AUTHOR);
			put("A3  -", RISTag.AUTHOR);
			put("A4  -", RISTag.AUTHOR);
			put("AU  -", RISTag.AUTHOR);
			put("Y1  -", RISTag.YEAR);
			put("PY  -", RISTag.YEAR);
			put("JF  -", RISTag.JOURNAL_TITLE);
			put("JO  -", RISTag.JOURNAL_TITLE);
			put("T2  -", RISTag.JOURNAL_TITLE);
			put("JA  -", RISTag.JOURNAL_ABBREVIATION);
			put("VL  -", RISTag.VOLUME);
			put("IS  -", RISTag.ISSUE);
			put("CP  -", RISTag.ISSUE);
			put("SP  -", RISTag.PAGE_FROM);
			put("EP  -", RISTag.PAGE_TO);
			put("M3  -", RISTag.DOI);
			put("DO  -", RISTag.DOI);
			put("N1  -", RISTag.RIS_COMMENT);
			put("N2  -", RISTag.RIS_COMMENT);
			put("AB  -", RISTag.ABSTRACT);
			put("SN  -", RISTag.ISSN);
			put("UR  -", RISTag.URL);
			
			put("N1  -", RISTag.RIS_COMMENT);
			put("N2  -", RISTag.RIS_COMMENT);
			put("DB  -", RISTag.RIS_SOURCE);
		}
	};

	Journal jrn;
	Article art;
	
	public RISTag isTag(String piece)
	{
		RISTag tag = tagMap.get(piece);
		if (tag == null)
			if (piece.matches("[A-Z]{2,2}\\s{2,2}-"))
				return RISTag.UNKNOWN;
			else
				return null;
		else
			return tag;
	}
	
	public void processTag(RISTag tag, String value) throws Exception
	{
		switch (tag)
		{
			case UNKNOWN:
				return;
			case CITATION_TYPE:
				//Obsolete
				return;
			case CITATION_TITLE://Citation title
				art.setTitle(value);
				return;
			case AUTHOR://Author
				art.addAuthor(Author.getFromRIS(value));
				return;
			case YEAR://Year
				if(value.trim().equals(""))
					art.setCurrentPublicationDate();
				else
					art.setPublicationDate(value);
				return;
			case JOURNAL_TITLE: //Title
				Journal tmp = Journal.getByTitle(value);
				if (tmp != null)
					jrn = tmp;
				else
				{
					jrn = new Journal();
					jrn.setTitle(value);
					if (Globals.userSession().user != null)
					{
	    	    		jrn.owner = Globals.userSession().user;
						if(jrn.introducer == null)
							jrn.introducer = Globals.userSession().user;
					}
				}
				return;
			case JOURNAL_ABBREVIATION:
				if (jrn.id == null)
					jrn.setAbbreviation(value);
				return;
			case VOLUME://volume
				art.volume = value;
				return;
			case ISSUE:
				art.issue = value;
				return;
			case PAGE_FROM:
				art.pageNumbers = value;
				return;
			case PAGE_TO:
				art.pageNumbers+=("-"+value);
				return;
			case DOI:
				if(value.contains("dx.doi.org"))
					value = value.replaceAll("http://dx.doi.org/", "").trim();
				if(value.contains("doi:"))
					value = value.replaceAll("doi:", "").trim();
				if(value.contains("DOI:"))
					value = value.replaceAll("DOI:", "").trim();
				art.doi = value.trim();
				return;
			case ABSTRACT:
				art.articleAbstract+=(value+" ");
				return;
			case ISSN:
				if (jrn.id == null)
					jrn.setISSN(value);
				return;
			case END_OF_REFERENCE:
				return;
			case URL:
				art.setLink(value);
				return;
			case RIS_COMMENT:
				art.addComment(value + "\n");
				return;
			case RIS_SOURCE:
				art.addComment("Source of export: " + value + "\n");
				return;
			default:
				return;
		}
		
	}
	
	public Article parse(InputStream is) throws Exception
	{
		jrn = new Journal();
		art = new Article();
		art.articleAbstract = "";
		art.setComment("");
		
		BufferedReader in = new BufferedReader(new InputStreamReader(is, "UTF-8"));
		String line;
		String curValue = "";
		RISTag curTag = RISTag.UNKNOWN, oldTag = RISTag.UNKNOWN;
		while ((line = in.readLine()) != null)
		{
			System.out.println(line);
			
			if(line.length() >= 5)
				curTag = isTag(line.substring(0, 5)); 
						
			if (curTag != null)
			{
				processTag(oldTag, curValue);
				oldTag = curTag;
				if(line.length() > 6)
					curValue = line.substring(6);
				else
					curValue = "";
			} else
				curValue+=line;
			
			if (curTag == RISTag.END_OF_REFERENCE)
				break;
		}
		in.close();
		
		Journal journal = Journal.getByISSN(jrn.getISSN());
		if (journal == null)
			journal = Journal.getByTitle(jrn.getTitle());

		// There is no such journal in our database, create one
		if (journal == null)
		{
			journal = new Journal();
			journal.setISSN(jrn.getISSN());
			journal.setTitle(jrn.getTitle());
			journal.setAbbreviation(jrn.getAbbreviation());
			if (Globals.userSession().user != null)
			{
				journal.owner = Globals.userSession().user;
				if (journal.introducer == null)
					journal.introducer = Globals.userSession().user;
			}
			Globals.session().saveOrUpdate(journal);
		}
		
		art.journal = journal;
		//art.journal = jrn;
		return art; // return a newly created article 
	}
	
	public Article parse(File f) throws Exception
	{
		return parse(new FileInputStream(f));
	}
}
