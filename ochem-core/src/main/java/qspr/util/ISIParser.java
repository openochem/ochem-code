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

import qspr.Globals;
import qspr.entities.Article;
import qspr.entities.Author;
import qspr.entities.Journal;

// TODO: Consider cases of books

public class ISIParser
{

	Journal jrn;
	
	Article art;
	
	public int isTag(String piece)
	{
		if (piece.equals("ER "))
			return -2;
		else if (piece.equals("TI ") || piece.equals("T1 ") || piece.equals("BT ") || piece.equals("CT "))
			return 1;
		else if (piece.equals("A1 ") || piece.equals("AU "))
			return 2;
		else if (piece.equals("Y1 ") || piece.equals("PY "))
			return 3;
		else if (piece.equals("SO "))
			return 4;
		else if (piece.equals("JA ") || piece.equals("JI "))
			return 5;
		else if (piece.equals("VL "))
			return 6;
		else if (piece.equals("IS ") || piece.equals("CP "))
			return 7;
		else if (piece.equals("BP ") || piece.equals("SP "))
			return 8;
		else if (piece.equals("EP "))
			return 9;
		else if (piece.equals("M3 ") || piece.equals("DI "))
			return 10;	
		else if (piece.equals("N1 ") || piece.equals("N2 ") || piece.equals("AB "))
			return 11;	
		else if (piece.equals("SN "))
			return 12;
		else if (piece.equals("PU "))
			return 13;
		else if (piece.equals("   "))
			return 0;
		return -1;
	}
	
	public void processTag(int tag, String value) throws Exception
	{
		switch (tag)
		{
			case -1:
			case  0: return;
			case 1://Citation title
				art.setTitle(value.replaceAll("#@#", ""));
				return;
			case 2://Author
				if(!value.trim().equals(""))
				{
					List<Author> auth = Author.getFromISI(value.trim());
					for (Author author : auth) {
						art.addAuthor(author);
					}	
				}
				return;
			case 3://Year
				if (!value.equals(""))
					art.setPublicationDate(value.replaceAll("#@#", ""));
				//RESOLVE DATE ISSUES
				return;
			case 4: //Title
				Journal tmp = Journal.getByTitle(value.replaceAll("#@#", ""));
				if (tmp != null)
					jrn = tmp;
				else
				{
					jrn = new Journal();
					jrn.setTitle(value.replaceAll("#@#", ""));
					if (Globals.userSession().user != null)
					{
	    	    		jrn.owner = Globals.userSession().user;
						if(jrn.introducer == null)
							jrn.introducer = Globals.userSession().user;
					}
				}
				return;
			case 5://no more in isi parser
				if (jrn.id == null)
				{
					jrn.setAbbreviation(value.replaceAll("#@#", ""));
				}
				return;
			case 6://volume
				art.volume = value.replaceAll("#@#", "");
				return;
			case 7:
				art.issue = value.replaceAll("#@#", "");
				return;
			case 8:
				art.pageNumbers = value.replaceAll("#@#", "");
				return;
			case 9:
				art.pageNumbers+=("-"+value.replaceAll("#@#", ""));
				return;
			case 10:
				if (!value.equals(""))
					art.doi = value.replaceAll("#@#", "").trim();
				return;
			case 11:
				art.articleAbstract+=(value.replaceAll("#@#", "")+" ");
				return;
			case 12:
				if (jrn.id == null)
					jrn.setISSN(value.replaceAll("#@#|-", ""));
				return;
			case 13:
				if (jrn.id == null)
					jrn.publisher = value.replaceAll("#@#", "").trim();
		}
		
	}
	
	public Article parse(File f) throws Exception
	{
		FileInputStream fis = new FileInputStream(f);
		BufferedReader in = new BufferedReader(new InputStreamReader(fis));
		String line;
		String curValue = "";
		int curTag = 0, oldTag = 0;
		
		jrn = new Journal();
		art = new Article();
		art.articleAbstract = "";
		while ((line = in.readLine()) != null)
		{
			if(line.length() >= 3)
				curTag = isTag(line.substring(0, 3)); 
						
				if (curTag != 0)
				{
					processTag(oldTag, curValue);
					oldTag = curTag;
					if(line.length() > 3)
						curValue = line.substring(3).trim();
					else
						curValue = "";
				}else
				{
					curValue+="#@# "+line.trim();
				}
				
				if (curTag == -2)
				{
					break;
				}
		}
		art.journal = jrn;
		Globals.session().save(jrn);
		in.close();
		
		return art; // return a newly created article 
	}
}
