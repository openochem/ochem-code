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
import java.util.List;

import qspr.Globals;
import qspr.entities.Article;
import qspr.entities.Author;
import qspr.entities.Journal;

// TODO: Consider cases of books


public class EndNoteParser
{

	Journal jrn;
	
	Article art;
	
	public int isTag(String piece)
	{
//		if (piece.equals("ER  - "))
//			return -2;
		if (piece.equals("%0 "))
			return 1;		
		if (piece.equals("%T "))
			return 2;
		if (piece.equals("%A "))
			return 3;
		if (piece.equals("%D "))
			return 4;
		if (piece.equals("%J ") || piece.equals("%B "))
			return 5;
		if (piece.equals("%V "))
			return 7;
		if (piece.equals("%N "))
			return 8;
		if (piece.equals("%P "))
			return 9;
		if (piece.equals("%R ") || piece.equals("%1 "))
			return 11;	
		if (piece.equals("%X "))
			return 12;
		if (piece.equals("%+ "))
			return 13;		
		if (piece.startsWith("%"))
			return -1;
		return 0;
	}
	
	public void processTag(int tag, String value) throws Exception
	{
		switch (tag)
		{
			case -1:
			case  0: return;
			case 1: //Citation type
				return;
			case 2://Citation title
				art.setTitle(value);
				return;
			case 3://Author
				if(!value.trim().equals(""))
				{
					List<Author> auth = Author.getFromEN(value.trim());
					for (Author author : auth) {
						art.addAuthor(author);
					}	
				}
				return;
			case 4://publication date
				if(value.trim().equals(""))
					art.setCurrentPublicationDate();
				else
					art.setPublicationDate(value);
				return;
			case 5: //Title
//				Journal tmp = Journal.getByTitle(value);
//				if (tmp != null)
//					jrn = tmp;
//				else
//				{
//					jrn = new Journal();
					jrn.setTitle(value);
					if (Globals.userSession().user != null)
					{
						jrn.owner = Globals.userSession().user;
						if(jrn.introducer == null)
							jrn.introducer = Globals.userSession().user;
					}
//				}
				return;
			case 7://volume
				art.volume = value;
				return;
			case 8:
				art.issue = value;
				return;
			case 9:
				art.pageNumbers = value;
				return;
			case 11:
				art.doi = value.replaceAll("doi:", "");
				return;
			case 12:
				art.articleAbstract=(value);
				return;
			case 13:
				art.affiliation = value;
				return;	
		}
		
	}
	
	public Article parse(InputStream fis) throws Exception
	{
		BufferedReader in = new BufferedReader(new InputStreamReader(fis));
		String line;
		int curTag = 0;
		String curValue;

		
		jrn = new Journal();
		art = new Article();

		while ((line = in.readLine()) != null)
		{
			int pos = line.indexOf("%");
			if (pos!=-1)
			{	
				line = line.substring(pos, line.length());
				if (line.length() > 3)
				{
					curTag = isTag(line.substring(0, 3)); 
					if (curTag != 0)
					{
						curValue = line.substring(3);
						processTag(curTag, curValue);
					}
				}
			}
		}
		art.journal = jrn;
		//Globals.session().saveOrUpdate(jrn);
		in.close();
		
		return art; // return a newly created article
	}
	
	public Article parse(File f) throws Exception
	{
		return parse(new FileInputStream(f));
	}
}
