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

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import qspr.entities.Article;
import qspr.entities.Basket;
import qspr.entities.Molecule;

@XmlRootElement
public class UploadContext
{
	//Settings
	@XmlAttribute
	public boolean allowMoleculePubchemSearch = true;
	@XmlAttribute
	public boolean allowArticlePubmedSearch = true;
	@XmlAttribute
	public boolean hiddenByDefault = true;
	//
	public boolean preview = false;
	//
	public String fileName;
	//
	public boolean interruptRequested = false;
	// Cache
	@XmlTransient
	public CacheMap<Molecule> moleculeCache = new CacheMap<Molecule>();
	@XmlTransient
	public CacheMap<String> recordHashCache = new CacheMap<String>(); //Record MD5 -> Some numerical identifier of a record (i.e. row number in a sheet)
	@XmlTransient
	public CacheMap<Article> articleCache = new CacheMap<Article>();
	@XmlTransient
	public CacheMap<Basket> basketCache = new CacheMap<Basket>();

	public void clearCaches(boolean clearInternalDublicatesCache)
	{
		moleculeCache.clear();
		articleCache.clear();
		basketCache.clear();
		
		if (clearInternalDublicatesCache)
			recordHashCache.clear();
	}
}