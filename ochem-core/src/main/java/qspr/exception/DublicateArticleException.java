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

package qspr.exception;

import qspr.entities.Article;

import com.eadmet.exceptions.UserFriendlyException;

public class DublicateArticleException extends UserFriendlyException
{
	private static final long serialVersionUID = 1L;
	
	private Article article;
	
	public DublicateArticleException(Article art)
	{
		super("This article is duplicate of another article introduced by " + art.introducer + " (Article identifier A" + art.id + ").");
		this.article = art;
	}
	
	public Object getObject() {
		return article.id;
	}
}
