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

package qspr.schedule;

import java.sql.Timestamp;
import java.util.Calendar;
import java.util.List;

import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.entities.Article;
@DatabaseMaintenanceJob
public class CleanArticleTask extends OchemCronjobTask
{

	@SuppressWarnings("unchecked")
	public void executeTask() throws Exception 
	{
		Globals.startMainTransaction();
		List<Article> articlecands = Globals.session().createCriteria(Article.class).add(Restrictions.isNotNull("pmid")).list();
		for (Article article : articlecands)
		{
			long time = Calendar.getInstance().getTimeInMillis() - 48 * 60 * 60 * 1000; // 48h
			Timestamp ts = new Timestamp(time);
			if (article.time != null) 
				if ((article.getRecordSizeDummy() == 0) && (article.getRecordSizeNondummy() == 0) && (ts.compareTo(article.time) >= 1)) 
				{
					log("Deleting an article: A" + article.id);
					Globals.session().delete(article);
				}
		}
		Globals.commitMainTransaction();
	}

}
