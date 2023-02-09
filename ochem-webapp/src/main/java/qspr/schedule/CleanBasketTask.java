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

import java.util.List;

import qspr.Globals;

@DatabaseMaintenanceJob
public class CleanBasketTask extends OchemCronjobTask {

	@SuppressWarnings("unchecked")
	public void executeTask() throws Exception 
	{
		Globals.startMainTransaction();
		// NoS 12.01.10 - Delete anonymous baskets older then 24 hours and ->not involved in models<-.
		// TODO: Should we allow anonymous baskets to be involved in model building? Now it's allowed 

		String sqlSelect = 
				"select b.basket_id from Basket b left join Model m1 on (m1.training_set_id = b.basket_id) "+
						"left join Model m2 on (m2.validation_set_id = b.basket_id) "+
						"left join ValidationSet vs using (basket_id) "+
						"left join Session s on (s.session_id = b.session_id) "+
						"where m1.model_id is null " +
						"and m2.model_id is null " +
						"and b.user_id is null " +
						"and vs.validation_set_id is null "+
						"and timestampdiff(HOUR, s.activity_time, now()) > 24";

		List<Integer> ids = Globals.session().createSQLQuery(sqlSelect).list();

		for (Number id : ids)
		{
			log("Basket deleted: "+id);
			Globals.session().createSQLQuery("delete from BasketEntry where basket_id = :basket_id").setInteger("basket_id", id.intValue()).executeUpdate();
			Globals.session().createSQLQuery("delete from Basket where basket_id = :basket_id").setInteger("basket_id", id.intValue()).executeUpdate();
		}

		Globals.commitMainTransaction();
	}

	public static void main(String[] args)
	{
		new CleanBasketTask().executeInternal();
	}
}
