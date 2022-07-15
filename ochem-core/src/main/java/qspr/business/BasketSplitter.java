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

package qspr.business;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import qspr.Globals;
import qspr.entities.Basket;
import qspr.util.WrapperThread;

/**
 * An abstraction of a splitter of a basket into two sets.
 * @author midnighter
 *
 */

public abstract class BasketSplitter extends WrapperThread
{
	public Basket basket;
	public Basket set1;
	public Basket set2;

	public String basketName1;
	public String basketName2;

	protected abstract List<Integer> getSplittedIndices(List<Long> epIDs);

	public void split()
	{
		basket = Basket.getBasket(Globals.userSession(), basket.id);

		setStatus("Fetching basket entries");
		List<Long> epIDs = basket.getRecordIDs();

		Collections.shuffle(epIDs); // random shuffling to avoid ordering effect

		setStatus("Performing the split");
		List<Integer> indices = getSplittedIndices(epIDs);

		if (basketName1 == null)
			basketName1 = basket.name + " (training)";
		if (basketName2 == null)
			basketName2 = basket.name + " (test)";

		set1 = new Basket();
		set1.name = Basket.getFreeBasketName(basketName1);
		set2 = new Basket();
		set2.name = Basket.getFreeBasketName(basketName2);

		set1.session = Globals.userSession();
		set2.session = Globals.userSession();
		set1.user = Globals.userSession().user;
		set2.user = Globals.userSession().user;

		Globals.session().save(set1);
		Globals.session().save(set2);

		Collections.sort(indices);

		List<Long> selectedEPIDs = new ArrayList<Long>();
		for (int i = indices.size() - 1; i >= 0; i--)
		{
			selectedEPIDs.add(epIDs.get(indices.get(i)));
			epIDs.remove(0 + indices.get(i));
		}

		setStatus("Adding entries to set 1");
		// Add the seleceted entries to set 1
		int i = 0;
		while (i < selectedEPIDs.size())
		{
			int k = Math.min(selectedEPIDs.size(), i + 2000);
			List<Long> batchEPIDs = selectedEPIDs.subList(i, k);
			i = k;
			Globals.session().createSQLQuery("insert into BasketEntry(basket_id, exp_property_id, exclude) select " + set1.id + ", exp_property_id, exclude from BasketEntry where basket_id=:basketID and exp_property_id in (:epIDs)")
			.setParameter("basketID", basket.id)
			.setParameterList("epIDs", batchEPIDs)
			.executeUpdate();
			setStatus("Adding entries to set 1 - " + i + " records added");
		}

		// Add the remaining entries to set 2
		setStatus("Adding entries to set 2");
		i = 0;
		while (i < epIDs.size())
		{
			int k = Math.min(epIDs.size(), i + 2000);
			List<Long> batchEPIDs = epIDs.subList(i, k);
			i = k;
			Globals.session().createSQLQuery("insert into BasketEntry(basket_id, exp_property_id, exclude) select " + set2.id + ", exp_property_id, exclude from BasketEntry where basket_id=:basketID and exp_property_id in (:epIDs)")
			.setParameter("basketID", basket.id)
			.setParameterList("epIDs", batchEPIDs)
			.executeUpdate();
			setStatus("Adding entries to set 2 - " + i + " records added");
		}

		setStatus("Finished");

	}

	@Override
	public void wrapped()
	{
		split();
	}
}


