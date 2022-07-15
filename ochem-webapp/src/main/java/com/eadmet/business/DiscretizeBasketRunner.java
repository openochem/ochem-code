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

package com.eadmet.business;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.dao.Repository;
import qspr.entities.Basket;
import qspr.entities.BasketEntry;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Predicate;
import qspr.entities.Property;
import qspr.entities.PropertyOption;
import qspr.entities.UnitCategory;
import qspr.util.WrapperThread;
import qspr.util.unitconversion.UnitConversion;
import qspr.workflow.utils.QSPRConstants;

public class DiscretizeBasketRunner extends WrapperThread
{
	private static transient final Logger logger = LogManager.getLogger(DiscretizeBasketRunner.class);

	DiscretizeBasketOptions opts;

	public DiscretizeBasketRunner(DiscretizeBasketOptions opts)
	{
		this.opts = opts;
	}

	@SuppressWarnings("unchecked")
	@Override
	public void wrapped() throws Exception 
	{
		Basket basket = Basket.getBasket(Globals.userSession(), opts.basketId);

		PropertyOption poptions[] = new PropertyOption[opts.options.length];

		assert opts.options.length == opts.thresholds.length + 1;

		for (int i = 0; i < opts.thresholds.length; i++)
			opts.thresholds[i] = new Double(opts.strThresholds[i]);

		Basket newBasket = Basket.getBasket(opts.newBasketName);
		Property newProperty = Repository.property.getProperty(opts.newPropertyName, true);
		newProperty.type = Property.TYPE_QUALITATIVE;
		newProperty.unitCategory = UnitCategory.getByName(QSPRConstants.CLASS);
		newProperty.defaultUnit = newProperty.unitCategory.getDefaultUnit();

		Property oldProperty = Property.getById(opts.propertyId);
		if (!oldProperty.obligatoryConditions.isEmpty())
			for (Property condition : oldProperty.obligatoryConditions)
				newProperty.obligatoryConditions.add(condition);

		int i = 0;
		for (String option : opts.options) 
		{
			PropertyOption poption = newProperty.getOptionByName(option);
			if (poption == null)
			{
				if(newProperty.id == null)Globals.session().saveOrUpdate(newProperty);
				poption = Repository.option.getPropertyOptionByName(option, newProperty.id, true, false);
				newProperty.options.add(poption);
				Globals.session().saveOrUpdate(newProperty);
			}
			poptions[i++] = poption;
		}
		Globals.session().saveOrUpdate(newProperty);


		int dublicates = 0;
		List<Long> beIDs = Globals.session().createCriteria(BasketEntry.class).setProjection(Projections.groupProperty("id")).add(Restrictions.eq("basket", basket)).list();
		int batchSize = 100;
		int cnt = 0;
		Predicate equals = Predicate.get("=");

		Set<String> md5 = new HashSet<String>();
		for(BasketEntry entry: newBasket.entries)
			md5.add(entry.ep.md5); 

		// Fetch and create the entries in batches
		for (int k = 0; k < beIDs.size(); k += batchSize)
		{
			if (cancelRequested)
				break;
			List<Long> batchIDs = beIDs.subList(k, Math.min(k + batchSize, beIDs.size()));
			List<BasketEntry> batchEntries = Globals.session().createCriteria(BasketEntry.class).add(Restrictions.in("id", batchIDs)).list();

			//int cnt = 0;
			for (BasketEntry be : batchEntries)
			{
				if (cancelRequested)
					break;
				setStatus("Processed " + (cnt++) + " out of " + beIDs.size() + (dublicates > 0 ? "\n ("+dublicates+" dublicates skipped)" : ""));
				logger.info("Processed " + ((cnt)) + " out of " + beIDs.size() + (dublicates > 0 ? "\n ("+dublicates+" dublicates skipped)" : ""));
				ExperimentalProperty ep = be.ep;
				ExperimentalProperty newEp = be.ep.cloneForReference();
				newEp.moleculenames = null;
				newEp.article = ep.article;
				newEp.property = newProperty;
				newEp.introducer = newEp.owner = Globals.userSession().user;
				newEp.rights = Globals.RIGHTS_NONE; // created by default as hidden
				newEp.unit = newEp.property.defaultUnit;
				newEp.connectedProperty = be.ep;
				newEp.predicate = equals;
				newEp.externalId = be.ep.externalId;
				newEp.conditions = be.ep.conditions;
				newEp.artLineNum = ep.artLineNum;
				newEp.artMolId = ep.artMolId;
				newEp.artPageNum = ep.artPageNum;
				newEp.artParagraph = ep.artParagraph;

				double value = UnitConversion.convert(ep.value, ep.unit, ep.property.defaultUnit, ep.molecule.molWeight);
				i = 0;
				while (i < opts.thresholds.length && value > opts.thresholds[i])
					i++;
				newEp.option = poptions[i];
				newEp.other = "Converted from a numerical value using the condition: ";
				if (i == 0)
					newEp.other += "" + value + " < " + opts.thresholds[0];
				else if (i == opts.thresholds.length)
					newEp.other += "" + value + " > " + opts.thresholds[opts.thresholds.length - 1];
				else
					newEp.other += "" + opts.thresholds[i - 1] + " < " + value + " < " + opts.thresholds[i];
				newEp.other += " " + ep.property.defaultUnit.getName();

				newEp.updateHash();
				if (!newEp.hasConflicts())
					Globals.session().save(newEp);
				else
				{
					logger.info("A record already exists, using it.");
					dublicates++;
					if(md5.contains(newEp.md5))continue; // it was already stored 
					newEp = newEp.duplicate;
				}

				md5.add(newEp.md5);

				BasketEntry newBe = new BasketEntry(newEp);
				newBe.basket = newBasket;
				newBe.ep = newEp;
				try{
					Globals.session().saveOrUpdate(newBe); // could be already there, just ignore exception
				}catch(Exception e){

				}
			}

			//setStatus("Processed " + k + " out of " + basket.entries.size());
			Globals.restartAllTransactions(true);
			equals = Predicate.get("=");
		}
		setStatus("Finished");
	}

}