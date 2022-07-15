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

package com.eadmet.batchupload.main;

import java.lang.annotation.Retention;
import java.lang.annotation.Target;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;

import org.hibernate.Hibernate;

import qspr.Globals;
import qspr.dao.Repository;
import qspr.entities.Article;
import qspr.entities.Basket;
import qspr.entities.BasketEntry;
import qspr.entities.ConditionSet;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Predicate;
import qspr.entities.Property;
import qspr.entities.PropertyValue;
import qspr.entities.Unit;
import qspr.metaserver.util.MixtureAttachment;
import qspr.util.CASRN;
import qspr.util.MoleculePeer;
import qspr.util.UploadContext;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.batchupload.main.RecordStub.ColumnValue;
import com.eadmet.batchupload.main.RecordStub.NameValueStub;
import com.eadmet.batchupload.main.RecordStub.ValueStub;
import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.OCHEMUtils;

public class RecordBuilder 
{
	UploadContext context;

	private ExperimentalProperty ep;
	private RecordStub stub;

	public RecordBuilder(UploadContext context)
	{
		this.context = context;
	}

	public void build(RecordPreview preview, RecordStub stub)
	{
		this.stub = stub;

		ep = new ExperimentalProperty();
		ep.id = preview.id;

		int maxOrder = -1;
		for (Method m :  RecordBuilder.class.getDeclaredMethods()) 
			if (m.isAnnotationPresent(Setter.class))
				maxOrder = Math.max(maxOrder, m.getAnnotation(Setter.class).order());

		for (int currentOrder = 1; currentOrder <= maxOrder; currentOrder++)
			for (Method m :  RecordBuilder.class.getDeclaredMethods()) 
			{
				if (!m.isAnnotationPresent(Setter.class))
					continue;
				if (m.getAnnotation(Setter.class).order() != currentOrder)
					continue;

				try
				{
					m.invoke(this);
				}
				catch (InvocationTargetException e)
				{
					preview.messages.add(BatchUploadMessage.newError(e.getTargetException()));
				}
				catch (Exception e)
				{
					preview.messages.add(BatchUploadMessage.newError(e));
				}
			}
		preview.ep = ep;
	}

	@Setter(order=1)
	void setMolecule() throws Exception {
		ep.molecule = null;

		try {

			if (stub.molecules.size() > 0) // first try to get molecule from structure
				for(ColumnValue mol :stub.molecules)

					switch(mol.column.name.toUpperCase()) {

					case  QSPRConstants.MIXTURE_CONDITION:
						MixtureAttachment at = ExperimentalProperty.createMixtureAttachment(mol.value);
						ep.molecule = MoleculePeer.fetchFromString(at.smiles(), context);
						break;

					default: 
						ep.molecule = MoleculePeer.fetchFromString(mol.value, context);
						if(ep.molecule != null) break; // the first found, no further processing
					}

			if(ep.molecule == null && stub.name.size() > 0) // trying to get molecule from name
				for(ColumnValue mol :stub.name){
					ep.molecule = MoleculePeer.fetchFromString(mol.value, context);
					if(ep.molecule != null){
						if(ep.molecule.mapping1.id != 1l) // not empty one!
							break;  // molecule was found!
					}
				}

			if(ep.molecule == null){
				ep.molecule = Repository.molecule.getEmptyMolecule(); // to prevent failure for empty molecules

				if(stub.molecules.size() > 0 && stub.name.size() > 0)
					throw new Exception("Could not create a molecule from the provided molecular structure and names; empty molecule is used");

				if(stub.molecules.size() > 0)
					throw new Exception("Could not create a molecule from the provided molecular structure; empty molecule is used");

				if(stub.name.size() > 0)
					throw new Exception("Could not create a molecule from the provided names; empty molecule is used");
			}

		}
		catch(Exception e) { // should not be but ...
			if(ep.molecule == null || ep.molecule.id == null) ep.molecule = Repository.molecule.getEmptyMolecule();
			throw e;
		}

	}

	@Setter(order=2)
	void setProperty() throws Exception
	{
		if (stub.property.value == null || stub.property.name == null)
			throw new Exception("No property set for this record");
		PropertyValue pv = new PropertyValue();
		isetProperty(pv, stub.property);
		if (pv.property.isTextual())
			throw new Exception("Textual properties are not allowed");
		ep.property = pv.property;
	}

	@Setter(order=3)
	void setValue() throws Exception
	{
		PropertyValue pv = new PropertyValue();
		pv.property = ep.property;
		isetValue(pv, stub.property.value.getFirst());
		ep.predicate = pv.predicate;
		if (ep.property.isNumeric())
		{
			ep.value = pv.value;
			ep.secondValue = pv.secondValue;
		} else
			ep.option = pv.option;
	}

	@Setter(order=4)
	void setSecondValue() throws Exception
	{
		PropertyValue pv = new PropertyValue();
		pv.property = ep.property;
		isetSecondValue(pv, stub.property.value.getFirst());
		if (pv.predicate != null)
			ep.predicate = pv.predicate;
		if (pv.secondValue != null)
			ep.secondValue = pv.secondValue;
	}

	@Setter(order=5)
	void setConditions() throws Exception
	{

		ConditionSet conditions = null;

		if (stub.property.conditions.size() != 0) {
			conditions = new ConditionSet();
			for (NameValueStub nvs : stub.property.conditions) 
			{
				PropertyValue pv = new PropertyValue();
				isetProperty(pv, nvs);
				isetValue(pv, nvs.value.getFirst());
				isetSecondValue(pv, nvs.value.getFirst());
				isetUnit(pv, nvs.value.getFirst());
				conditions.values.add(pv);
			}
		}

		for(ColumnValue mol :stub.molecules)
			switch(mol.column.name.toUpperCase()) {
			case  QSPRConstants.MIXTURE_CONDITION:

				MixtureAttachment at = ExperimentalProperty.createMixtureAttachment(mol.value);
				if(at.fractions == null) break;
				if(conditions ==null)conditions = new ConditionSet();
				conditions.addMixtures(at.toString());
				break;
			}

		ep.conditions = conditions == null ? null: conditions.get();
	}

	@Setter(order=6)
	void setUnit() throws Exception
	{
		PropertyValue pv = new PropertyValue();
		pv.property = ep.property;
		isetUnit(pv, stub.property.value.getFirst());
		ep.unit = pv.unit;
	}

	@Setter(order=7)
	void setArticle() throws Exception
	{
		if (stub.article.size() > 0)
		{
			if (Globals.userSession().user == null)
				throw new UserFriendlyException("Creation of article is only allowed to registered users");

			String name = stub.article.getFirst().value;
			if (name.equals("unpublished")) // upload of data from a basket
				ep.article = Article.getDefaultArticle(Globals.userSession().user, context);
			else
				ep.article = Article.getArticle(name, context);

			if (ep.article == null)
				throw new Exception("Could not find article "+stub.article.getFirst().value);
		}
		else
			throw new Exception("No article found");
	}

	@Setter(order=8)
	void setBasket() throws Exception
	{
		if (ep.basketEntries == null)
			ep.basketEntries = new ArrayList<BasketEntry>();
		else
			ep.basketEntries.clear();

		for (ColumnValue basket : stub.basket) 
		{
			String name = basket.value;
			if (name.equals("default"))
				name = OCHEMUtils.getFilteredBasketName(context.fileName);
			Basket b = Basket.getBasket(name, context);
			Hibernate.initialize(b.entries);

			BasketEntry be = new BasketEntry();
			be.ep = ep;
			be.basket = b;

			ep.basketEntries.add(be);
			b.cachedCount = null;
		} 
	}

	@Setter(order=9)
	void setNames() throws Exception
	{
		for (ColumnValue cv : stub.casrn) 
		{
			String value = CASRN.checkCasrnSyntax(cv.value);
			if (value != null)
				ep.addName(value);
		}

		for (ColumnValue cv : stub.name)
			ep.addName(cv.value);
	}

	@Setter(order=10)
	void setComments() throws Exception
	{
		for (ColumnValue cv : stub.comments)
		{
			if (ep.other == null)
				ep.other = cv.value;
			else
				ep.other += ("\t" + cv.value);
		}
	}

	@Setter(order=11)
	void setAttributes()
	{
		if (stub.hidden.size() > 0)
			if (Boolean.valueOf(stub.hidden.getFirst().value))
				//				if (Globals.userSession().user != null)
				ep.rights = Globals.RIGHTS_NONE;

		if (context.hiddenByDefault)
			//			if (Globals.userSession().user != null)
			ep.rights = Globals.RIGHTS_NONE;

		if (stub.line.size() > 0)
			ep.artLineNum = Double.valueOf(stub.line.getFirst().value).intValue();

		if (stub.n.size() > 0)
			ep.artMolId = stub.n.getFirst().value;

		if (stub.externalid.size() > 0)
			ep.externalId = stub.externalid.getFirst().value;

		if (stub.page.size() > 0)
			ep.artPageNum = Double.valueOf(stub.page.getFirst().value).intValue();

		if (stub.table.size() > 0)
			ep.artTableNum = stub.table.getFirst().value;

		if (stub.evidence.size() > 0)
			if ("measured".equals(stub.evidence.getFirst().value))
				ep.connectedProperty = ep;
	}

	@Setter(order=12)
	void setOwner()
	{
		ep.owner = Globals.userSession().user;

		if (ep.id == null || ep.id < 0)
			ep.introducer = Globals.userSession().user;
	}

	// Reused for ExperimentalProperty and Condition
	private void isetProperty(PropertyValue pv, NameValueStub stub)
	{
		pv.property = Property.getByName(stub.name.value);
	}

	// Reused for ExperimentalProperty and Condition
	private void isetValue(PropertyValue pv, ValueStub stub) throws Exception
	{
		String value = stub.value.value;
		pv.predicate = Predicate.get("=");
		if (pv.property.isNumeric())
		{
			pv.setValueWithPredicate(value);

			if (pv.value == null)
				throw new Exception("Value "+value+" could not be properly set for property "+pv.property.getName());

			if (stub.predicate.size() > 0)
				pv.predicate = Predicate.get(stub.predicate.get(0).value);
		}
		else if (pv.property.isQualitative())
		{
			pv.option = pv.property.getOptionByName(value);
			if (pv.option == null)
				throw new Exception("Option "+value+" not found for property "+pv.property.getName());
		}
		else
			pv.textualValue = value;
	}

	// Reused for ExperimentalProperty and Condition
	private void isetUnit(PropertyValue pv, ValueStub stub) throws Exception
	{
		if (!pv.property.isNumeric())
			return;

		if (stub.unit.size() > 0)
		{
			String unit = stub.unit.getFirst().value;
			if (unit.equals("default"))
				pv.unit = pv.property.defaultUnit;
			else
				pv.unit = Unit.getByNameAndCategory(unit, pv.property.defaultUnit.category.name, false);
			if (pv.unit == null)
				throw new Exception("Unit" +unit+" not found for property "+pv.property.getName());
		}
		else
			throw new Exception("No unit found");
	}

	// Reused for ExperimentalProperty and Condition
	private void isetSecondValue(PropertyValue pv, ValueStub stub) throws Exception
	{
		if (pv.property.isQualitative())
			return;

		if (stub.accuracy.size() > 0 && stub.interval.size() > 0)
			throw new Exception("Value can not have both accuracy and interval columns set simultaneously");

		if (stub.accuracy.size() > 0)
		{
			ep.secondValue = Double.valueOf(stub.accuracy.getFirst().value);
			ep.predicate = Predicate.get("+-");
		} else
			if (stub.interval.size() > 0)
			{
				ep.secondValue = Double.valueOf(stub.interval.getFirst().value);
				ep.predicate = Predicate.get("-");
			}
	}

}

@Retention(value=java.lang.annotation.RetentionPolicy.RUNTIME)
@Target(value={java.lang.annotation.ElementType.METHOD})
abstract @interface Setter {
	public abstract int order() default 1;
}

