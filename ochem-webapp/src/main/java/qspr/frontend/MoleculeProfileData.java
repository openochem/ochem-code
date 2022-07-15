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

package qspr.frontend;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.annotation.XmlRootElement;

import org.hibernate.Criteria;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.business.BasketFilter;
import qspr.business.BasketPeer;
import qspr.business.ModelsFilter;
import qspr.dao.MetalBondParserSdf;
import qspr.dao.Various;
import qspr.entities.Basket;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Mapping2;
import qspr.entities.Molecule;
import qspr.entities.Property;
import qspr.entities.Tag;

import com.eadmet.exceptions.UserFriendlyException;

/**
 * Data for the molecule profile page
 * @author midnighter
 *
 */
@XmlRootElement(name = "molecule-profile")
public class MoleculeProfileData
{
	public Molecule molecule;
	public String canonicalName;
	public String smiles;
	public String inchi;

	public List<Basket> baskets;

	public List<Tag> tags;

	public int recordsCount;

	public long modelsCount;
	public long pendingModelsCount;
	public String formula;

	public List<Property> properties;

	@SuppressWarnings("unchecked")
	public MoleculeProfileData(Mapping2 mapping2) throws IOException
	{
		molecule = mapping2.getMolecule();

		if (molecule.isEmptyMolecule())
			throw new UserFriendlyException("This is an empty molecule. The molecule profile is unavailable.");

		String mol = MetalBondParserSdf.eliminateMetalBond(molecule.getData());
		canonicalName = Various.molecule.convertToCanonicalName(mol);
		smiles = Various.molecule.convertToSmilesOrSmart(mol);
		inchi = mapping2.inchi2;
		formula = Various.molecule.getFormula(mol);

		// Baskets
		BasketFilter filter = new BasketFilter();
		filter.mapping2 = mapping2;
		baskets = BasketPeer.getListCriteria(filter).list();

		// Models
		if (!baskets.isEmpty())
		{
			ModelsFilter modelFilter = new ModelsFilter();
			modelFilter.trainingSets = baskets;
			modelsCount = (Long) modelFilter.createCriteria().setProjection(Projections.countDistinct("id")).uniqueResult();

			modelFilter.pendingTasks = true;
			pendingModelsCount = (Long) modelFilter.createCriteria().setProjection(Projections.countDistinct("id")).uniqueResult();
		}

		// Properties
		getProperties();
	}

	@SuppressWarnings("unchecked")
	private void getProperties()
	{
		properties = new ArrayList<Property>();
		recordsCount = 0;

		Criteria c = Globals.session().createCriteria(ExperimentalProperty.class).createAlias("molecule", "mol")
				.add(Restrictions.eq("mol.mapping2", molecule.mapping2));
		ExperimentalProperty.addAccessRestrictions(c, Globals.RIGHTS_FREELY_AVAILABLE, Globals.userSession().user, false, true);
		c.setProjection(Projections.projectionList().add(Projections.groupProperty("property")).add(Projections.countDistinct("id")));
		List<Object[]> rows = c.list();

		for (Object[] row : rows)
		{
			Property property = (Property) row[0];
			property.count = (Long) row[1];
			properties.add(property);
			recordsCount += property.count;
		}
	}

	// Make JAXB happy
	public MoleculeProfileData()
	{

	}
}
