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

package qspr.dao;

import java.io.IOException;
import java.util.List;

import org.hibernate.Query;
import org.hibernate.type.IntegerType;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.entities.ColoredName;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Molecule;

public class RecordDAOImp implements RecordDAO {

	@Override
	public ExperimentalProperty getRecord(Long recordId) {
		return (ExperimentalProperty)Globals.session().get(ExperimentalProperty.class, recordId);
	}

	@Override
	public double getWeight(Long recordId) throws IOException {
		Molecule m = ((ExperimentalProperty) Globals.session().get(ExperimentalProperty.class, recordId)).molecule;
		return m == null?0:m.getMolWeight();
	}

	@Override
	public String getSolventName(int id2) {
		System.out.println("Mapping2: " + id2);
		Query q = Globals.session().createQuery(" from ExperimentalProperty ep where ep.article.id = :ar and ep.molecule.mapping2.id = :id ")
				.setParameter("ar", OCHEMConfiguration.solvent).setParameter("id", id2);
		@SuppressWarnings("unchecked")
		List<ExperimentalProperty> entities = q.list();
		if(entities.size() == 0) return null;

		Query q1 = Globals.session().createQuery(" from ColoredName c where ep.id = :id")
				.setParameter("id", entities.get(0).id);

		@SuppressWarnings("unchecked")
		List<ColoredName> names = q1.list();

		if(names.size() >= 1)
			return names.get(0).name;
		return null;
	}

	@Override
	public int getUnpublishedRecordsCountForUser(Long user) {
		String query = "select count(*) c from ExperimentalProperty where " +
				" introducer_id = "+user + " and deleted is null and " +
				" rights != " + Globals.RIGHTS_FREELY_AVAILABLE;
		System.out.println(query);
		
		return  (Integer)Globals.session().createSQLQuery(query)
				.addScalar("c", IntegerType.INSTANCE)
				.uniqueResult();
	}
}
