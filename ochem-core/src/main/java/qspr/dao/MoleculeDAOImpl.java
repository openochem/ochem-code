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

import org.hibernate.Criteria;
import org.hibernate.Query;
import org.hibernate.criterion.Restrictions;

import com.eadmet.utils.OCHEMUtils;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.entities.ColoredName;
import qspr.entities.Mapping1;
import qspr.entities.Mapping2;
import qspr.entities.Molecule;
import qspr.entities.MoleculeFormat;
import qspr.entities.OriginalMolecule;
import qspr.workflow.utils.QSPRConstants;

public class MoleculeDAOImpl implements MoleculeDAO {

	@Override
	public OriginalMolecule getOriginalMoleculeByStructure(String data) {

		data = OriginalMolecule.removeTimestamps(data); //Avoiding  the same molecule with different timestamps
		String md5 = OCHEMUtils.getMD5(data);

		return
				(OriginalMolecule) Globals.session().createCriteria(OriginalMolecule.class)
				.add(Restrictions.eq("md5", md5))
				.setMaxResults(1).uniqueResult();
	}

	@Override
	public Molecule getMoleculeByMD5(String md5) {
		return
				(Molecule) Globals.session().createCriteria(Molecule.class).add(Restrictions.eq("md5", md5)).setMaxResults(1).uniqueResult();
	}

	@Override
	public Molecule getMoleculeByInvariants(String picMd5, String fullInchi2) {
		return 
				(Molecule) Globals.session().createCriteria(Molecule.class)
				.add(Restrictions.eq("pictureMd5", picMd5))
				.createCriteria("mapping2")
				.add(Restrictions.eq("inchi2", fullInchi2))
				.setMaxResults(1).uniqueResult();
	}

	@SuppressWarnings("unchecked")
	public Molecule getEmptyMolecule()
	{
		List<Molecule> list = Globals.session().createCriteria(Molecule.class)
				.add(Restrictions.eq("md5", QSPRConstants.EMPTY_MD5))
				.list();
		if (list.size() > 0)
			return list.get(0);
		else
			throw new RuntimeException("Empty molecule not found in database - critical error");
	}

	@Override
	@SuppressWarnings("unchecked")
	public MoleculeFormat getFormatByName(String _name) {
		List<MoleculeFormat> Format = 
				Globals.session().createQuery("from MoleculeFormat where name=:name")
				.setString("name", _name)
				.list();

		if(Format.size() > 0)
			return Format.get(0);
		else{
			MoleculeFormat format = new MoleculeFormat();
			format.name = _name;
			Globals.session().save(format);
			return format;
		}
	}

	@Override
	@SuppressWarnings("unchecked")
	public Mapping2 getMapping2(String inchi1, String inchi2, String data) 
	{
		String fullInchi2 = "";
		if (inchi2.length() == 32)
			fullInchi2 = inchi2;
		else
			fullInchi2 = inchi1+"-"+inchi2;

		Criteria criteria = Globals.session().createCriteria(Mapping2.class)
				.add(Restrictions.eq("inchi2", fullInchi2));

		List<Mapping2> mapping2s = criteria.list();

		if (mapping2s.size() == 0)
		{
			Mapping2 mapping2 = new Mapping2();
			mapping2.inchi2 = fullInchi2;


			qspr.fragmententities.Mapping2 externalMapping = qspr.fragmententities.Mapping2.get(inchi1, inchi2, data);			//Get Mapping2 ID from Central Repository...
			mapping2.id = externalMapping.id;                                //Theoretically garanties some level of consistency

			mapping2.mapping1 = Mapping1.get(inchi1, data);
			Globals.session().saveOrUpdate(mapping2);
			return mapping2;
		}
		else
			return mapping2s.get(0);
	}

	@Override
	public Mapping2 getMapping2(Integer mp2Id) {
		return (Mapping2) Globals.session().get(Mapping2.class, mp2Id);
	}

	@Override
	public Molecule getMolecule(long id)
	{
		return (Molecule) Globals.session().get(Molecule.class, id);
	}

	@Override
	public Molecule getBySolventName(String name) throws IOException {
	    Query q = Globals.session().createQuery("from ColoredName c where c.ep.article.id = :ar and c.name=:tn").setParameter("ar", OCHEMConfiguration.solvent).setParameter("tn", name);
	    @SuppressWarnings("unchecked")
		List<ColoredName> entities = q.list();

	    if(entities.size()>1)
	    	throw new IOException("For " + name + " two solvents were found.");

	    if(entities.size()==1)
	    	return entities.get(0).ep.molecule;
		return null;
	}
}
