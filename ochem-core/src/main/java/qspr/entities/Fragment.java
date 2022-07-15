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

package qspr.entities;

import java.sql.Timestamp;
import java.util.Calendar;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.persistence.CascadeType;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.ManyToMany;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.Query;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import cern.colt.Timer;

@Entity
//@XmlRootElement(name = "fragment")
public class Fragment
{
	private static transient final Logger logger = LogManager.getLogger(Fragment.class);
	private static final int TIMEOUT_IN_MINUTES = 15;

	@Id
	@GeneratedValue
	@Column(name = "fragment_id")
	@XmlAttribute
	public Integer id;

	@Column
	@XmlElement
	public String inchi1;

	@Column(name = "lastused_time")
	public Timestamp time;

	@Column(name = "search_md5")
	public String searchHash;

	@ManyToMany(targetEntity = Mapping1.class, mappedBy = "moleculeFragment", cascade = { CascadeType.PERSIST,
		CascadeType.MERGE })
	@XmlElement
	public Set<Mapping1> mappings = new HashSet<Mapping1>();

	public static Fragment get(String inchi1)
	{
		return get(inchi1, null, null);
	}

	public static Fragment get(String inchi1, List<Long> mp1ids, String searchHash)
	{
		//////////////////////////////////////////////////////////////////
		// cleaning of temporary Fragment table and MoleculeFragment table
		Timestamp timeOut = new Timestamp(Calendar.getInstance().getTimeInMillis() - TIMEOUT_IN_MINUTES * 60 * 1000);

		Criteria clean = Globals.session().createCriteria(Fragment.class).add(Restrictions.lt("time", timeOut));
		@SuppressWarnings("unchecked")
		List<Fragment> fragmentsToDelete = clean.list();
		for (Fragment fragment : fragmentsToDelete)
		{
			Integer fragment_id = fragment.id;
			Globals.session().createSQLQuery("delete from Mapping1Fragment where fragment_id = :fragment_id")
			.setInteger("fragment_id", fragment_id).executeUpdate();
			Globals.session().delete(fragment);
			Globals.session().flush();
		}

		/////////////////////////////////////////
		// getting or filling of temporary tables
		logger.info("searchHash: " + searchHash);
		Criteria c = Globals.session().createCriteria(Fragment.class).add(Restrictions.eq("inchi1", inchi1)).add(
				Restrictions.eq("searchHash", searchHash));
		@SuppressWarnings("unchecked")
		List<Fragment> list = c.list();

		if (list.size() > 0)
		{
			Fragment frag = list.get(0);
			frag.time = new Timestamp(Calendar.getInstance().getTimeInMillis());
			Globals.session().saveOrUpdate(frag);
			return frag;
		}
		else
		{
			Fragment frag = new Fragment();
			frag.inchi1 = inchi1;
			frag.searchHash = searchHash;
			frag.time = new Timestamp(Calendar.getInstance().getTimeInMillis());
			Globals.session().saveOrUpdate(frag);
			Timer t = new Timer();
			t.start();
			logger.info("Getting the list of molecules containing a fragment");
			List<Integer> fragList = qspr.fragmententities.Mapping1.getMappingIdListByFragmentInchi(inchi1, mp1ids);
			logger.info("Got the list of " + fragList.size() + " molecules in " + t.seconds() + "s");
			t.stop();
			t = new Timer();
			t.start();
			logger.info("Inserting molecules to a temporary table");
			while (fragList.size() > 0)
			{
				int count = 0;
				StringBuffer sb = new StringBuffer();
				while (fragList.size() > 0 && count < 5000)
				{
					sb.append(",(");
					sb.append(frag.id);
					sb.append(",");
					sb.append(fragList.remove(0));
					sb.append(")");
					count++;
				}
				Query query = Globals.session().createSQLQuery(
						"insert ignore into Mapping1Fragment(fragment_id, mapping1_id) values "
								+ sb.deleteCharAt(0).toString());
				query.executeUpdate();
				System.out.print("..." + t.seconds());
			}
			logger.info("\n");
			logger.info("Fragment insert performed in " + t.seconds() + "s");
			Globals.session().flush();
			t.stop();
			return frag;
		}
	}


}
