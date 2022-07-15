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

package qspr.fragmententities;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.FetchType;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.OneToMany;
import javax.xml.bind.annotation.XmlType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.Query;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;

@Entity
@XmlType(name="FragMapping1" )
@SuppressWarnings("unchecked")
public class Mapping1 
{
	private static transient final Logger logger = LogManager.getLogger(Mapping1.class);

	@Id
	@GeneratedValue
	@Column(name = "mapping1_id")
	@javax.persistence.OrderBy
	public Long id;

	@Column
	public String inchi1;

	/*
	 * = 0: Molecule has been fragmented
	 * > 0: Molecule is assigned as fragmentationTask
	 * 
	 * normal: 
	 * null Molecule will be sent as 1m Task
	 * - 1: Molecule will be sent as 1h Task
	 * - 2: Molecule will be sent as 24h Task
	 * - 3: Molecule timed out as 24h Task
	 * - 4: not used
	 * - 5: ??? restart as -1
	 * - 6: old -4 restart as -1
	 * - 7: old -3 restart as -1
	 * 
	 * errors:
	 * - 9: Molecule has no structure
	 * -10: Molecule has md5 as inchi
	 * -11: ECommerce Do Not Fragment
	 * 
	 * ( < 0: deprecated: Time in mins when task timed out)
	 * 
	 * -100001: some error during calculation
	 */
	@Column
	public Integer frag_status;

	//	@ManyToMany
	//	(
	//			targetEntity = Fragment.class,
	//			cascade={CascadeType.PERSIST, CascadeType.MERGE}
	//	)
	//	@JoinTable
	//	(
	//			name = "Mapping1Fragment",
	//			joinColumns = {@JoinColumn(name="mapping1_id")},
	//			inverseJoinColumns ={@JoinColumn(name="fragment_id")}
	//	)
	//@XmlTransient
	@OneToMany(mappedBy="mapping1", fetch=FetchType.LAZY)
	public Set<Mapping1Fragment> mapping1Fragments;


	@OneToMany(mappedBy="mapping1", fetch=FetchType.LAZY)
	public List<Mapping2> mapping2;

	public static Mapping1 getNonFragmentable(String inchi1, String molData)
	{
		Mapping1 result = null;
		Criteria criteria = Globals.alternateSession().createCriteria(Mapping1.class)
				.add(Restrictions.eq("inchi1", inchi1));

		List<Mapping1> mapping1s = criteria.list();

		if (mapping1s.size() == 0)
		{
			result = new Mapping1();
			result.inchi1 = inchi1;
			result.frag_status = -11;
			Globals.alternateSession().saveOrUpdate(result);	
		}
		else
		{
			result =  mapping1s.get(0);
		}
		return result;
	}


	public static Mapping1 get(String inchi1, String molData)
	{
		Mapping1 result = null;
		Criteria criteria = Globals.alternateSession().createCriteria(Mapping1.class)
				.add(Restrictions.eq("inchi1", inchi1));

		List<Mapping1> mapping1s = criteria.list();

		if (mapping1s.size() == 0)
		{
			result = new Mapping1();
			result.inchi1 = inchi1;
			result.frag_status = null;
			Globals.alternateSession().saveOrUpdate(result);
		} else
		{
			result =  mapping1s.get(0);

			if (Integer.valueOf(-11).equals(result.frag_status)) //ECommerce
			{
				result.frag_status = null;
				Globals.alternateSession().saveOrUpdate(result);
			}
		}
		return result;
	}

	public static List<Integer> getMappingIdListByFragmentInchi(String inchi1, List<Long> mapping1s)
	{

		if (mapping1s != null) 
		{ 
			if (mapping1s.size() > 11000)
			{
				List<Integer> result = new ArrayList<Integer>();
				java.util.Iterator<Long> iter = mapping1s.iterator();
				// Split by 10000 when making a query, otherwise we have an exponential time growth
				while (iter.hasNext())
				{
					List<Long> mp1ids = new ArrayList<Long>();
					int k = 0;
					while (iter.hasNext() && k++ < 10000)
						mp1ids.add(iter.next());
					result.addAll(getMappingIdListByFragmentInchi(inchi1, mp1ids));
				}
				return result;
			}
			else
			{
				logger.info("Sending " + mapping1s.size() + " to the fragments database");
				Query query = Globals.alternateSession()
						.createSQLQuery("select mf.mapping1_id from Fragment f natural join Mapping1Fragment mf where f.inchi1 = :inchi1 and mf.mapping1_id in (:mapping1s) order by mf.fragment_id desc, mf.mapping1_id desc")
						.addScalar("mapping1_id")
						.setParameterList("mapping1s", mapping1s)
						.setParameter("inchi1", inchi1);

				return query.list();
			}
		}
		else 
		{	
			List<Integer> result = new ArrayList<Integer>();
			Query query = Globals.alternateSession()
					.createSQLQuery("select mf.mapping1_id from Fragment f natural join Mapping1Fragment mf where f.inchi1 = :inchi1 order by mf.fragment_id desc, mf.mapping1_id desc")
					.addScalar("mapping1_id")
					//.setParameterList("mapping1s", mapping1s)
					.setParameter("inchi1", inchi1);

			result = query.list();

			return result;
		}
	}

	public void addFragment(String inchi1, String data, Integer inchi1_count)
	{
		Fragment fragment = Fragment.get(inchi1, data);
		//logger.info(this.id + " - " + fragment.id);
		Mapping1Fragment.get(this, fragment, inchi1_count);
	}


	public String getData() throws Exception 
	{
		if (mapping2.size() > 0)
		{
			Mapping2 map2 = mapping2.get(0);
			return map2.data;
		} else 
			throw new Exception("No structure data for " + inchi1);

		//		if (mapping1Fragment.size() == 1)
		//		{
		//			Fragment frag = mapping1Fragment.iterator().next();
		//			return frag.fragment_data;
		//		} else if (mapping1Fragment.size() > 1) 
		//		{
		//			throw new Exception("Molecule is already fragmented!");	
		//		} else
		//			throw new Exception("No structure data for " + inchi1);


	}

	public static List<Integer> getMapping1IdListByFragmentID(Integer id)
	{
		Query query = Globals.alternateSession()
				.createSQLQuery("select mf.mapping1_id from Mapping1Fragment mf where mf.fragment_id = :id order by mf.mapping1_id desc")
				.addScalar("mapping1_id")
				.setParameter("id", id);

		List<Integer> result = query.list();

		return result;
	}



}
