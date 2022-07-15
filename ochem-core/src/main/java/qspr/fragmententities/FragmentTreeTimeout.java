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

import java.io.Serializable;
import java.util.List;

import javax.persistence.Entity;
import javax.persistence.Id;

import org.hibernate.Criteria;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;

@Entity
public class FragmentTreeTimeout implements Serializable {

	/**
	 * create table FragmentTreeTimeout ( inchi1 varchar(32) primary key not null, seconds int not null) 
	 */
	private static final long serialVersionUID = 1342966780843953968L;

	@Id
	String inchi1;

	Long seconds;

	@SuppressWarnings("unchecked")
	public static FragmentTreeTimeout get(String inchi1) {

		Criteria c = Globals.alternateSession().createCriteria(FragmentTreeTimeout.class)
				.add(Restrictions.eq("inchi1", inchi1) );

		List<FragmentTreeTimeout> list = c.list();
		FragmentTreeTimeout ft;

		if (list.size() > 0)
			ft = list.get(0);
		else return null;	

		return ft;
	}

	@SuppressWarnings("unchecked")
	public static FragmentTreeTimeout set(String inchi1, Long seconds) {

		Criteria c = Globals.alternateSession().createCriteria(FragmentTreeTimeout.class)
				.add(Restrictions.eq("inchi1", inchi1) );

		List<FragmentTreeTimeout> list = c.list();
		FragmentTreeTimeout ft;
		if (list.size() > 0)
		{
			ft = list.get(0);
			ft.seconds = seconds;			
		} else {
			ft = new FragmentTreeTimeout();
			ft.inchi1 = inchi1;
			ft.seconds = seconds;
		}
		Globals.alternateSession().saveOrUpdate(ft);

		return ft;
	}

	public String toString() {
		return inchi1;
	}

}
