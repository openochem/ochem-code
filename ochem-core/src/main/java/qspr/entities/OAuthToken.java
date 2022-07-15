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

import java.util.List;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

import org.hibernate.Criteria;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.annotations.Loggable;

@Entity
@XmlRootElement(name = "token")
@Loggable
public class OAuthToken {
	
	@Id
	@Column(name = "token_id")
	@GeneratedValue
	private Long id;

	@ManyToOne
	@JoinColumn(name = "app_id")
	public OAuthApp app;
	
	@Column(name= "token", unique = true)
	@XmlElement(name = "value")
	private String token;
	
	@ManyToOne
	@JoinColumn(name = "user_id")
	public User user;
	
	public static OAuthToken create(OAuthApp app, User user) {
		OAuthToken token = new OAuthToken();
		token.app = app;
		token.token = User.generatePassword(64).replaceAll("[^a-zA-Z\\d\\s:]", "");
		token.user = user;
		return token;
	}
	
	public void save() {
		Globals.session().save(this);
	}
	
	public void delete() {
		Globals.session().delete(this);
	}
	
	@Override
	public String toString() {
		return this.token;
	}
	
	public static OAuthToken get(String token) {
		Criteria criteria = Globals.session().createCriteria(OAuthToken.class).add(Restrictions.eq("token", token));
		List ret = criteria.list();
		if (ret.isEmpty()) {
			return null;
		}
		return (OAuthToken) ret.get(0);
	}
}
