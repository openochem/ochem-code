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

import org.hibernate.Criteria;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.annotations.Loggable;

@Entity
@Loggable
public class OAuthCode {
	
	@Id
	@Column(name = "code_id")
	@GeneratedValue
	private Long id;

	@ManyToOne
	@JoinColumn(name = "app_id")
	public OAuthApp app;
	
	@Column(name= "code", unique = true)
	private String code;
	
	@Column(name= "redirect_uri")
	public String redirectURI;
	
	@ManyToOne
	@JoinColumn(name = "user_id")
	public User user;
	
	public static OAuthCode create(OAuthApp app, String redirect_uri, User user) {
		OAuthCode code = new OAuthCode();
		code.user = user;
		code.app = app;
		code.redirectURI = redirect_uri;
		code.code = User.generatePassword(32).replaceAll("[^a-zA-Z\\d\\s:]", "");
		return code;
	}
	
	public void save() {
		Globals.session().save(this);
	}
	
	@Override
	public String toString() {
		return this.code;
	}
	
	public static OAuthCode get(String code) {
		Criteria criteria = Globals.session().createCriteria(OAuthCode.class).add(Restrictions.eq("code", code));
		List ret = criteria.list();
		if (ret.isEmpty()) {
			return null;
		}
		return (OAuthCode) ret.get(0);
	}
	
	public void delete() {
		Globals.session().delete(this);
	}
}
