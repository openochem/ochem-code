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

import java.net.MalformedURLException;
import java.net.URL;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.UUID;  

import javax.crypto.IllegalBlockSizeException;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.hibernate.Criteria;
import org.hibernate.criterion.Restrictions;

import com.eadmet.utils.CryptUtils;

import qspr.Globals;
import qspr.annotations.Loggable;

@Entity
@XmlRootElement(name = "app")
@Loggable
public class OAuthApp {

	@Id
	@Column(name = "app_id")
	@GeneratedValue
	@XmlAttribute
	public Long id;
	
	@ManyToOne
	@JoinColumn(name = "user_id")
	@XmlTransient
	public User user;
	
	@Column(name = "client_id", unique = true)
	public String clientID;
	
	@Column(name= "domain_list", unique = false)
	public String domainList;
	
	@Column(name= "client_secret", unique = true)
	@XmlTransient
	private String encodedSecret;
	
	@XmlElement(name = "secret")
	public String getSecret() throws Exception
	{
		try
		{
			return CryptUtils.desDecode(encodedSecret);
		}
		catch (IllegalBlockSizeException e)
		{
			return encodedSecret;
		}
	}
	
	public void generateSecret() throws Exception 
	{
		String generated_secret = User.generatePassword(32).replaceAll("[^a-zA-Z\\d\\s:]", "");
		this.encodedSecret = CryptUtils.desEncode(generated_secret);
	}
	
	public void generateClientID() {
		this.clientID = UUID.randomUUID().toString();
	}
	
	public void init() throws Exception {
		generateClientID();
		generateSecret();
	}
	
	public void save() {
		Globals.session().save(this);
	}
	
	public void delete() {
		Globals.session().delete(this);
	}
	
	public static List<OAuthApp> getAppsForUser(User user) {
		Criteria criteria = Globals.session().createCriteria(OAuthApp.class).add(Restrictions.eq("user", user));
		return criteria.list();
	}

	public static OAuthApp create(User user) throws Exception {
		OAuthApp app = new OAuthApp();
		app.user = user;
		app.init();
		return app;
	}
	
	public static OAuthApp get(String id) {
		Criteria criteria = Globals.session().createCriteria(OAuthApp.class).add(Restrictions.eq("clientID", id));
		List ret = criteria.list();
		if (ret.isEmpty()) {
			return null;
		}
		return (OAuthApp) ret.get(0);
	}
	
	public Set<String> getDomainSet() {
		Set<String> ret = new HashSet<>();
		ret.addAll(Arrays.asList(domainList.split(",")));
		return ret;
	}

	public boolean checkDomain(String redirect_uri) throws MalformedURLException {
		URL url = new URL(redirect_uri);
		String host = url.getHost();
		return getDomainSet().contains(host);
	}
}
