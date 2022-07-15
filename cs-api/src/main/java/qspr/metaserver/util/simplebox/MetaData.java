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

package qspr.metaserver.util.simplebox;

import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;

import org.json.simple.JSONValue;

/**
 * 
 * 
 * 
 * @author pantelis
 */

public class MetaData extends AbstractSBObject {

	private List<String> description;
	private List<String> identifiers;
	private List<String> comments;
	private List<String> rights;
	private List<String> seeAlso;
	private List<String> title;
	private List<String> subject;

	public MetaData() {
	}

	// <editor-fold defaultstate="collapsed" desc="Getters and Setters">
	public List<String> getComments() {
		return comments;
	}

	public void setComments(List<String> comments) {
		this.comments = comments;
	}

	public List<String> getDescription() {
		return description;
	}

	public void setDescription(List<String> description) {
		this.description = description;
	}

	public List<String> getIdentifiers() {
		return identifiers;
	}

	public void setIdentifiers(List<String> identifiers) {
		this.identifiers = identifiers;
	}

	public List<String> getRights() {
		return rights;
	}

	public void setRights(List<String> rights) {
		this.rights = rights;
	}

	public List<String> getSeeAlso() {
		return seeAlso;
	}

	public void setSeeAlso(List<String> seeAlso) {
		this.seeAlso = seeAlso;
	}

	public List<String> getSubject() {
		return subject;
	}

	public void setSubject(List<String> subject) {
		this.subject = subject;
	}

	public List<String> getTitle() {
		return title;
	}

	public void setTitle(List<String> title) {		this.title = title;

	}
	// </editor-fold>


	// <editor-fold defaultstate="collapsed" desc="Adders">
	public MetaData addDescriptions(String... desc) {
		if (description == null) {
			description = new ArrayList<String>();
		}
		description.addAll(Arrays.asList(desc));
		return this;
	}

	public MetaData addIdentifiers(String... desc) {
		if (identifiers == null) {
			identifiers = new ArrayList<String>();
		}
		identifiers.addAll(Arrays.asList(desc));
		return this;

	}

	public MetaData addComments(String... desc) {
		if (comments == null) {
			comments = new ArrayList<String>();
		}
		comments.addAll(Arrays.asList(desc));
		return this;

	}

	public MetaData addRights(String... desc) {
		if (rights == null) {
			rights = new ArrayList<String>();
		}
		rights.addAll(Arrays.asList(desc));
		return this;
	}

	public MetaData addSeeAlso(String... see) {
		if (seeAlso == null) {
			seeAlso = new ArrayList<String>();
		}
		seeAlso.addAll(Arrays.asList(see));
		return this;

	}

	public MetaData addTitles(String... tit) {
		if (title == null) {
			title = new ArrayList<String>();
		}
		title.addAll(Arrays.asList(tit));
		return this;
	}

	public MetaData addSubjects(String... desc) {
		if (subject == null) {
			subject = new ArrayList<String>();
		}
		subject.addAll(Arrays.asList(desc));
		return this;
	}
	// </editor-fold>

	public void writeJSONString(Writer writer) throws IOException {
		HashMap<String,Object> obj = new LinkedHashMap<String,Object>();
		obj.put("object_type", "MetaData");
		if (title != null && !title.isEmpty())
			obj.put("title", title);

		if (subject != null && !subject.isEmpty())
			obj.put("subject", subject);

		if (description != null && !description.isEmpty())
			obj.put("description", description);

		if (identifiers != null && !identifiers.isEmpty())
			obj.put("identifier", identifiers);
		if (comments != null && !comments.isEmpty())
			obj.put("comment", comments);
		if (rights != null && !rights.isEmpty())
			obj.put("rights", rights);
		if (seeAlso != null && !seeAlso.isEmpty())
			obj.put("seeAlso", seeAlso);
		JSONValue.writeJSONString(obj, writer);
	}

}