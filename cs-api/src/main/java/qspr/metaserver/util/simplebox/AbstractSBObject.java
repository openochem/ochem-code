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
import java.io.NotSerializableException;
import java.io.StringWriter;
import java.util.UUID;

import org.json.simple.JSONAware;
import org.json.simple.JSONObject;

/**
 * 
 * 
 * 
 * @author Pantelis Sopasakis
 */

public abstract class AbstractSBObject implements ISBObject {

	private UUID uuid;

	private JSONObject json;

	public AbstractSBObject() {
		uuid = UUID.randomUUID();
	}

	public JSONAware getJson() {
		return json;
	}

	@Override
	public UUID getUuid() {
		return uuid;
	}

	@Override
	public void setUUID(UUID uuid) {
		this.uuid = uuid;
	}

	public String asJSONString() throws NotSerializableException {
		StringWriter out = new StringWriter();
		try {
			this.writeJSONString(out);
		} catch (final IOException ex) {
			throw new NotSerializableException(
					"The current object is not serializable to JSON");
		}
		return out.toString();

	}

}
