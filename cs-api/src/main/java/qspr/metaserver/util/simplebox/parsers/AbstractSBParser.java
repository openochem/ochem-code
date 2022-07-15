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

package qspr.metaserver.util.simplebox.parsers;

import java.util.ArrayList;
import java.util.List;
import java.util.UUID;

import org.json.simple.JSONObject;

import qspr.metaserver.util.simplebox.IMeta;
import qspr.metaserver.util.simplebox.ISBObject;
import qspr.metaserver.util.simplebox.MetaData;

public abstract class AbstractSBParser<T extends ISBObject> implements
		IJsonParser<T> {

	protected final JSONObject json;

	public AbstractSBParser(final JSONObject json) {
		this.json = json;
	}

	protected String parseProperty(String propName) {
		Object o = json.get(propName);
		if (o != null) {
			return o.toString();
		}
		return null;
	}

	protected List<String> parseList(String variableName) {
		Object o = json.get(variableName);
		if (o != null) {
			List<?> vals = List.class.cast(o);
			List<String> list = new ArrayList<String>(vals.size());
			for (Object element : vals) {
				list.add(element.toString());
			}
			return list;
		}
		return null;
	}

	protected void assertObjectTypePresent(String expectedType)
			throws JSONParsingException {
		String obj_type = parseProperty("object_type");
		if (obj_type == null) {
			throw new JSONParsingException("No obj_type field found");
		}
		if (!expectedType.equals(obj_type)) {
			throw new JSONParsingException("proper object_type not found. "
					+ "Object Type was : '" + obj_type + "'");
		}
	}
	
	/**
	 * Returns the metadata object that is attached to 
	 * the current object using the 'metadata' property
	 * @throws JSONParsingException 
	 */
	protected void attachMetadata(IMeta root) throws JSONParsingException{
		Object metadata = json.get("metadata");
		if (metadata != null) {
			JSONObject metaJson = (JSONObject) metadata;
			MetaDataParser metaParser = new MetaDataParser(metaJson);
			MetaData meta = metaParser.parse();
			root.setMeta(meta);

			if (meta.getIdentifiers() != null
					&& !meta.getIdentifiers().isEmpty()) {
				root.setUUID(UUID.fromString(meta.getIdentifiers().get(0)));
			}
		}
	}
}
