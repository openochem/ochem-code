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

package qspr.depiction;

import java.io.IOException;

import qspr.dao.ChemInfEngine;

abstract public class Depiction {

	protected ChemInfEngine engine = null;

	public abstract double getHeight();

	public abstract void setHeight(double height);

	public abstract double getWidth();

	public abstract void setWidth(double width);

	public abstract String getFormat();

	public abstract void setFormat(String format);

	public byte[] getImage() throws IOException{

		try {
			byte [] img = getDefaultImp();
			if(img != null) return img;
			return getSecondImp(this);
		}catch(Exception e) {
			return getSecondImp(this);
		}
	}

	public abstract byte[] getSecondImp(Depiction e) throws IOException;

	public abstract byte[] getDefaultImp() throws IOException;

	public abstract void setDims(double width, double height);
}
