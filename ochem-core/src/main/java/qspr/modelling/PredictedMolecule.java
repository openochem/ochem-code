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

package qspr.modelling;

import java.io.Serializable;

import qspr.entities.BasketEntry;
import qspr.entities.Molecule;

/**
 * A simple class for serialization of molecule references
 * @author midnighter
 *
 */
public class PredictedMolecule implements Serializable
{
	public Integer mp2ID;
	public Long molId;
	
	/**
	 * Link to the experimental property (optional)
	 */
	public Long epId;
	
	public PredictedMolecule(BasketEntry be)
	{
		epId = be.ep.id;
		if (be.ep.molecule == null)
		{
			be.ep.molecule = new Molecule();
			be.ep.molecule.id = -1L;
		}
		molId = be.ep.molecule.id;
		if (be.ep.molecule.mapping2 != null)
			mp2ID = be.ep.molecule.mapping2.id;
	}
	
	private static final long serialVersionUID = 1L;
}
