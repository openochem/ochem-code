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

package qspr.dao;

import java.io.IOException;

import qspr.entities.Mapping2;
import qspr.entities.Molecule;
import qspr.entities.MoleculeFormat;
import qspr.entities.OriginalMolecule;

public interface MoleculeDAO {
	public OriginalMolecule getOriginalMoleculeByStructure(String md5);
	public Molecule getMoleculeByMD5(String md5);
	public Molecule getMoleculeByInvariants(String picMd5, String fullInchi2);
	public Molecule getEmptyMolecule();
	public Molecule getMolecule(long id);
	public MoleculeFormat getFormatByName(String name);
	public Mapping2 getMapping2(String inchi1, String inchi2, String data);
	public Mapping2 getMapping2(Integer mp2Id);
	public Molecule getBySolventName(String name) throws IOException;
}
