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

package com.eadmet.mmpa;

import java.awt.image.BufferedImage;
import java.io.ByteArrayOutputStream;
import java.io.IOException;

import javax.imageio.ImageIO;

import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import com.eadmet.exceptions.UserFriendlyException;

import qspr.util.CDKUtils;

public class MMPFragDepictionCDK extends MMPFragDepiction {

	IAtomContainer ac;

	@Override
	public byte[] getDefaultImp() throws IOException {
		DepictionGenerator dg = new DepictionGenerator().withSize(width, height)
                .withAtomColors().withFillToFit();
		try {
			BufferedImage img = dg.depict(ac).toImg();
			ByteArrayOutputStream baos = new ByteArrayOutputStream();
			ImageIO.write(img, format, baos);
			return baos.toByteArray();
		} catch (Exception e) {
			e.printStackTrace();
			throw new UserFriendlyException(e);
		}
	}

	@Override
	public void setMolecule(String mol) throws IOException {
		ac = CDKUtils.readOneMoleculeInAnyFormat(mol.replaceAll("\\[Al\\]", "*"));
		if (hideHydrogens) {
			AtomContainerManipulator.suppressHydrogens(ac);
		}
	}

}
