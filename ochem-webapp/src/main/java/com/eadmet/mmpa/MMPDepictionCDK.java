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
import java.util.Collection;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import javax.imageio.ImageIO;

import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.geometry.GeometryUtil;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.smsd.Isomorphism;
import org.openscience.cdk.smsd.interfaces.Algorithm;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import com.eadmet.exceptions.UserFriendlyException;
import com.google.common.collect.Iterables;
import com.google.common.collect.Sets;

import qspr.util.CDKUtils;

public class MMPDepictionCDK extends MMPDepiction {

		private IAtomContainer molecule;
		private IAtomContainer template;
		private Isomorphism isomorphism;
		private Set<IAtom> atomset = new HashSet<>();
		private Set<IBond> bondset = new HashSet<>();
		
		public static IAtomContainer prepare(IAtomContainer mol) throws CDKException {
			mol = AtomContainerManipulator.removeHydrogens(mol);
			
			if (!GeometryUtil.has2DCoordinates(mol)) {
				new StructureDiagramGenerator().generateCoordinates(mol);
			}
			
			return mol;
		}

		@Override
		public void initialize(String molecule, String template) {
			try {
//				IAtomContainer mol = prepare(CDKUtils.readOneMoleculeInAnyFormat(molecule));
//				IAtomContainer tmp = prepare(CDKUtils.readOneMoleculeInAnyFormat(template));
//				if (Iterables.size(mol.atoms()) >= Iterables.size(tmp.atoms())) {
//					this.molecule = mol;
//					this.template = tmp;
//				} else {
//					this.molecule = tmp;
//					this.template = mol;
//				}
				this.molecule = prepare(CDKUtils.readOneMoleculeInAnyFormat(molecule));
				this.template = prepare(CDKUtils.readOneMoleculeInAnyFormat(template));
			} catch (Exception e) {
				e.printStackTrace();
				throw new UserFriendlyException(e);
			}
			this.isomorphism = new Isomorphism(Algorithm.MCSPlus, false);
			try {
				isomorphism.init(this.template, this.molecule, true, true);
				isomorphism.setChemFilters(false, false, false);
				isomorphism.setBondInSensitiveTimeOut(0.2);
				isomorphism.setBondSensitiveTimeOut(0.2);
			} catch (CDKException e) {
				e.printStackTrace();
				throw new UserFriendlyException(e);
			}
		}

		@Override
		public void search() {
			template = isomorphism.getReactantMolecule();
			molecule = isomorphism.getProductMolecule();
			if (!(isomorphism.getAllMapping() == null)) {
				atomset = Sets.newHashSet(molecule.atoms());
				bondset = Sets.newHashSet(molecule.bonds());
				atomset.removeAll(isomorphism.getFirstAtomMapping().values());
				bondset.removeAll(isomorphism.getFirstBondMap().values());
			}
		}

		@Override
		public byte[] getDefaultImp() {
			DepictionGenerator dg = new DepictionGenerator().withSize(width, height)
                    .withAtomColors().withFillToFit();
			try {
				BufferedImage img = dg.withHighlight(atomset, mmpColor).withHighlight(bondset, mmpColor).depict(molecule).toImg();
				ByteArrayOutputStream baos = new ByteArrayOutputStream();
				ImageIO.write(img, format, baos);
				return baos.toByteArray();
			} catch (Exception e) {
				e.printStackTrace();
				throw new UserFriendlyException(e);
			}
		}
}
