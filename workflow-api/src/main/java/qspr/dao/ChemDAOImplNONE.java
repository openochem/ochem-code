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
import java.io.InputStream;
import java.util.List;
import java.util.concurrent.TimeoutException;

import org.apache.commons.lang.NotImplementedException;

import qspr.workflow.datatypes.DataTable;

public class ChemDAOImplNONE extends ChemDAO{

	@Override
	public ChemInfEngine engine() {
		return ChemInfEngine.NONE;
	}

	@Override
	public String aromatize(String molecule, Aromatisation type) throws IOException {
		throw new NotImplementedException("is not supposed to be called with ChemInfEngine.NONE");
	}

	@Override
	public String convertToSmilesOrSmart(String molecule, String format) throws IOException {
		throw new NotImplementedException("is not supposed to be called with ChemInfEngine.NONE");
	}

	@Override
	public String convertToSDFUnblockedImp(String molecule) throws IOException, TimeoutException {
		throw new NotImplementedException("is not supposed to be called with ChemInfEngine.NONE");
	}

	@Override
	public List<String> readSDFMolsFromFile(String filePath) throws IOException {
		throw new NotImplementedException("is not supposed to be called with ChemInfEngine.NONE");
	}

	@Override
	public String convertToFormat(String molecule, String format) throws IOException {
		throw new NotImplementedException("is not supposed to be called with ChemInfEngine.NONE");
	}

	@Override
	public String guessFormat(String molecule) throws IOException {
		throw new NotImplementedException("is not supposed to be called with ChemInfEngine.NONE");
	}

	@Override
	public double getMassImp(String sdfData) {
		throw new NotImplementedException("is not supposed to be called with ChemInfEngine.NONE");
	}

	@Override
	public DataTable getAllDataInSDF(InputStream inp) throws IOException {
		throw new NotImplementedException("is not supposed to be called with ChemInfEngine.NONE");
	}

	@Override
	protected String getMaxComponent(String molecule) throws IOException {
		throw new NotImplementedException("is not supposed to be called with ChemInfEngine.NONE");
	}

	@Override
	public String getFormulaImpl(String molecule) throws IOException {
		throw new NotImplementedException("is not supposed to be called with ChemInfEngine.NONE");
	}

	@Override
	public String convertToCanonicalName(String mol) {
		throw new NotImplementedException("is not supposed to be called with ChemInfEngine.NONE");
	}

	@Override
	public String convertToKekuleSMILES(String mol) throws IOException {
		throw new NotImplementedException("is not supposed to be called with ChemInfEngine.NONE");
	}

	@Override
	public String getAtomProperty(Properties property, String sdf) throws IOException {
		throw new NotImplementedException("is not supposed to be called with ChemInfEngine.NONE");
	}

	@Override
	public int getAtomCount(String sdf) throws IOException {
		throw new NotImplementedException("is not supposed to be called with ChemInfEngine.NONE");
	}

	@Override
	public int getBreakableBoundCount(String sdf) throws Exception {
		throw new NotImplementedException("is not supposed to be called with ChemInfEngine.NONE");
	}

	@Override
	public String[] getHydrogenFragmentsReplacedWithAl(String molecule) throws Exception {
		throw new NotImplementedException("is not supposed to be called with ChemInfEngine.NONE");
	}

	@Override
	public String[] splitOrderedByChargeAndSize(String data) throws IOException {
		throw new NotImplementedException("is not supposed to be called with ChemInfEngine.NONE");
	}

	@Override
	public boolean isFullyMappedSMIRKS(String reaction) {
		throw new NotImplementedException("is not supposed to be called with ChemInfEngine.NONE");
	}

	@Override
	public List<String> getInChIComponents(String mol, ION type) {
		throw new NotImplementedException("is not supposed to be called with ChemInfEngine.NONE");
	}

	@Override
	public double getTanimoto(byte[] fp_a, byte[] fp_b) throws Exception {
		throw new NotImplementedException("is not supposed to be called with ChemInfEngine.NONE");
	}

	@Override
	public double getTanimoto(String sdf_a, byte[] fp_b) throws Exception {
		throw new NotImplementedException("is not supposed to be called with ChemInfEngine.NONE");
	}

	@Override
	public byte[] getFingerprint(String sdf_a) throws Exception {
		throw new NotImplementedException("is not supposed to be called with ChemInfEngine.NONE");
	}

	@Override
	public String addPropertyToSDF(String sdf, String propertyName, String propertyValue) throws IOException {
		throw new NotImplementedException("is not supposed to be called with ChemInfEngine.NONE");
	}

}
