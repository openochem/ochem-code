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

package qspr.metaserver.configurations;

import javax.servlet.http.HttpServletRequest;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

@XmlRootElement(name = "inductivedescriptors-configuration")
public class DescriptorsInductiveDescriptorsConfiguration extends DescriptorsAbstractConfiguration
{
	private static final long serialVersionUID = -982127847571951884L;

	////////////////////////////
	//  Molecule descriptors  //
	///////////////////////////

	@XmlElement public boolean calcMElectronegativity = true;
	@XmlElement public boolean calcMHardness          = true;
	@XmlElement public boolean calcMSoftness          = true;
	@XmlElement public boolean calcMPartialCharges    = true;
	@XmlElement public boolean calcMInductiveParam    = true;
	@XmlElement public boolean calcMStericParam       = true;

	////////////////////////
	//  Atom descriptors  //
	////////////////////////

	@XmlElement public boolean calcAElectronegativity = false;
	@XmlElement public boolean calcAHardness          = false;
	@XmlElement public boolean calcASoftness          = false;
	@XmlElement public boolean calcAPartialCharges    = false;
	@XmlElement public boolean calcAInductiveParam    = false;
	@XmlElement public boolean calcAStericParam       = false;

	public String toString() { return ""; }

	// Just for debugging
	@XmlTransient
	public String getSummary()
	{
		StringBuffer buf = new StringBuffer();

		if (calcMElectronegativity) 
			buf.append("calcMElectronegativity\n");
		if (calcMHardness)         
			buf.append("calcMHardness\n");
		if (calcMSoftness)         
			buf.append("calcMSoftness\n");
		if (calcMPartialCharges)   
			buf.append("calcMPartialCharges\n");
		if (calcMInductiveParam)   
			buf.append("calcMInductiveParam\n");;
			if (calcMStericParam)
				buf.append("calcMStericParam\n");

			if (calcAElectronegativity) 
				buf.append("calcAElectronegativity\n");
			if (calcAHardness)         
				buf.append("calcAHardness\n");
			if (calcASoftness)         
				buf.append("calcASoftness\n");
			if (calcAPartialCharges)   
				buf.append("calcAPartialCharges\n");
			if (calcAInductiveParam)   
				buf.append("calcAInductiveParam\n");;
				if (calcAStericParam)
					buf.append("calcAStericParam\n");

				return buf.toString();
	}

	@Override
	public boolean requires3D() {
		return true;
	}

	@Override
	public String getDefaultTypeName() {
		return DescriptorsConfiguration.INDUCTIVE;
	}

	@Override
	DescriptorsAbstractConfiguration setConfiguration(HttpServletRequest request) {
		calcMElectronegativity = request.getParameter("inductiveMElectoneg") != null;
		calcMHardness = request.getParameter("inductiveMHardness") != null;
		calcMSoftness = request.getParameter("inductiveMSoftness") != null;
		calcMPartialCharges = request.getParameter("inductiveMCharge") != null;
		calcMInductiveParam = request.getParameter("inductiveMSigma") != null;
		calcMStericParam = request.getParameter("inductiveMRs") != null;
		return this;
	}
}
