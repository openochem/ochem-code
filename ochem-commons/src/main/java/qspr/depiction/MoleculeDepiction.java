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

import java.awt.Color;
import java.io.IOException;

import javax.servlet.http.HttpServletRequest;

import com.eadmet.exceptions.UserFriendlyException;

import qspr.dao.ChemDAO;
import qspr.dao.ChemInfEngine;

public abstract class MoleculeDepiction extends Depiction {
	// image props
	protected double height = 400;
	protected double width = 400;
	protected String format = "png"; 
	protected boolean hideHydrogens = true; 
	protected boolean aromatize = false; 
	protected boolean isError = false; 
	protected Color color = new Color(255, 255, 255);
	protected Color errcolor = new Color(245, 169, 169);
	protected double alpha = 1.0;
	protected String alphaPrefix = "FF";
	protected String additionalStringConfig = "";

	// choose implementation
	//	protected final static ChemInfEngine defaultImpl = OCHEMConfiguration.getCheminfEngine();
	// FIXME: use this to determine the engine to use

	// molecule data
	String moleculeData;

	protected MoleculeDepiction() {
		// no action, hidden
	}

	public MoleculeDepiction(MoleculeDepiction config) {
		super();
		configure(config);
	}

	public static MoleculeDepiction get(ChemInfEngine choice) {

		MoleculeDepiction imp = null;

		try {
			switch (choice) {
			case CDK:
				imp = (MoleculeDepiction) Class.forName("qspr.depiction.MoleculeDepictionCDK").newInstance();
				break;    
			case CHEMAXON:
				imp = (MoleculeDepiction) Class.forName("qspr.depiction.MoleculeDepictionChemAxon").newInstance();
				break;
			default:
				throw new UserFriendlyException("Molecule Depiction implementation unavailable: " + choice);
			}
		} catch (InstantiationException | IllegalAccessException | ClassNotFoundException e) {
			throw new UserFriendlyException(e);
		}

		imp.engine = choice;
		return imp;
	};

	public static MoleculeDepiction get(MoleculeDepiction config, ChemInfEngine choice) {
		MoleculeDepiction ret = get(choice);
		ret.configure(config);
		return ret;
	}

	public byte[] getSecondImp(Depiction e) throws IOException{
		for(ChemInfEngine eng: ChemInfEngine.values())
			if(!ChemDAO.ignoreEngine(eng)) {
				MoleculeDepiction imp = get(eng);
				if(imp == null)return null;
				imp.configure(this);
				imp.moleculeData = this.moleculeData;
				return imp.getDefaultImp();
			}
		return null;			
	}


	//	public static MoleculeDepiction get() {
	//		return get(defaultImpl);
	//	}

	public void configure(MoleculeDepiction config) {
		this.height = config.getHeight();
		this.width = config.getWidth();
		this.format = config.getFormat();
		this.hideHydrogens = config.isHideHydrogens();
		this.aromatize = config.isAromatize();
		this.color = config.getColor();
		this.errcolor = config.getErrorColor();
		this.alpha = config.getAlpha();
		this.additionalStringConfig = config.additionalStringConfig;
	}

	public void setMolecule(String mol) throws IOException {
		this.moleculeData = mol;
	}

	//	public void setMolecule(Molecule mol) throws IOException {
	//		// TODO: we could use a mol specified like this (with a call to mol.getImage()) to save a default image with the default settings to the db to no longer recalculate it
	////		if (!customFormat)
	//// 	    	{
	//// 	    		imgBinary = dbMolecule.getImage(); // see this method of Molecule to see the default format (should be the same as the default attribute values of this class)
	//// 	    	}
	//		if(mol.mapping1.inchi1.length() != 14) {
	//    		this.setError(true);
	//    	}
	//		this.moleculeData = mol.getData();
	//	}

	public String getMolecule() throws IOException {
		return moleculeData;
	}

	public Color getErrorColor() {
		return errcolor;
	}

	public void setErrorColor(Color errcolor) {
		this.errcolor = errcolor;
	}

	public void setErrorColor(String errcolor) {
		this.errcolor = Color.decode("#" + errcolor);
	}

	public Color getColor() {
		return color;
	}

	public void setColor(Color color) {
		this.color = color;
	}

	public void setColor(String color) {
		this.color = Color.decode("#" + color);
	}

	public double getAlpha() {
		return alpha;
	}

	public void setAlpha(double alpha) {
		this.alpha = alpha;
		initAlpha();
	}

	private String colorToHex(Color color) {
		return Integer.toHexString(color.getRed()) + Integer.toHexString(color.getGreen()) + Integer.toHexString(color.getBlue());
	}

	private void initAlpha() {
		alphaPrefix = Integer.toHexString((int)Math.round(255 * alpha));
		if (alphaPrefix.length() < 2)
			alphaPrefix = "0" + alphaPrefix;
	}

	public boolean isHideHydrogens() {
		return hideHydrogens;
	}

	public void setHideHydrogens(boolean hideHydrogens) {
		this.hideHydrogens = hideHydrogens;
	}

	@Override
	public double getHeight() {
		return height;
	}

	@Override
	public void setHeight(double height) {
		this.height = height;
	}

	@Override
	public double getWidth() {
		return width;
	}

	@Override
	public void setWidth(double width) {
		this.width = width;
	}

	@Override
	public String getFormat() {
		return format;
	}

	@Override
	public void setFormat(String format) {
		this.format = format;
	}

	public boolean isAromatize() {
		return aromatize;
	}

	public void setAromatize(boolean aromatize) {
		this.aromatize = aromatize;
	}

	public void setDims(double width, double height) {
		setWidth(width);
		setHeight(height);
	}

	public boolean isError() {
		return isError;
	}

	public void setError(boolean isError) {
		this.isError = isError;
	}

	public String toStringConfig() {
		String format = getFormat() + ":w" + getWidth() + ",h" + getHeight();
		if (hideHydrogens)
			format += ",H_off";
		if (isError) {
			format += ",#" + colorToHex(errcolor);
		} else {
			format += ",#" + alphaPrefix.toLowerCase() + colorToHex(color);
		}
		if (!additionalStringConfig.equals("")) {
			format += additionalStringConfig;
		}
		return format;
	}

	public void appendToStringConfig(String toAppend) {
		this.additionalStringConfig += toAppend;
	}

	public void configure(HttpServletRequest request) {
		// Retrieving GET/POST parameters
		if (request.getParameter("alpha") != null)
		{
			this.setAlpha(Double.valueOf(request.getParameter("alpha")));
		}

		String type = (request.getParameter("type") != null) ? request.getParameter("type") : "png";
		this.setFormat(type);

		String w = (request.getParameter("w") != null) ? request.getParameter("w") : "150";
		String h = (request.getParameter("h") != null) ? request.getParameter("h") : "150";
		this.setDims(Double.parseDouble(w), Double.parseDouble(h));

		String validColor = (request.getParameter("color") != null) ?  request.getParameter("color") : "FFFFFF";
		this.setColor(validColor);

		String errorColor = (request.getParameter("ecolor") != null) ?  request.getParameter("ecolor") : "F5A9A9";
		this.setErrorColor(errorColor);

		boolean hideHydrogens = true;
		if (request.getParameter("hh") != null && request.getParameter("hh").equals("false")) {
			// only turn on hydrogens if hiding is explicitly off
			hideHydrogens = false;
		}
		this.setHideHydrogens(hideHydrogens);

		boolean aromatize = request.getParameter("aromatize") != null;
		this.setAromatize(aromatize);
	}
}
