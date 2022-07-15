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

import java.io.IOException;

import javax.servlet.http.HttpServletRequest;
import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement(name="rdkitdescriptors-configuration")
public class DescriptorsRDKITConfiguration extends DescriptorsAbstractDragonConfiguration {

	private static final long serialVersionUID = 1L;
	public Double WHIM_THRESHOLD = 0.1;
	public Integer TOPOLOGICAL_NBITS = 1024;
	public Integer MORGAN_NBITS = 1024;
	public Integer MORGAN_RADIUS = 2;
	public Boolean MORGAN_FCFP = false;
	public Boolean MORGAN_COUNTS = false;
	public Integer AVALON_NBITS = 1024;
	public Boolean AVALON_COUNTS = null;

	// N.B.! block numbers are by by one larger, i.e. block for MORGAN is 11 and not 10

	final static public int  TOPOLOGICAL = 4,  WHIM = 8, MORGAN = 9, AVALON = 16;

	public DescriptorsRDKITConfiguration()
	{
		super(51199); // by default - selection of some and not all block!
	}

	public DescriptorsRDKITConfiguration(int dragonBlocks)
	{
		super(dragonBlocks);
	}

	@Override
	public int getBlocks() {
		return 18;
	}

	@Override
	public String getExtension() {
		return "rdkitb";
	}

	@Override
	public String toString()
	{
		String s=super.toString();
		switch(s){
		case "2D blocks: 10":
			return (MORGAN_FCFP != null && MORGAN_FCFP?"FCFP":"ECFP")+2*MORGAN_RADIUS + (MORGAN_COUNTS != null && MORGAN_COUNTS?"c":"");
		case "2D blocks: 17":
			return "AVALON"+ (AVALON_COUNTS != null && AVALON_COUNTS?"c":"");
		case "2D blocks: 18":
			return "Fingerprint";
		case "2D blocks: 11":
			return "MACCS";
		}

		return s;
	}

	@Override
	protected long get3DMask() {
		return 488;
	}

	@Override
	public String getDefaultTypeName() {
		return DescriptorsConfiguration.RDKIT;
	}

	@Override
	public DescriptorsAbstractConfiguration setAllOn() throws IOException{
		throw new IOException("setAllOn operation is not defined for " + toString());
	}

	@Override
	protected void setAllOn2D() throws IOException {
		throw new IOException("setAllOn operation is not defined for " + toString());
	}

	@Override
	public boolean isCachable() {
		return false;
	}

	@Override
	DescriptorsAbstractConfiguration setConfiguration(HttpServletRequest request) {

		if(super.setConfiguration(request) == null) return null;

		TOPOLOGICAL_NBITS = request.getParameter("topological_nbits")!=null? Integer.valueOf(request.getParameter("topological_nbits")) : null;
		WHIM_THRESHOLD = request.getParameter("whim_thresh")!=null? Double.valueOf(request.getParameter("whim_thresh")) : null;
		MORGAN_NBITS = request.getParameter("morgan_nbits")!=null? Integer.valueOf(request.getParameter("morgan_nbits")):null;
		MORGAN_RADIUS = request.getParameter("morgan_radius")!=null? Integer.valueOf(request.getParameter("morgan_radius")):null;
		MORGAN_FCFP = request.getParameter("morgan_fcfp")!=null ? true : null;
		MORGAN_COUNTS = request.getParameter("morgan_counts")!=null ? true : false;
		AVALON_NBITS = request.getParameter("avalon_nbits")!=null? Integer.valueOf(request.getParameter("avalon_nbits")):null;
		AVALON_COUNTS = request.getParameter("avalon_counts")!=null ? true : false;

		if(!bit(DescriptorsRDKITConfiguration.TOPOLOGICAL))TOPOLOGICAL_NBITS = null;
		if(!bit(DescriptorsRDKITConfiguration.WHIM))WHIM_THRESHOLD = null;
		if(!bit(DescriptorsRDKITConfiguration.MORGAN)) {
			MORGAN_NBITS = null;
			MORGAN_RADIUS = null;
			MORGAN_FCFP = null;
			MORGAN_COUNTS = null;
		}

		if(!bit(DescriptorsRDKITConfiguration.AVALON)) {
			AVALON_NBITS = null;
			AVALON_COUNTS = null;
		}

		if(MORGAN_FCFP  != null && !MORGAN_FCFP) MORGAN_FCFP = null;
		if(MORGAN_COUNTS  != null && !MORGAN_COUNTS) MORGAN_COUNTS = null;

		return this;
	}

}

