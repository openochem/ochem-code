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
import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement(name="spectrophores-configuration")
public class DescriptorsSpectrophoresConfiguration extends DescriptorsAbstractConfiguration {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public static final int cageTypeNo=0;
	public static final int cageTypeUnique=1;
	public static final int cageTypeMirror=2;
	public static final int cageTypeAll=3;

	int accuracy=20;
	int cage=cageTypeNo;
	float resolution=0.1f;

	public void setAccuracy(int accuracy){
		this.accuracy=accuracy;
	}

	public int getAccuracy(){
		return accuracy;
	}

	public void setCageType(int cageType){
		this.cage=cageType;
	}

	public int getCageType(){
		return cage;
	}

	public String getCageTypeString(){
		return cage==cageTypeNo?"No":
			cage==cageTypeUnique?"Unique":
				cage==cageTypeMirror?"Mirror":
					cage==cageTypeAll?"All":"Error"
						;
	}

	public void setResolution(float resolution){
		this.resolution=resolution;
	}

	public float getResolution(){
		return resolution;
	}

	public String toString()
	{
		return "accuracy="+accuracy;
	}

	@Override
	public boolean requires3D() {
		return true;
	}

	@Override
	public String getDefaultTypeName() {
		return DescriptorsConfiguration.Spectrophores;
	}

	@Override
	DescriptorsSpectrophoresConfiguration setConfiguration(HttpServletRequest request) {
		setAccuracy(Integer.valueOf(request.getParameter("spectrophores-accuracy")));
		setCageType(Integer.valueOf(request.getParameter("spectrophores-cage")));
		setResolution(new Float(request.getParameter("spectrophores-resolution")));
		return this;
	}	
}
