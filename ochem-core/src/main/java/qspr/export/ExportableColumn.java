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

package qspr.export;

import javax.xml.bind.annotation.XmlEnum;

/**
 * Enumerates item types exportable from OCHEM
 * 1 item can correspond to several columns in the exported file
 * 
 * @author midnighter
 *
 */

@XmlEnum
public enum ExportableColumn
{
	SMILES, 
	CASRN, 
	RECORDID, 
	MOLECULEID, 
	EXTERNAL_ID,  // was renames, since it is not clear in "which house" it is "inhouse". 
	N,
	NAMES, 
	INTRODUCER,
	MODIFIER,
	ARTICLE,
	ERROR,
	EXP_VALUE,
	PREDICTED_VALUE, 
	EXP_VALUE_CONVERTED, 
	DM_VALUE, 
	ACCURACY,
	APPLICABILITY_DOMAIN,
	CONDITIONS, 
	DESCRIPTORS,
	DESCRIPTORSNAMES,
	COMMENTS,
	INCHI_KEY,
	COMPRESSED
	;

	public String toString()
	{
		switch (this)
		{
		case EXP_VALUE:
			return "Experimentaly measured values";
		case EXP_VALUE_CONVERTED:
			return "Experimentaly measured values (in converted units)";				
		case PREDICTED_VALUE:
			return "Predicted values (in converted units)";
		case DM_VALUE:
			return "DM (distance to model) values";
		case CONDITIONS:
			return "Conditions of experiments";
		case INTRODUCER:
			return "Introducers of the records";
		case MODIFIER:
			return "Last modifiers of the records";
		case ARTICLE:
			return "Publication IDs";
		case ERROR:
			return "Error messages";
		case EXTERNAL_ID:
			return "External unique identifier";
		case INCHI_KEY:
			return "Inchi-key";
		case COMMENTS:
			return "Comments";
		case SMILES:
			return "Structure (SMILES or SDF)";
		case N:
			return "Identifier in article (N)";
		case APPLICABILITY_DOMAIN:
			return "Applicability Domain (FALSE if predictions are outside of the AD)";
		case COMPRESSED:
			return "Merge information for the same molecule";
		default:
			return this.name();
		}
	}

	//	final static Set<ExportableColumn> DEFAULT_EXPORT_SET = 
	//			new HashSet<ExportableColumn>(Arrays.asList(new ExportableColumn[]{
	//					SMILES, CASRN, NAMES, N, EXP_VALUE, PREDICTED_VALUE}));;

	public boolean isRestricted()
	{
		//List<ExportableColumn> restrictedAccessColumns = Arrays.asList(new ExportableColumn[]{SMILES, ARTICLE, CONDITIONS, NAMES, CASRN, INCHI_KEY});
		//return restrictedAccessColumns.contains(this);
		return false;
	}

}
