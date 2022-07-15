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

package qspr;



// We still have to move many hard-coded variable names from the code here

/**
 * Enumeration for all variable names stored in HTTP sessions
 * See also Globals.setSessionVariable and Globals.getSessionVariable
 * @author midnighter
 * 
 */
// Midnighter on Sep 14, 2011
public enum SessionVariable 
{
	IMPORTED_TEMPLATE,
	MODEL_CONFIGURATOR,
	MODEL,
	MODEL_PROPERTY,
	MODEL_APPLIER,
	MODEL_APPLIER_INDEXMAP,
	BATCH_READER,
	MULTIPLE_MODELS_STARTER,
	MULTIPLE_MODELS_REPORT,
	BASKET_SELECT, // basket which is selected in the records browser
	ALERT_SCREENING_PROCESSOR,
	DESCRIPTORS_CALCULATOR_PROCESSOR,
	SETCOMPARE_PROCESSOR,
	EXPDESIGN_PROCESSOR,
	PROFILE,
	FILTERED_SELECTION_IDSET, 	// BatchEdit
	COMPRESSED_EP_MAP,			// BatchEdit,
	SHOPPING_CART, // oChemTrade shopping cart
	MOLOPTIMISER_PROCESSOR,
	MOLOPTIMISER_SETUP,
	MODEL_UPLOAD_RESULT,
	MMPOPTIMISER
}
