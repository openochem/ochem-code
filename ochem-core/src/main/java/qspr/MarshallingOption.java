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

/**
 * By default, JAXB marshalling invokes the annotated methods of java classes
 * This results into unnecessary calls and delays when this information is not really needed
 * This class contains the list of "marshalling options" used by these methods to prevent unnecessary invocations
 * See Globals.setMarshallingOption and Globals.getMarshallingOption
 * @author midnighter
 *
 */

// Midnighter on Sep 12, 2011

public enum MarshallingOption 
{
	 NO_BASKET_DETAILS,
	 MODEL_SET_SIZE,
	 ARTICLE_IN_MODEL_DOT,
	 ARTICLE_PENDING_TASKS,
	 PROPERTY_UNITCATEGORY,
	 PROPERTY_RECORD_COUNT,
	 PROPERTY_UNITS,
	 PROPERTY_OPTIONS,
	 PROPERTY_OPTIONS_FULL,
	 PROPERTY_OPTIONS_APPLIER,
	 UNIT_RECORDS_COUNT,
	 UNIT_PROPERTIES,
	 UNITCATEGORY_UNITS,
	 EXTENDED_BATCHUPLOAD_INFO,
	 USER_PUBLIC_MODELS,
	 PROPERTY_PREDICATES,
	 BASKET_CONDITIONS,
	 MODELTEMPLATE_FULLXML,
	 PROPERTY_MODERATOR,
	 BASKET_LIST_MODE,
	 EXPORTACTION_BONUSES,
	 USER_RECORDS,
	 PAIR_TRANSFORMATION,
	 TRANSFOMATION_ANNOTATIONS
}
