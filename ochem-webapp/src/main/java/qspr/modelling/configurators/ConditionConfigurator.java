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

package qspr.modelling.configurators;

import java.util.ArrayList;
import java.util.List;

import javax.servlet.http.HttpServletRequest;

import qspr.modelling.configurations.ExternalCondition;

/**
 * Helper Class to keep logic
 * @author itetko
 *
 */

public class ConditionConfigurator {

	static public List<ExternalCondition> addConditions(HttpServletRequest request) {

		// Conditions
		String[] conditions = request.getParameterValues("condition");
		if (conditions == null || conditions.length == 0) return null; // nothing is there !

		List<ExternalCondition> externalDescriptors =  new ArrayList<ExternalCondition>();

		for (String conditionId : conditions)
		{
			ExternalCondition eDesc = new ExternalCondition(Long.valueOf(conditionId));
			eDesc.defaultValue = Double.valueOf(request.getParameter((eDesc.getProperty().isQualitative() ? "option-" : "value-") + conditionId));
			if (eDesc.getProperty().isNumeric())
				eDesc.unitId = Long.valueOf(request.getParameter("unit-" + conditionId));
			externalDescriptors.add(eDesc);

			String[] mappingsLeft = request.getParameterValues("mapping-" + conditionId + "-1");
			String[] mappingsRight = request.getParameterValues("mapping-" + conditionId + "-2");

			if (eDesc.getProperty().isQualitative())
				if (mappingsLeft != null && mappingsLeft.length > 0)
					// User has specified some rules to merge options of a
					// qualitative condition
					for (int i = 0; i < mappingsLeft.length; i++)
						eDesc.addMerging(Long.valueOf(mappingsLeft[i]), Long.valueOf(mappingsRight[i]));
		}

		return externalDescriptors;

	}


}
