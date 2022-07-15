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

package qspr.toxicity;

import org.springframework.web.servlet.mvc.support.ControllerClassNameHandlerMapping;

import qspr.controllers.AlertsController;

// The class takes into account some customized URLs
public class CustomOCHEMClassNameHandlerMapping extends ControllerClassNameHandlerMapping
{
	@SuppressWarnings("rawtypes")
	protected String[] generatePathMappings(Class beanClass) 
	{
		if (beanClass.equals(AlertsController.class))
			return new String[]{"/substructurealert", "/substructurealert/*", "/alerts", "/alerts/*"};
		return super.generatePathMappings(beanClass);
	}
}
