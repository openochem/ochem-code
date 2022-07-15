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

import org.springframework.web.servlet.mvc.multiaction.AbstractUrlMethodNameResolver;

public class SophisticatedMethodNameResolver extends AbstractUrlMethodNameResolver 
{
	private boolean allToOne = true;
	
	public boolean isAllToOne()
	{
		return allToOne;
	}

	public void setAllToOne(boolean allToOne)
	{
		this.allToOne = allToOne;
	}

	protected String getHandlerMethodNameForUrlPath(String urlPath) {
		
		if (allToOne)
			return "show";
		else
		{
			String parts[] = urlPath.split("/");
			if (parts.length == 2)
				return "home";
			else if (parts.length == 3)
				return parts[2].substring(0, parts[2].length() - 3);
			else
				return parts[2];
		}
	}
}

