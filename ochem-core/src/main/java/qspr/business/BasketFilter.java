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

package qspr.business;

import qspr.entities.Mapping2;

// Midnighter on Mar 22, 2012
/**
 * Filters for basket queries
 * @author midnighter
 *
 */
public class BasketFilter 
{
	public Long id;
	public String name;
	public boolean showPublicBaskets = false;
	public boolean showSystemBaskets = false;
	public boolean showGroupBaskets = false;

	/**
	 * Baskets containing particular molecule
	 */
	public Mapping2 mapping2;
}
