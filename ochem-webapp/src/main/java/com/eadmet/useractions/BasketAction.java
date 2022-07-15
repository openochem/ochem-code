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

package com.eadmet.useractions;

import javax.xml.bind.annotation.XmlRootElement;

import qspr.entities.Basket;

@XmlRootElement
public class BasketAction extends AbstractUserAction
{
	Basket basket;
	int recordsCount;
	BasketActionType type;
	
	@Override
	public String getLogLine()
	{
		switch (type)
		{
		case ADDENTRY:
			return String.format(" has added %d records to basket \"%s\"", recordsCount, basket.name);
		case REMOVEENTRY:
			return String.format(" has removed %d records from basket \"%s\"", recordsCount, basket.name);
		case DELETE_BASKET:
			return "has deleted a basket \"" + basket.name + "\"";
		case RENAME:
			return "has renamed a basket to \"" + basket.name + "\"";
			
		}
		return null;
	}
	
	public BasketAction()
	{
		
	}
	
	public BasketAction(Basket basket, int recordsCount, BasketActionType type)
	{
		this.basket = new Basket();
		this.basket.id = basket.id;
		this.basket.name = basket.name;
		this.recordsCount = recordsCount;
		this.type = type;
	}
	
	public static enum BasketActionType {ADDENTRY, REMOVEENTRY, DELETE_BASKET, RENAME} ; 
	
}
