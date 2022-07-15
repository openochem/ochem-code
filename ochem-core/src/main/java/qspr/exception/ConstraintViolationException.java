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

package qspr.exception;

import java.util.StringTokenizer;

import com.eadmet.exceptions.UserFriendlyException;


public class ConstraintViolationException extends UserFriendlyException 
{
	private static final long serialVersionUID = 1L;
	
	public String table;
	public String constraint;
	public String column;
	
	public String referencedTable;
	public String referencedColumn;
	
	public ConstraintViolationException(org.hibernate.exception.ConstraintViolationException e)
	{
		super(e.getCause().getMessage());
		try{
			StringTokenizer tokenizer = new StringTokenizer(e.getCause().getMessage(), "(");
			tokenizer.nextToken("/");
			table = tokenizer.nextToken("`").substring(1);
			tokenizer.nextToken();
			constraint = tokenizer.nextToken();
			tokenizer.nextElement();
			column = tokenizer.nextToken();
			tokenizer.nextToken();
			referencedTable = tokenizer.nextToken();
			tokenizer.nextToken();
			referencedColumn = tokenizer.nextToken();
		} catch (Exception exc)
		{
			
		}
	}
	
	public String getUserMessage()
	{
		if ("ExperimentalProperty".equals(table) && "exp_property_id_conn".equals(column))
			return "This record is referenced by another record";
		else if ("ArticleUserPdf".equals(table))
			return "Some other user has PDFs attached to target article";
		return getMessage();
	}
}
