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

package qspr.entities;

import java.io.Serializable;
import java.util.List;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlRootElement;

import qspr.Globals;
import qspr.annotations.Loggable;

/**
 * 
 * Alert substitution variable (aka "vector binding").
 * Intended to simplify complex SMARTS patterns and make them more interpretable
 *  
 * @author midnighter
 *
 */
@Entity
@Loggable
@XmlRootElement(name = "alert-variable")
public class AlertSubstitutionVariable implements Serializable
{
	private static final long serialVersionUID = 1L;

	@Id
	@GeneratedValue
	@Column(name = "asv_id")
	@XmlAttribute
	public Long id;

	@Column
	public String name;

	@Column(name = "variable_name")
	public String variableName;

	@Column
	public String substitution;

	@Column
	public String description;

	private static List<AlertSubstitutionVariable> allVariables = null;
	@SuppressWarnings("unchecked")
	public static String substituteVariables(String smarts)
	{
		if (allVariables == null)
			allVariables = Globals.session().createCriteria(AlertSubstitutionVariable.class).list();

		for (AlertSubstitutionVariable variable : allVariables)
		{
			// Allow the user to use both $Hal and [$Hal], sometimes its convenient to have brackets around the variable
			if (variable.substitution.startsWith("["))
				smarts = smarts.replace((CharSequence)("[$" + variable.variableName + "]"), (CharSequence)(variable.substitution));
			smarts = smarts.replace((CharSequence)("$" + variable.variableName), (CharSequence)(variable.substitution));
		}

		return smarts;
	}
}
