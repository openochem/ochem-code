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

package qspr.metaserver.util;
import java.lang.reflect.InvocationTargetException;

import org.apache.commons.lang.NotImplementedException;

import com.eadmet.exceptions.UserFriendlyException;

import qspr.dao.ChemInfEngine;
import qspr.dao.Various;

/**
 * An extended OCHEM SMART pattern with NOT and AND keywords 
 * Example: CCC AND NOT N
 * 
 *  * @author midnighter
 */
public abstract class ExtendedSMART
{
	public boolean invalid;
	private MatchablePattern pattern;
	private String singleSmart;
	
	public static ExtendedSMART create(String smarts, ChemInfEngine engine) {
		try {
			switch (engine) {
			case CHEMAXON:
				return (ExtendedSMART) Class.forName("qspr.metaserver.util.ExtendedSMARTChemaxon").getConstructor(String.class).newInstance(smarts);
			case CDK:
				return (ExtendedSMART) Class.forName("qspr.metaserver.util.ExtendedSMARTCDK").getConstructor(String.class).newInstance(smarts);
			default:
				throw new NotImplementedException("Engine not found:" + engine);
			}
		} catch (InstantiationException | IllegalAccessException | IllegalArgumentException
				| InvocationTargetException | NoSuchMethodException | SecurityException
				| ClassNotFoundException e) {
			throw new UserFriendlyException(e);
		}
	}

	public boolean match(String sdf) throws Exception
	{
		if (invalid)
			return false;

		return pattern.matches(sdf);
	}

	public int getMatchCount(String sdf) throws Exception
	{
		if (invalid)
			return 0;

		if (singleSmart == null)
			return (pattern.matches(sdf)) ? 1 : 0;

		return getMatchCount(sdf, singleSmart);
	}

	protected abstract int getMatchCount(String molecule, String smart) throws Exception;

	public abstract boolean match(final String molecule, final String smart) throws Exception;

	public ExtendedSMART(String exSmart)
	{
		Junction disjunction = new Junction(this, true);
		try
		{
			disjunction.parse(exSmart);
			if (disjunction.operands[0] instanceof SimpleSMART && disjunction.operands.length == 1)
				singleSmart = ((SimpleSMART) disjunction.operands[0]).smart;
			pattern = disjunction;
		}
		catch (Exception e)
		{
			invalid = true;
			//e.printStackTrace();
		}
	}
}

interface MatchablePattern
{
	public boolean matches(String molecule) throws Exception;
}

class Inversion implements MatchablePattern
{
	public MatchablePattern operand;

	@Override
	public boolean matches(String molecule) throws Exception 
	{
		return !operand.matches(molecule);
	}

	public Inversion(MatchablePattern pattern)
	{
		operand = pattern;
	}
}

class Junction implements MatchablePattern
{
	public MatchablePattern[] operands;

	/**
	 * Is it a disjunction or a conjunction?
	 */
	boolean or = false;
	
	ExtendedSMART impl;
	

	public Junction(ExtendedSMART impl, boolean or)
	{
		this.or = or;
		this.impl = impl;
	}

	@Override
	public boolean matches(String molecule) throws Exception 
	{
		for (int i = 0; i < operands.length; i++)
		{
			if (or && operands[i].matches(molecule))
				return true;

			if (!or && !operands[i].matches(molecule))
				return false;
		}

		return !or;
	}

	public void parse(String st)
	{
		String[] parts = st.split(or ? " OR " : " AND ");
		operands = new MatchablePattern[parts.length];
		int i = 0;
		for (String part : parts) {
			part = part.trim();
			if (part.contains(" OR "))
			{
				Junction disj = new Junction(impl, true);
				disj.parse(part);
				operands[i] = disj;
			}
			else if (part.contains(" AND "))
			{
				Junction conj = new Junction(impl, false);
				conj.parse(part);
				operands[i] = conj;
			}
			else if (part.toLowerCase().startsWith("not"))
				operands[i] = new Inversion(new SimpleSMART(impl, part.substring(3).trim()));
			else if (part.contains("<") || part.contains(">"))
			{
				Predicate predicate = new Predicate(impl);
				predicate.parse(part);
				operands[i] = predicate;
			}
			else
				operands[i] = new SimpleSMART(impl, part);

			i++;
		}
	}
}

class Predicate implements MatchablePattern
{
	public boolean mw;
	public String smart;
	public String predicate;
	public Integer value;
	
	ExtendedSMART impl;
	
	public Predicate(ExtendedSMART impl) {
		this.impl = impl;
	}

	@Override
	public boolean matches(String molecule) throws Exception 
	{
		Integer count;
		if (mw)
			count = new Double(Various.molecule.getMass(molecule)).intValue();
		else
			count = impl.getMatchCount(molecule, smart);
		if ("<".equals(predicate))
			return count < value;
		else if (">".equals(predicate))
			return count > value;
			else if ("=".equals(predicate))
				return count == value;
		return false;
	}

	public void parse(String st)
	{
		String[] parts = null;
		String[] predicates = new String[]{"<", ">"};
		for (int i = 0; i < predicates.length; i++)
			if (st.contains(predicates[i]))
			{
				predicate = predicates[i];
				parts = st.split(predicate);
			}

		String o1 = parts[0].trim();
		String o2 = parts[1].trim();
		if (o1.equals("MW"))
			mw = true;
		else
			smart = o1;
		value = Integer.valueOf(o2.trim());
	}
}

class SimpleSMART implements MatchablePattern
{
	public String smart;
	ExtendedSMART impl;

	@Override
	public boolean matches(String molecule) throws Exception 
	{
		return impl.match(molecule, smart);
	}

	public SimpleSMART(ExtendedSMART impl, String smartStr)
	{
		smart = smartStr;
		this.impl = impl;
	}
}
