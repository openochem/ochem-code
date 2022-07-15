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

import java.util.ArrayList;
import java.util.List;

import qspr.dao.ChemInfEngine;

/**
 * An OCHEM interface for matching multiple extended SMARTS
 * @author midnighter
 * 
 */
public class SMARTMatcher 
{
	public List<ExtendedSMART> extendedSmarts = new ArrayList<ExtendedSMART>();
	private ChemInfEngine engine;

	public static ExtendedSMART getExtended(String pattern, ChemInfEngine engine) {
		return ExtendedSMART.create(pattern, engine);
	};

	public SMARTMatcher(List<String> patterns, ChemInfEngine engine)
	{
		this.engine = engine;
		setPatterns(patterns);
	}

	public SMARTMatcher(String pattern, ChemInfEngine engine)
	{
		this.engine = engine;
		List<String> patterns = new ArrayList<String>();
		patterns.add(pattern);
		setPatterns(patterns);
	}

	public void setPatterns(List<String> patterns)
	{
		String pattern = null;
		for (int i = 0; i < patterns.size(); i++) 
		{
			pattern = patterns.get(i);
			extendedSmarts.add(getExtended(pattern, engine));
		}
	}

	public void addPattern(String pattern)
	{
		extendedSmarts.add(getExtended(pattern, engine));
	}

	public boolean matchPattern(String sdf, int patternNum) throws Exception 
	{
		return extendedSmarts.get(patternNum).match(sdf);
	}

	public int getMatchCount(String sdf, int patternNum) throws Exception 
	{
		return extendedSmarts.get(patternNum).getMatchCount(sdf);
	}
}


