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

package qspr.oscript;

import java.io.IOException;

import org.antlr.runtime.RecognitionException;
import org.junit.Assert;
import org.junit.Test;

import com.eadmet.exceptions.UserFriendlyException;

public class OScriptTest 
{
	@Test
	public void basicTest() throws IOException, RecognitionException
	{
		TransformationRule rule = TransformationRule.compile("// Comment\n*[Pressure=10 atm] -> logPow\n //One more comment\n*[Pressure = 10 atm] -> logPow");
		checkBasicRule(rule.rules.get(0));
		checkBasicRule(rule.rules.get(1));
		
		
	}
	
	@Test(expected = UserFriendlyException.class)
	public void incorrectSyntaxTest() throws RecognitionException, IOException
	{
		TransformationRule.compile("aa - bb");
	}
	
	@Test(expected = UserFriendlyException.class)
	public void incorrectSyntaxTest2() throws RecognitionException, IOException
	{
		TransformationRule.compile("*[Pressure=] -> logPow");
	}
	
	private void checkBasicRule(AtomicRule rule)
	{
		Assert.assertEquals("*", rule.left.propertyName);
		Assert.assertEquals("Pressure", rule.left.conditionName);
		Assert.assertEquals("atm", rule.left.unitName);
		Assert.assertEquals("10", rule.left.value);
		Assert.assertEquals("logPow", rule.right.propertyName);
	}
}
