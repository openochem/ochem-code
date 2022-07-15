
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

package qspr.util.unitconversion;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class UnitConversionTest {

	private void doTest(String formula, double value, double expected, double delta) 
	{
		assertEquals(expected, UnitConversion.parse(formula, value, 0), delta);
	}
	
	private static final double nn = -9999;
	
	//  ///////////////////////////////////
	//  //  Simple one-term expressions  //
	//  ///////////////////////////////////
	
	// Numerical literals.
	@Test public void testA001() { doTest("1"   , nn,  1.0, 0); }
	@Test public void testA002() { doTest("1.1" , nn,  1.1, 0); }
	@Test public void testA003() { doTest("-1"  , nn, -1.0, 0); }
	@Test public void testA004() { doTest("-1.7", nn, -1.7, 0); }
	
	// Constants.
	@Test public void testA101() { doTest("E"   , nn, Math.E       , 0); }
	@Test public void testA102() { doTest("PI"  , nn, Math.PI      , 0); }
	@Test public void testA103() { doTest("NA"  , nn, 6.02214179e23, 0); }
	
	// Functions.
	@Test public void testA201() { doTest("cos(0)"   , nn,  1., 1e-12); }
	@Test public void testA202() { doTest("cos(PI)"  , nn, -1., 1e-12); }
	@Test public void testA203() { doTest("cos(PI/2)", nn,  0., 1e-12); }
	
	@Test public void testA211() { doTest("sin(0)"   , nn,  0., 1e-12); }
	@Test public void testA212() { doTest("sin(PI)"  , nn,  0., 1e-12); }
	@Test public void testA213() { doTest("sin(PI/2)", nn,  1., 1e-12); }
	
	@Test public void testA221() { doTest("sqrt(0)", nn,  0., 0); }
	@Test public void testA222() { doTest("sqrt(1)", nn,  1., 0); }
	@Test public void testA223() { doTest("sqrt(9)", nn,  3., 0); }
	
	@Test public void testA231() { doTest("exp(0)",  nn,  1.                , 1e-12); }
	@Test public void testA232() { doTest("exp(1)",  nn,  Math.E            , 1e-12); }
	@Test public void testA233() { doTest("exp(-1)", nn,  1./Math.E         , 1e-12); }
	@Test public void testA234() { doTest("exp(10)", nn,  22026.465794806717, 1e-12); }
	
	@Test public void testA241a() { doTest("lg(1)"  , nn,  0., 0); }
	@Test public void testA242a() { doTest("lg(10)" , nn,  1., 0); }
	@Test public void testA243a() { doTest("lg(100)", nn,  2., 0); }
	
	@Test public void testA241b() { doTest("log10(1)"  , nn,  0., 0); }
	@Test public void testA242b() { doTest("log10(10)" , nn,  1., 0); }
	@Test public void testA243b() { doTest("log10(100)", nn,  2., 0); }
	
	@Test public void testA251() { doTest("logn(1)"  , nn,  0., 0); }
	@Test public void testA252() { doTest("logn(E)"  , nn,  1., 0); }
	@Test public void testA253() { doTest("logn(E*E)", nn,  2., 0); }
	
	@Test public void testA255() { doTest("ln(1)"  , nn,  0., 0); }
	@Test public void testA256() { doTest("ln(E)"  , nn,  1., 0); }
	@Test public void testA257() { doTest("ln(E*E)", nn,  2., 0); }

	@Test public void testA261() { doTest("pow(1,0)" , nn,  1., 0); }
	@Test public void testA262() { doTest("pow(1,1)" , nn,  1., 0); }
	@Test public void testA263() { doTest("pow(1,2)" , nn,  1., 0); }
	@Test public void testA264() { doTest("pow(2,1)" , nn,  2., 0); }
	@Test public void testA265() { doTest("pow(-1,0)", nn,  1., 0); }
	@Test public void testA266() { doTest("pow(-1,1)", nn, -1., 0); }
	@Test public void testA267() { doTest("pow(-1,2)", nn,  1., 0); }
	@Test public void testA268() { doTest("pow(-1,3)", nn, -1., 0); }
	@Test public void testA269() { doTest("pow(2,3)" , nn,  8., 0); }
	@Test public void testA270() { doTest("pow(10,#)", 3,  1000, 0); }
	
	@Test public void testA281() { doTest("log2(1)" , nn, 0., 0); }
	@Test public void testA282() { doTest("log2(2)" , nn, 1., 0); }
	@Test public void testA283() { doTest("log2(4)" , nn, 2., 0); }
	
	@Test public void testA291() { doTest("log(2,4) "     , nn, 2., 0); }
	@Test public void testA292() { doTest("log(10,1)"     , nn, 0., 0); }
	@Test public void testA293() { doTest("log(PI,PI)"    , nn, 1., 0); }
	@Test public void testA294() { doTest("log(37,37*37)" , nn, 2., 0); }
	
	//  ////////////////////////////
	//  //  Variable expressions  //
	//  ////////////////////////////
	
	@Test public void testB001() { doTest("#"    ,  0,  0., 0); }
	@Test public void testB002() { doTest("#"    ,  1,  1., 0); }
	@Test public void testB003() { doTest("#"    , -1, -1., 0); }
	
	@Test public void testB011() { doTest("2*#"  ,  0,  0., 0); }
	@Test public void testB012() { doTest("2*#"  ,  1,  2., 0); }
	@Test public void testB013() { doTest("2*#"  , -1, -2., 0); }
	
	@Test public void testB021() { doTest("#+#"  ,  0,  0., 0); }
	@Test public void testB022() { doTest("#+#"  ,  1,  2., 0); }
	@Test public void testB023() { doTest("#+#"  , -1, -2., 0); }
	
	@Test public void testB031() { doTest("-#"   ,  1, -1., 0); }
	@Test public void testB032() { doTest("#-#"  ,  3,  0., 0); }
	@Test public void testB033() { doTest("-#+#" , -7,  0., 0); }
	
	//  /////////////////////////////
	//  // Arithmetic expressions  //
	//  /////////////////////////////
	
	// Unary expressions.
	@Test public void testC001() { doTest("-1"     ,  nn, -1.     , 0); }
	@Test public void testC002() { doTest("-PI"    ,  nn, -Math.PI, 0); }
	@Test public void testC003() { doTest("-lg(10)",  nn, -1.     , 0); }
	@Test public void testC004() { doTest("-#"     ,  -1,  1.     , 0); }
	
	// Addition.
	@Test public void testC101() { doTest("1+1"  , nn, 2. , 0); }
	@Test public void testC102() { doTest("1-1"  , nn, 0. , 0); }
	@Test public void testC103() { doTest("1+1+1", nn, 3. , 0); }
	@Test public void testC104() { doTest("1-2+3", nn, 2. , 0); }
	@Test public void testC105() { doTest("1+2"  , nn, 3. , 0); }
	@Test public void testC106() { doTest("2+1"  , nn, 3. , 0); }
	
	// Multiplication.
	@Test public void testC201() { doTest("1*0"    , nn,  0.  , 1e-12); }
	@Test public void testC202() { doTest("1*1"    , nn,  1.  , 1e-12); }
	@Test public void testC203() { doTest("-1*1"   , nn, -1.  , 1e-12); }
	@Test public void testC204() { doTest("2*3"    , nn,  6.  , 1e-12); }
	@Test public void testC205() { doTest("3*2"    , nn,  6.  , 1e-12); }
	@Test public void testC206() { doTest("2*3*0.1", nn,  0.6 , 1e-12); }
	
	// Division.
	@Test public void testC301() { doTest("8/4" , nn,   2.0 , 1e-12); }
	@Test public void testC302() { doTest("-8/2", nn,  -4.0 , 1e-12); }

	
	//  ////////////////
	//  //  Brackets  //
	//  ////////////////
	
	@Test public void testD001() { doTest("2+(2*3)"  , nn,  8.0 , 1e-12); }
	@Test public void testD002() { doTest("-2*(8/4.)", nn, -4.0 , 1e-12); }

	//  /////////////////////////
	//  //  Mixed expressions  //
	//  /////////////////////////
	
	@Test public void testE001() { doTest("log10(pow(10,7))", nn,  7.0 , 1e-12); }
	@Test public void testE002() { doTest("exp(logn(33))"   , nn, 33.0 , 1e-12); }
	@Test public void testE003() { doTest("1+2*3+4/5-(-6)"  , nn, 69./5, 1e-12); }
	@Test public void testE004() { doTest("(1+3)/(3-1)"     , nn,  2.0 , 1e-12); }
	@Test public void testE005() { doTest("pow(pow(2,2),2)" , nn, 16.0 , 1e-12); }
	@Test public void testE006() { doTest("#*#/#"           , 2,  2.0 , 1e-12); }
	@Test public void testE007() { doTest("77*(1/2)"        , nn, 77./2, 1e-12); }

	//  ///////////////////
	//  //  Conversions  //
	//  ///////////////////
	
	@Test public void testF001() { doTest("1000*#"   ,     3, 3000.             , 1e-12); }  // 3km -> 3000m.
	@Test public void testF002() { doTest("-log10(#)", 25e-6, 4.6020599913279624, 1e-12); }  // 25microMolar EC50 -> 4.6microMolar pEC50.
}
