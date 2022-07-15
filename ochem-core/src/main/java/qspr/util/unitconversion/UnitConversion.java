
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

import java.io.ByteArrayInputStream;

import org.antlr.runtime.ANTLRInputStream;
import org.antlr.runtime.CommonTokenStream;

import qspr.entities.Unit;
import qspr.exception.UnitConversionException;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.unitconversion.UnitConversionLexer;
import com.eadmet.unitconversion.UnitConversionParser;
import com.eadmet.utils.NumericalValueStandardizer;

// Just a helper wrapper around antlr parser for unit conversion

public class UnitConversion 
{
	/**
	 * Converts a value into another unit according to a conversion formula.
	 * 
	 * Example: parse('1000*#', 7) converts 7km into 7000m.
	 * 
	 * @param formula The conversion formula, e.g., '-logn(#)'.
	 * @param value The value, given in the original unit.
	 * @return The value, converted to the new unit.
	 */
	public static double parse(String formula, double value, double molweight)
	{
		try
		{
			final ANTLRInputStream input = 
					new ANTLRInputStream(new ByteArrayInputStream(formula.getBytes()));

			// Implement hashing of parsers, if there are performance problems / Midnighter
			UnitConversionParser parser = new UnitConversionParser(
					new CommonTokenStream(new UnitConversionLexer(input)));

			parser.varvalue = value;
			parser.varmw = molweight;

			return parser.start();
		}
		catch (Exception e)
		{
			throw new RuntimeException(e);
		}
	}

	public static double convert(double value, Unit unit, Unit targetUnit, Double molWeight) throws UnitConversionException{
		return convert(value, unit, targetUnit, molWeight, NumericalValueStandardizer.SIGNIFICANT_DIGITS);
	}


	/**
	 *  Converts to default or to provided targetUnit otherwise
	 * @param value
	 * @param unit
	 * @param targetUnit
	 * @param molWeight
	 * @param digits
	 * @return
	 * @throws UnitConversionException
	 */
	public static double convert(double value, Unit unit, Unit targetUnit, Double molWeight, int digits) throws UnitConversionException
	{
		if (targetUnit != null && targetUnit.equals(unit))
			return value;

		if (targetUnit != null && !unit.category.equals(targetUnit.category))
			throw new UserFriendlyException("Trying to convert a unit of " + unit.category.name + "(" + unit.getName() + ") to a unit of " + targetUnit.category.name + "(" + targetUnit.getName() + ")");

		Double res = value;
		double mw = molWeight != null ? molWeight : 1;
		if (unit.isdefault == 0)
			if (unit.toDefaultConversion != null)
				res = UnitConversion.parse(unit.toDefaultConversion, res, mw); // converting to the default unit
			else 
				throw new UnitConversionException(unit, true);

		if (targetUnit != null)
			if (targetUnit.isdefault == 0)
				if (targetUnit.fromDefaultConversion != null)
					res = UnitConversion.parse(targetUnit.fromDefaultConversion, res, mw);
				else
					throw new UnitConversionException(targetUnit, false);

		if (Double.isInfinite(res) || Double.isNaN(res)){
			System.out.println("Infinity or NaN received while converting "+value+" from "+unit.getName() + " to "+targetUnit.getName()+ " for mw="+molWeight);
			return Double.NaN;
		}

		return NumericalValueStandardizer.getSignificantDigitsDouble(res, digits);
	}
}
