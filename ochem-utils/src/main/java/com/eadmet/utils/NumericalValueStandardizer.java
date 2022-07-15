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

package com.eadmet.utils;

import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;

public class NumericalValueStandardizer
{

	public static final int SIGNIFICANT_DIGITS = 4; // sets the default number of significant digits
	public static final int SIGNIFICANT_DIGITS_CANONICAL = 3; // sets the default number of significant digits

	private static Map<Integer, DecimalFormat> formats = new HashMap<Integer,DecimalFormat>();

	/*
	 * By default we will have SIGNIFICANT_DIGITS digits
	 */
	public static String getSignificantDigits(double value)
	{
		return value == 0? "0" : value == 1? "1": getSignificantDigitsStr(value, SIGNIFICANT_DIGITS);
	}

	public static String getWithConfidenceIntervals(double value, double precision)
	{
		if(value == 0)return "0";
		if(value == 1)return "1";

		if (precision == 0)
			return getSignificantDigitsStr(value, 2);

		long n1 = Math.round(0.4999 + Math.log10(Math.abs(value)));
		long n2 = Math.round(0.4999 + Math.log10(Math.abs(precision)));

		int digits = (int) (1 + n1 - n2);

		return getSignificantDigitsStr(value, digits);
	}

	/**
	 * Get a string representation of a number using the provided number of significant digits
	 */
	public static String getSignificantDigitsStr(double value, int significant_digits)
	{
		if(value == 0)return "0";
		if(value == 1)return "1";

		if (Double.isNaN(value))
			return "NaN";

		if (Double.isInfinite(value))
			return "Inf";

		DecimalFormat DoubleFormat = getFormat(significant_digits);

		double val = Double.parseDouble(DoubleFormat.format(value));
		String s = "" + val;
		if (s.contains("E"))
			return s;
		if (val == (Math.round(val)))
			return "" + Math.round(val);
		while (s.endsWith("0"))
			s = s.substring(0, s.length() - 1);

		return s;
	}

	/**
	 * Creates and caches formats instead of creating them again and again
	 * @param significant_digits
	 * @return
	 */
	static DecimalFormat getFormat(Integer significant_digits){
		if(!formats.containsKey(significant_digits)){

			String format = "0.";
			for (int i = 1; i < significant_digits; i++)
				format += "0";
			format += "E0";

			Locale.setDefault(Locale.US);
			DecimalFormat DoubleFormat = new DecimalFormat(format);
			formats.put(significant_digits, DoubleFormat);
		}

		return formats.get(significant_digits);
	}

	/**
	 * Round a double value using the specified number of significant digits.
	 * Uses a string conversion method and converts it back to a double value.
	 * 
	 */
	public static double getSignificantDigitsDouble(double value, int significant_digits)
	{
		if(value == 0)return 0;
		if(value == 1)return 1;

		try
		{
			return Double.parseDouble(getSignificantDigitsStr(value, significant_digits));
		} catch (Exception e)
		{
			return value;
		}
	}

	// NoS 28.08.2011
	// Keeps two non-zero signs after comma, if they are there, more - if needed
	// 10.002 -> 10, 10.023 -> 10.02, 0.0001234 -> 0.00012
	// Performance - 250K floats per second. Could be better...
	// Used in WEKA and LibSVM servers
	public static String formattedFloatValuesForWEKA(double value)
	{
		if (value == 0) return "0";
		if (value == 1) return "1";

		int pow = 2 - Math.min(0, (int) Math.floor(Math.log10(Math.abs(value))));
		String format = "#.";
		for (int i = 0; i < pow; i++)
			format += "#";
		DecimalFormat myFormatter = new DecimalFormat(format);
		return myFormatter.format(value);
	}

	public static void main(String[] args)
	{
		for (int i = 0; i < 15; i++)
		{
			double val = 44.09 * 100000000 / Math.pow(10, i);
			double prec = 3 * 100000000 / Math.pow(10, i);
			System.out.println("" + val + " " + prec + " --> " + NumericalValueStandardizer.getWithConfidenceIntervals(val, prec) + "  +- "
					+ NumericalValueStandardizer.getSignificantDigitsStr(prec, 0));
		}

		double val = 44.89;
		double prec = 0.3;

		System.out.println("" + val + " " + prec + " --> " + NumericalValueStandardizer.getWithConfidenceIntervals(val, prec) + "  +- "
				+ NumericalValueStandardizer.getSignificantDigitsStr(prec, 0));

		for (int i = 100; i < 1000; i+=10)
		{
			val = i + 0.4;
			System.out.println("" + val + " --> " + NumericalValueStandardizer.getSignificantDigits(val)
);
		}
		
	}
}
