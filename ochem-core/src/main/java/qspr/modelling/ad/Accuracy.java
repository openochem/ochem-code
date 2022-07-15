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

package qspr.modelling.ad;

import qspr.entities.Property;

// A class for accuracy averaging / Midnighter
public abstract class Accuracy
{
	public double sum = 0;
	public int count = 0;
	abstract double getValue(double real, double predicted);
	abstract public double getAverageValue();
	abstract public boolean isBetterThan(Accuracy accuracy);
	
	public void addValue(double real, double predicted)
	{
		sum += getValue(real, predicted);
		count++;
	}
	
	public void removeValue(double real, double predicted)
	{
		sum -= getValue(real, predicted);
		count--;
	}
	
	public static Accuracy forProperty(Property p)
	{
		return p.isNumeric() ? new RMSE() : new ClassificationAccuracy();
	}
}

class RMSE extends Accuracy
{
	@Override
	public double getAverageValue() {
		return Math.sqrt(sum / count);
	}

	@Override
	double getValue(double real, double predicted) {
		return (real - predicted) * (real - predicted);
	}

	@Override
	public boolean isBetterThan(Accuracy accuracy) {
		return getAverageValue() < accuracy.getAverageValue();
	}
}

class ClassificationAccuracy extends Accuracy
{
	@Override
	public double getAverageValue() {
		return (count - sum) / count;
	}

	@Override
	double getValue(double real, double predicted) {
		return 0 == Math.round(real - predicted) ? 0 : 1;
	}
	
	@Override
	public boolean isBetterThan(Accuracy accuracy) {
		return getAverageValue() > accuracy.getAverageValue();
	}
}
