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

package qspr.dao;

import java.util.ArrayList;
import java.util.List;

import qspr.Globals;
import qspr.entities.Property;
import qspr.entities.PropertyOption;

public class PropertyOptionDAOImpl implements PropertyOptionDAO {

	@Override
	public PropertyOption getPropertyOptionById(Long optionId) {
		return (PropertyOption) Globals.session().get(PropertyOption.class, optionId);
	}

	@Override
	public PropertyOption getPropertyOptionByName(String name, Long propertyId, boolean createIfMissing, boolean save) {
		if(name == null) return null;
		name = name.trim();
		if(name.length() == 0) return null;
		Property property = Repository.property.getPropertyById(propertyId);
		PropertyOption option = property.getOption(name);
		if(option != null) return option;
		if(createIfMissing) {
			option = new PropertyOption(name, property);
			if(save)Globals.session().saveOrUpdate(option);
		}
		return option;
	}

	void heapPermutation(String a[], int size, int n, List<String[]> permutations)
	{
		// if size becomes 1 then prints the obtained
		// permutation
		if (size == 1) {
			permutations.add(a.clone());
			return;
		}

		for (int i = 0; i < size; i++) {
			heapPermutation(a, size - 1, n, permutations);

			// if size is odd, swap 0th i.e (first) and
			// (size-1)th i.e (last) element
			if (size % 2 == 1)
				swap(a, 0,size - 1);

			// If size is even, swap ith and
			// (size-1)th i.e (last) element
			else
				swap(a, i, size - 1);
		}
	}

	private void swap(String[] a, int i, int j) {
		String aa = a[i];
		a[i] = a[j];
		a[j] = aa;
	}

	@Override
	public PropertyOption searchPropertyOptionByPermutations(String name, Long propertyId) {
		if(name == null)return null;
		String names[] = name.split(";");
		if(names.length == 1) return getPropertyOptionByName(name, propertyId, false, false);
		for(int i=0;i<names.length;i++)
			names[i]=names[i].trim();
		List<String[]> permutations = new ArrayList<String[]>();
		heapPermutation(names, names.length, names.length, permutations);
		Property property = Repository.property.getPropertyById(propertyId);
		for(String nn[]:permutations){
			String n = "";
			for(int i=0;i<nn.length;i++)
				n = i == 0?nn[0]:n + ";" +nn[i];
			PropertyOption option = property.getOption(n);
			if(option != null) return option;
		}

		return null;
	}

	public static void main(String[] args){
	
		String name="a;b;c";
		String names[] = name.split(";");
		List<String[]> permutations = new ArrayList<String[]>();
		
		PropertyOptionDAOImpl aa = new PropertyOptionDAOImpl();
		aa.heapPermutation(names, names.length, names.length, permutations);
		
		for(String nn[]:permutations){
			String n = "";
			for(int i=0;i<nn.length;i++)
				n = i == 0?nn[0]:n + ";" +nn[i];
			System.out.println(n);
		}
	}

}
