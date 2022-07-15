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
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import Jama.Matrix;


public class DOptimal {

	private static transient final Logger logger = LogManager.getLogger(DOptimal.class);

	private List<Integer> fixed = new LinkedList<Integer>();
	private int compoundNr;
	private double[][] descriptors;

	private List<List<Integer>> subsets = new ArrayList<List<Integer>>();
	private List<Integer> optimalSet = new LinkedList<Integer>();
	private double score = 0;

	private int numberOfInitialSubsetsToCreate;
	private int maximumNumberOfOptimizationCycles = 20;


	public DOptimal(List<Integer> fixedCompounds, int numberOfCompoundsToSelect, double[][] principalProperties) {

		if (fixedCompounds == null) {
			this.fixed = new ArrayList<Integer>();
		}
		else {
			this.fixed = fixedCompounds;
		}

		this.compoundNr = numberOfCompoundsToSelect;
		this.descriptors = principalProperties;

		this.numberOfInitialSubsetsToCreate = this.descriptors.length * 100;

		addCombinations();

		double bestInitialScore = 0;
		List<Integer> bestInitialSet = null;

		for (List<Integer> subset : this.subsets) {
			double tempScore = calculateScore(subset);

			if (tempScore > bestInitialScore) {

				bestInitialScore = tempScore;
				bestInitialSet = new LinkedList<Integer>();

				for (Integer com : subset) {
					bestInitialSet.add(com);
				}
			}
		}

		logger.info("3. Scores calculated");

		List<Integer> optimizedSet = optimizedSet(bestInitialSet);
		double cs = calculateScore(optimizedSet);

		if (cs > this.score) {
			this.optimalSet = optimizedSet;
			this.score = cs;
		}
	}

	private void addCombinations() {

		while (this.subsets.size() < this.numberOfInitialSubsetsToCreate) {

			HashMap<Integer, Integer> subsetCompounds = new HashMap<Integer, Integer>();

			for (Integer fixcomp : this.fixed) {
				subsetCompounds.put(fixcomp, 0);
			}

			while(subsetCompounds.size() < this.compoundNr + this.fixed.size()) {
				Integer randomnumber =  (int) (Math.random() * this.descriptors.length);
				subsetCompounds.put(randomnumber, 0);
			}

			List<Integer> combination = new ArrayList<Integer>();
			for (int i : subsetCompounds.keySet()) {
				combination.add(i);
			}

			this.subsets.add(combination);
		}
	}

	/*
	private double calculateScore(List<Integer> sub) {
		double[][] matrix = new double[this.compoundNr + this.fixed.size()][this.descriptors[0].length];
		int ids[] = new int[this.compoundNr + this.fixed.size()];
		int i = 0;

		for (Integer compound : sub) {
			ids[i] = compound;

			for (int j = 0; j < this.descriptors[0].length; j++) {
				matrix[i][j] = this.descriptors[compound][j];
			}
			i++;
		}

		Matrix m = new Matrix(matrix);
		double det = m.transpose().times(m).det();

		return det;
	}
	 */

	private double calculateScore(List<Integer> sub) {

		double[][] matrix = new double[this.compoundNr + this.fixed.size()][this.descriptors[0].length + 1];
		int ids[] = new int[this.compoundNr + this.fixed.size()];

		int i = 0;

		for (Integer compound : sub) {
			ids[i] = compound;

			matrix[i][0] = 1;
			for (int j = 1; j < matrix[0].length; j++) {
				matrix[i][j] = this.descriptors[compound][j-1];
			}
			i++;
		}

		Matrix m = new Matrix(matrix);
		double det = m.transpose().times(m).det();

		return det;
	}

	private List<Integer> optimizedSet(List<Integer> initialSet) {

		List<Integer> disopt = new ArrayList<Integer>();

		for (int compound : initialSet) {
			disopt.add(compound);
		}

		HashMap<Integer, Integer> swaps = new HashMap<Integer, Integer>();

		boolean change = true;
		int init = 0;
		while (change && init < this.maximumNumberOfOptimizationCycles) {

			init++;

			List<Integer> temp = new ArrayList<Integer>();
			List<Integer> refer = new ArrayList<Integer>();

			for (int compound : disopt) {
				temp.add(compound);
				refer.add(compound);
			}

			if (swaps.size() != 1) {
				swaps = new HashMap<Integer, Integer>();
			}

			change = false;

			for (Integer c: refer) {

				if (!(this.fixed.contains(c))) {
					double sc = calculateScore(temp);
					Integer max = c;

					temp.remove(c);
					boolean updated = false;

					for (Integer i = 0; i < this.descriptors.length; i++) {
						if (!(disopt.contains(i))) {

							temp.add(i);
							double tempScore = calculateScore(temp);

							if (tempScore > sc) {
								sc = tempScore;
								max = i;
								updated = true;
							}
							temp.remove(i);
						}
					}

					if (updated && null == swaps.get(max)) {

						disopt.remove(c);
						disopt.add(max);
						temp.add(max);

						swaps.put(c, max);
						change = true;

					}
					else {
						temp.add(c);
					}
				}
			}
		}

		return disopt;
	}

	public List<Integer> getSelectedSet() {
		return this.optimalSet;
	}

	public double getScore() {
		return this.score;
	}

}
