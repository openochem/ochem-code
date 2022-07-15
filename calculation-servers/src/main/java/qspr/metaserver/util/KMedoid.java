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


public class KMedoid {
	
	private List<Integer> preselectedSet = new LinkedList<Integer>();
	private int compoundsToSelect;
	private double[][] descriptorMatrix;
	
	private List<Integer> optimalSet = new LinkedList<Integer>();
	
	private double distanceMatrix[][] = null;
	private int numberOfInitialAttempts = 15;
	private double density = 0;
	private HashMap<Integer, List<Integer>> center2Cluster = null;
	
	
	public KMedoid(List<Integer> fixedCompounds, int numberOfCompoundsToSelect, double[][] principalProperties) {
		
		if (fixedCompounds == null) {
			this.preselectedSet = new ArrayList<Integer>();
		}
		else {
			this.preselectedSet = fixedCompounds;
		}
		
		this.compoundsToSelect = numberOfCompoundsToSelect;
		this.descriptorMatrix = principalProperties;
		
		executeSelection();
		
	}
	
	private void executeSelection() {
		
		this.distanceMatrix = ExperimentalDesignHelper.calculateDistanceMatrix(this.descriptorMatrix);
		
		for (int z = 0; z < this.numberOfInitialAttempts; z++) {
			HashMap<Integer, Integer> tempSelection = new HashMap<Integer, Integer>();
			for (Integer compound : this.preselectedSet) {
				tempSelection.put(compound, 1);
			}
			
			while (tempSelection.size() < this.preselectedSet.size() + this.compoundsToSelect) {
				Integer randomnumber =  (int) (Math.random() * this.descriptorMatrix.length);
				if (!this.preselectedSet.contains(randomnumber)) {
					tempSelection.put(randomnumber, 0);
				}
			}
			
			boolean change = true;
			double tempDensity = 0;
			int cycles = 0;
			
			while (change && cycles < 50) {
				cycles++;
				this.center2Cluster = new HashMap<Integer, List<Integer>>();
				for (Integer compound : tempSelection.keySet()) {
					this.center2Cluster.put(compound, new ArrayList<Integer>());
				}
				
				for (Integer i  = 0; i < this.descriptorMatrix.length; i++) {
					double minDist = 1000000000;
					Integer cluster = 1000000000;
					
					for (Integer s : tempSelection.keySet()) {
						if (this.distanceMatrix[i][s] < minDist || cluster.intValue() == 1000000000) {
							if (!(this.center2Cluster.containsKey(i) && i.intValue() != s.intValue())) {
								minDist = this.distanceMatrix[i][s];
								cluster = s;
							}
						}
					}
					
					List<Integer> temp = this.center2Cluster.get(cluster);
					temp.add(i);
					this.center2Cluster.put(cluster, temp);
				}
				
				HashMap<Integer, Integer> newSelection = new HashMap<Integer, Integer>();
				
				change = false;
				
				List<Integer> toBeDumped = new ArrayList<Integer>();
				for (Integer center : tempSelection.keySet()) {
					toBeDumped.add(center);
					if (tempSelection.get(center).intValue() == 1) {
						newSelection.put(center, 1);
					}
					else {
						Integer newCenter = getClusterCenter(this.center2Cluster.get(center));
						if (newCenter.intValue() != center.intValue()) {
							change = true;
						}
						newSelection.put(newCenter, 0);
					}
				}
				
				tempSelection = newSelection;
				tempDensity = calculateDensity(this.center2Cluster);
				
			}
			
			if (this.density == 0 || tempDensity < this.density) {
				
				this.density = tempDensity;
				this.optimalSet = new ArrayList<Integer>();
				
				for (Integer compound : tempSelection.keySet()) {
					if (tempSelection.get(compound).intValue() != 1) {
						this.optimalSet.add(compound);
					}
				}
			}
		}
		System.out.print(this.density);
		
	}
	
	private Integer getClusterCenter(List<Integer> cluster) {
		Integer center = 100000;
		double minDistance = 1000000000;
		
		for (Integer i : cluster) {
			double tempDistance = 0;
			for (Integer j : cluster) {
				tempDistance += this.distanceMatrix[i][j];
			}
			
			if (tempDistance < minDistance || center.intValue() == 100000) {
				minDistance = tempDistance;
				center = i;
			}
		}
		
		return center;
	}
	
	private double calculateDensity(HashMap<Integer, List<Integer>> clustering) {
		double distance = 0;
		
		for (Integer c : clustering.keySet()) {
			List<Integer> cluster = clustering.get(c);
			for (Integer i : cluster) {
				distance += this.distanceMatrix[i][c];
			}
		}
		
		return distance;
	}
	
	public List<Integer> getSelectedSet() {
		return this.optimalSet;
	}
	
}
