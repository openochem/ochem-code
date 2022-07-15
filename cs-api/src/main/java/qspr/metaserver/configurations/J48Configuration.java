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

package qspr.metaserver.configurations;

import javax.xml.bind.annotation.XmlRootElement;

import qspr.workflow.utils.QSPRConstants;

@XmlRootElement(name = "j48-configuration")
public class J48Configuration extends AbstractWekaConfiguration
{
	private static final long serialVersionUID = 1L;

	// -U, Use unpruned tree.
	public Boolean useUnpruned = false;

	// -C confidence, set confidence threshold for pruning. (Default: 0.25)
	public Double confidence = 0.25;

	// -M number, set minimum number of instances per leaf. (Default: 2)
	public Integer instancesPerLeaf = 2;

	// -R, use reduced error pruning. No subtree raising is performed.
	public Boolean reducedPruning = false;

	// -N number, set number of folds for reduced error pruning. One fold is
	// used as the pruning set. (Default: 3)
	public Integer numReducedPruningFolds = 3;

	// -B, use binary splits for nominal attributes.
	public Boolean useBinarySplits = false;

	// -S, son't perform subtree raising.
	public Boolean dontPerformRaising = false;

	// -L, do not clean up after the tree has been built.
	public Boolean noCleanup = false;

	// -A, if set, Laplace smoothing is used for predicted probabilites.
	public Boolean useLaplaseSmoothing = false;

	public String toString() 
	{
		StringBuilder sb = new StringBuilder();
		if (useUnpruned != null && useUnpruned)
			sb.append(", use unpruned");

		if (confidence != null)
			sb.append(", confidence = "+confidence);

		if (instancesPerLeaf != null)
			sb.append(", inst per leaf = "+instancesPerLeaf);

		if (reducedPruning != null && reducedPruning)
			if (numReducedPruningFolds != null)
				sb.append(", reduced pruning with folds = "+numReducedPruningFolds);

		if (useBinarySplits != null && useBinarySplits)
			sb.append(", bin splits");

		if (dontPerformRaising != null && dontPerformRaising)
			sb.append(", no raising");

		if (noCleanup != null && noCleanup)
			sb.append(", no cleanup");

		if (useLaplaseSmoothing != null && useLaplaseSmoothing)
			sb.append(", Laplase smoothing");

		if (sb.length() > 2)
			return "\n"+sb.substring(2, sb.length());
		else
			return "";
	}

	@Override
	public String getDefaultName() {
		return QSPRConstants.J48;
	}
}
