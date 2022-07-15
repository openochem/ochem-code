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

package qspr.metaserver.cs;

import java.io.BufferedWriter;

import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.configurations.RFConfiguration;
import qspr.workflow.utils.QSPRConstants;

public class RFServer extends WekaAbstractServer
{

	@Override
	protected void writeMethodSpecificCfgXML(BufferedWriter out, ModelAbstractConfiguration conf) throws Exception
	{
		RFConfiguration receivedConf = (RFConfiguration) conf;
		out.write("<option type=\"single\" name=\"I\">" + receivedConf.numTrees + "</option>\n");
		out.write("<option type=\"single\" name=\"K\">" + receivedConf.numFeatures + "</option>\n");
		out.write("<option type=\"single\" name=\"depth\">" + receivedConf.maxDepth + "</option>\n");
		out.write("<option type=\"single\" name=\"S\">" + conf.getSeed() + "</option>\n");
	}

	public RFServer(){
		supportedTaskType = QSPRConstants.RF;
		className = "weka.classifiers.trees.RandomForest";
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);
	}

}
