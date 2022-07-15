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

import qspr.metaserver.configurations.J48Configuration;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.workflow.utils.QSPRConstants;

public class J48Server extends WekaAbstractServer
{

	@Override
	protected void writeMethodSpecificCfgXML(BufferedWriter out, ModelAbstractConfiguration conf) throws Exception
	{
		J48Configuration receivedConf = (J48Configuration) conf;
		if (receivedConf.useUnpruned != null && receivedConf.useUnpruned)
		{
			out.write("<option type=\"flag\" name=\"U\"/>\n");
		} else
		{
			if (receivedConf.reducedPruning != null && receivedConf.reducedPruning)
			{
				out.write("<option type=\"flag\" name=\"R\"/>\n");
				
				if (receivedConf.numReducedPruningFolds != null)
					out.write("<option type=\"single\" name=\"N\">" + receivedConf.numReducedPruningFolds + "</option>\n");
				
				out.write("<option type=\"single\" name=\"Q\">" + receivedConf.getSeed() + "</option>\n");
			} else
			{
				out.write("<option type=\"single\" name=\"C\">" + receivedConf.confidence + "</option>\n");
			}
			
			if (receivedConf.instancesPerLeaf != null)
				out.write("<option type=\"single\" name=\"M\">" + receivedConf.instancesPerLeaf + "</option>\n");
		}
		
		if (receivedConf.useBinarySplits != null && receivedConf.useBinarySplits)
			out.write("<option type=\"flag\" name=\"B\"/>\n");	
		
		if (receivedConf.useLaplaseSmoothing != null && receivedConf.useLaplaseSmoothing)
			out.write("<option type=\"flag\" name=\"A\"/>\n");		

		if (receivedConf.noCleanup != null && receivedConf.noCleanup)
			out.write("<option type=\"flag\" name=\"L\"/>\n");
		
		if (receivedConf.dontPerformRaising != null && receivedConf.dontPerformRaising)
			out.write("<option type=\"flag\" name=\"S\"/>\n");	
	}
	
	public J48Server(){
		supportedTaskType = QSPRConstants.J48;
		className = "weka.classifiers.trees.J48";
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);
	}

}
