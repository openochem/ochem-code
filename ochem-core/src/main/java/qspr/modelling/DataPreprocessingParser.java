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

package qspr.modelling;

import javax.servlet.http.HttpServletRequest;

import qspr.OCHEMConfiguration;
import qspr.dao.ChemInfEngine;
import qspr.entities.DataHandlingOptions;
import qspr.entities.ModelProtocol;
import qspr.metaserver.configurations.ModelAbstractConfiguration.ScalingType;
import qspr.workflow.utils.QSPRConstants;
import qspr.metaserver.configurations.StandartizationOptions;
import qspr.metaserver.configurations.ValidationConfiguration.MixtureValidation;

public class DataPreprocessingParser {
	
	public static void parseStandartizationUI(HttpServletRequest request, StandartizationOptions standartization, DataHandlingOptions dataHandling, ModelProtocol protocol)
	{
		if (standartization != null)
		{
			ChemInfEngine standardizer;
			if (request.getParameter("standardizer") != null) {
				String req = request.getParameter("standardizer");
				standardizer = ChemInfEngine.valueOf(req);
			} else {
				standardizer = OCHEMConfiguration.getCheminfEngine();
			}
			
			if (request.getParameter("desalt") != null)
				standartization.desaltWith = standardizer;
			else
			{
				standartization.desaltWith = null;
				if (protocol != null && 
						protocol.validationConfiguration != null && !protocol.validationConfiguration.isRecordValidated()) { // this protocol taken from previous step
					protocol.validationConfiguration.mixtureValidation = MixtureValidation.valueOf(request.getParameter("mixtureValidation").toUpperCase());
					if(protocol.validationConfiguration.mixtureValidation == MixtureValidation.MIXTURE) protocol.validationConfiguration.mixtureValidation = null;
				}
			}

			if (request.getParameter("neutralize") != null)
				standartization.neutralizeWith = standardizer;
			else
				standartization.neutralizeWith = null;

			if (request.getParameter("standardize") != null)
				standartization.standardizeWith = standardizer;
			else
				standartization.standardizeWith = null;

			if (request.getParameter("cleanstructure") != null)
				standartization.cleanStructureWith = standardizer;
			else
				standartization.cleanStructureWith = null;

			if (dataHandling != null)
			{
				if (request.getParameter("intervals") != null)
					dataHandling.intervals = QSPRConstants.USE;
				else
					dataHandling.intervals = null;

				if (request.getParameter("greaterless") != null)
					dataHandling.greaterless = QSPRConstants.USE;
				else
					dataHandling.greaterless = null;

				if (request.getParameter("approximateequals") != null)
					dataHandling.approximateequals = QSPRConstants.USE;
				else
					dataHandling.approximateequals = null;

				dataHandling.handleRanges = "true".equalsIgnoreCase(request.getParameter("handleRanges"));
			}

		}


	}

	/**
	 * returns array of two elements: x-scaling and y-scaling
	 */
	public static ScalingType[] parseNormalisationUI(HttpServletRequest request)
	{
		ScalingType xScaling = null;
		ScalingType yScaling = null;

		if (request.getParameter("normx") != null) {
			String normX = request.getParameter("normx");
			if (!normX.isEmpty())
				xScaling = ScalingType.valueOf(normX);
		}

		if (request.getParameter("normy") != null) {
			String normY = request.getParameter("normy");
			if (!normY.isEmpty())
				yScaling = ScalingType.valueOf(normY);
		}

		return new ScalingType[]{xScaling, yScaling};
	}
}
