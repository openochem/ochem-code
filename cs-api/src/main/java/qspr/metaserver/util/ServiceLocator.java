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

import javax.xml.rpc.ServiceException;

import ochem.eadmet.wsapi.ModelServiceLocator;
import ochem.eadmet.wsapi.ModelServicePortType;

import org.apache.axis.AxisProperties;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import qspr.metaserver.serv.ServerRunnerConfiguration;

/**
 * The locator of the OCHEM web service.
 * Uses the global configuration parameter to identify OCHEM URL
 * 
 * @author midnighter
 *
 */
public class ServiceLocator
{
	private static final Logger logger = LoggerFactory.getLogger(ServiceLocator.class);
	
	public static ModelServicePortType getModelService() throws ServiceException
	{
		// Ignore invalid SSL certificates. 
		AxisProperties.setProperty("axis.socketSecureFactory","org.apache.axis.components.net.SunFakeTrustSocketFactory");
		
		ModelServiceLocator locator = new ModelServiceLocator();
		if (ServerRunnerConfiguration.instance != null && ServerRunnerConfiguration.instance.ochemURL != null)
		{
			locator.setModelServiceHttpSoap11EndpointEndpointAddress(ServerRunnerConfiguration.instance.ochemURL + "/services/ModelService.ModelServiceHttpSoap11Endpoint");
			locator.setModelServiceHttpSoap12EndpointEndpointAddress(ServerRunnerConfiguration.instance.ochemURL + "/services/ModelService.ModelServiceHttpSoap12Endpoint");
			logger.info("Using web service URL " + locator.getModelServiceHttpSoap11EndpointAddress());
		}
		
		return locator.getModelServiceHttpSoap11Endpoint();
	}
}
