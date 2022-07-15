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

import java.io.BufferedWriter;
import java.io.StringWriter;

import org.junit.Test;

import org.junit.Assert;

import com.eadmet.utils.config.GlobalConfigurator;
import com.eadmet.utils.mailer.Mailer;

/**
 * @author midnighter
 */
public class GlobalConfiguratorTest
{
	// A simple test based on an external file
	@Test
	public void basicTest() throws Exception
	{
		GlobalConfigurator configurator = new GlobalConfigurator();
		configurator.addResource(GlobalConfiguratorTest.class.getClassLoader().getResource("test-config.cfg"));
		configurator.configure();
		
		Assert.assertEquals(Mailer.mailHost, "7.7.7.7");
		Assert.assertEquals(CommonConfig.authorizedIPs, "127.0.0.1");
	}
	
	@Test
	public void exportTest() throws Exception
	{
		GlobalConfigurator configurator = new GlobalConfigurator();
		configurator.addResource(GlobalConfiguratorTest.class.getClassLoader().getResource("test-config.cfg"));
		configurator.configure();
		
		StringWriter sw = new StringWriter();
		BufferedWriter bw = new BufferedWriter(sw);
		GlobalConfigurator.export(bw);
		
		String conf = sw.toString();
		System.out.println(conf);
		Assert.assertEquals(true, conf.contains("mailer.mail_host = 7.7.7.7"));
	}
}
