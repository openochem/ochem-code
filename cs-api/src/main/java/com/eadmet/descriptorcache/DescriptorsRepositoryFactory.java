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

package com.eadmet.descriptorcache;

import java.lang.reflect.Proxy;
import java.util.ArrayList;
import java.util.List;

import qspr.metaserver.transport.NoSqlTransport;
import qspr.metaserver.transport.ReattemptingProxy;


/**
 * A factory class to create instances for both clean and double-proxy-wrapped repository implementations
 * @author novserj
 *
 */
public class DescriptorsRepositoryFactory 
{
	public static DescriptorsRepository getRepository() throws Exception
	{
		DescriptorsRepositoryImpl r = new DescriptorsRepositoryImpl();
		DescriptorsRepository connection = (DescriptorsRepository)Proxy.newProxyInstance(r.getClass().getClassLoader(), new Class[]{DescriptorsRepository.class}, new ConnectionInjectingProxy(r));
		connection.initCollections();
		return connection;
	}

	public static DescriptorsRepository getReattemptingRepository() throws Exception
	{
		DescriptorsRepositoryImpl r = new DescriptorsRepositoryImpl();
		DescriptorsRepository connection = (DescriptorsRepository)Proxy.newProxyInstance(r.getClass().getClassLoader(), new Class[]{DescriptorsRepository.class}, new ConnectionInjectingProxy(r));
		DescriptorsRepository reattempt = (DescriptorsRepository)Proxy.newProxyInstance(r.getClass().getClassLoader(), new Class[]{DescriptorsRepository.class}, new ReattemptingProxy(connection));
		DescriptorsRepository portion = (DescriptorsRepository)Proxy.newProxyInstance(r.getClass().getClassLoader(), new Class[]{DescriptorsRepository.class}, new PortionRequestProxy(reattempt));

		portion.initCollections();
		return portion;
	}


	public static void main(String[] args) throws Exception
	{
		NoSqlTransport.host = "saria";
		int counter=1;
		//		for (int i=0; i<10; i++)
		//		{
		//			final int counter = i;
		//			Thread t = new Thread()
		//			{
		//				public void run()
		//				{
		//					try {
		DescriptorsRepository r = DescriptorsRepositoryFactory.getReattemptingRepository();
		for (int j=0; j<10; j++)
		{
			long timer = System.nanoTime();
			System.out.println("Starting query in thread "+counter+" iteration "+j);							
			List<DescriptorConfigEntry> list = r.getConfigurationsByType("ALogPS");
			List<CacheEntry> clist = r.getDescriptors(list.get(0));
			List<String> tmp = new ArrayList<String>();
			for (CacheEntry t : clist) 
				tmp.add(t.moleculeMD5);
			List<CacheEntry> nclist = r.getDescriptors(tmp.toArray(new String[0]), list.get(0));
			System.out.println("Finished query in thread "+counter+" iteration "+j+" in "+(System.nanoTime() - timer) / 1000000+"ms, got "+list.size()+" configs, "+clist.size()+" entries, "+nclist.size()+" entries");
		}
		//					} catch (Exception e) {
		//						e.printStackTrace();
		//					}
		//				}
		//			};
		//			t.start();
		//		}
	}
}