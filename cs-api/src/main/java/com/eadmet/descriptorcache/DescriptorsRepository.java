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

import java.util.List;

import qspr.metaserver.transport.Reattempt;

/**
 * The abstraction for repository layer for accessing the descriptors storage.
 * Now, the only implementation is MongoDB, but it can be transparently moved to other type of storage (e.g., SQL)
 * 
 * @author midnighter, novserj
 *
 */
public interface DescriptorsRepository 
{
	/**
	 * Get the cache entries given a set of filters.
	 * 
	 * @param md5 - the "md5" field of the cache entry
	 * @param config - the descriptors configuration entity (might be a saved one or a new one)
	 * @param user - the user requesting the descriptors (can be NULL for public cache)
	 * @param shared - do we request only own entries (false) or entries shared by someone else and accessible to the user (true)
	 * @return - the list of cache entries (can contain NULLs)
	 * @throws Exception
	 */
	@Reattempt(name="descriptor repository")
	@Portion(size=5000)
	public List<CacheEntry> getDescriptors(String[] md5, DescriptorConfigEntry config) throws Exception;

	@Reattempt(name="descriptor repository")
	@Portion(size=5000)
	public List<CacheEntry> getDescriptors(Integer[] mp2, DescriptorConfigEntry config) throws Exception;

	@Reattempt(name="descriptor repository")
	public List<CacheEntry> getDescriptors(DescriptorConfigEntry config) throws Exception;

	@Reattempt(name="descriptor repository")
	public void saveDescriptors(List<CacheEntry> cacheEntries) throws Exception;

	@Reattempt(name="descriptor repository")
	public void saveConfig(DescriptorConfigEntry config) throws Exception;

	@Reattempt(name="descriptor repository")
	public DescriptorConfigEntry getDescriptorConfig(String md5) throws Exception;

	@Reattempt(name="descriptor repository")
	public DescriptorConfigEntry getDescriptorConfigById(String objectId) throws Exception;

	@Reattempt(name="descriptor repository")
	public List<DescriptorConfigEntry> getConfigurations(String user) throws Exception;

	@Reattempt(name="descriptor repository")
	public List<DescriptorConfigEntry> getAllConfigurations() throws Exception;

	@Reattempt(name="descriptor repository")
	public List<DescriptorConfigEntry> getConfigurationsByType(String descType) throws Exception;

	@Reattempt(name="descriptor repository")
	public void clearCache(String descType, String user) throws Exception;

	@Reattempt(name="descriptor repository")
	public void clearCache(String descType) throws Exception;

	@Reattempt(name="descriptor repository")
	public void clearCache(DescriptorConfigEntry descType) throws Exception;

	@Reattempt(name="descriptor repository")
	public void initCollections() throws Exception;

	@Reattempt(name="descriptor repository")
	public void clearLostEntries() throws Exception;

	@Reattempt(name="descriptor repository")
	public void updateDescriptors(DescriptorConfigEntry oldConfig,  DescriptorConfigEntry newConfig) throws Exception;

}