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

import com.eadmet.exceptions.UserFriendlyException;

import qspr.metaserver.protocol.Task;

@Deprecated
@XmlRootElement(name="dragon-configuration")
public class DescriptorsDragonConfiguration extends DescriptorsAbstractDragonConfiguration
{
	private static final long serialVersionUID = 1L;

	public DescriptorsDragonConfiguration()
	{
	}

	public DescriptorsDragonConfiguration(int dragonBlocks)
	{
		super(dragonBlocks);
		throw new UserFriendlyException("This Configuration of Dragon is obsolete and should not be used. Something is wrong. Contact " + Task.EMAIL_OCHEM);
	}

	@Override
	public int getBlocks() {
		return 0;
		//throw new UserFriendlyException("This Configuration should not be used. Something is wrong. Contact " + Task.EMAIL_OCHEM);
	}


	@Override
	public String getExtension() {
		throw new UserFriendlyException("This Configuration of Dragon is obsolete and should not be used. Something is wrong. Contact " + Task.EMAIL_OCHEM);
	}

	@Override
	protected long get3DMask() {
		return 0x4FF80;
	}

}
