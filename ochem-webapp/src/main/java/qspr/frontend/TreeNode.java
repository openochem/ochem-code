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

package qspr.frontend;

import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.annotation.XmlAnyElement;
import javax.xml.bind.annotation.XmlElementWrapper;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

/**
 * A marshallable tree node for representing a "table of contents"
 * 
 * @author midnighter
 *
 */
@XmlRootElement(name = "tree-node")
public class TreeNode
{
	/**
	 * The object associated with this node
	 */
	@XmlAnyElement
	public Object object;
	
	/**
	 * Children nodes
	 */
	@XmlElementWrapper(name = "children")
	@XmlAnyElement
	public List<TreeNode> children = new ArrayList<TreeNode>();
	
	/**
	 * Parent nodes
	 */
	@XmlTransient
	public List<TreeNode> parents = new ArrayList<TreeNode>();
}
