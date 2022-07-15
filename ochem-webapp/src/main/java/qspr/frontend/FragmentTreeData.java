
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
import java.util.TreeSet;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement
public class FragmentTreeData 
{
	//@XmlElementWrapper(name = "nodes")
	//@XmlElement(name = "node")
	public TreeSet<FTTreeNode> nodes = new TreeSet<FTTreeNode>();
	
	//@XmlElementWrapper(name = "edges")
	//@XmlElement(name = "edge")
	public List<TreeEdge> edges = new ArrayList<TreeEdge>();
	
	public boolean addNode(long x, long y, long fragId, String info)
	{
		FTTreeNode node = new FTTreeNode();
		node.x = x;
		node.y = y;
		node.fragId = fragId;
		node.info = info;		
		if (nodes.contains(node))
		{
			return false;
		} else {
			nodes.add(node);
			return true;
		}
	}
	
	public void addEdge(long fragFrom, long fragTo, String info)
	{
		TreeEdge edge = new TreeEdge();
		edge.fragFrom = fragFrom;
		edge.fragTo = fragTo;
		edge.info = info;
		edges.add(edge);
	}
	
	public FragmentTreeData()
	{
	}
}

class FTTreeNode implements java.lang.Comparable<FTTreeNode>
{
	@XmlAttribute
	long x;
	@XmlAttribute
	long y;
	@XmlAttribute
	long fragId;
	@XmlElement
	String info;
	
	public FTTreeNode()
	{		
	}
	
	public int compareTo(FTTreeNode n)
	{
		if (this.fragId == n.fragId) return 0;
		if (this.fragId < n.fragId) return -1; else return 1;
	}
		
}

class TreeEdge
{
	@XmlAttribute
	long fragFrom;
	@XmlAttribute
	long fragTo;
	@XmlElement
	String info;
	
	public TreeEdge()
	{
		
	}
}