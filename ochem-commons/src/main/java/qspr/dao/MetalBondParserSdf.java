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

package qspr.dao;

import java.util.HashSet;
import java.util.Set;

/**
 * 
 * @author itetko
 * Parse SDF file to eliminate Metalic Bonds
 */

public class MetalBondParserSdf {

	static Set<Integer> pt=null;

	public static String eliminateMetalBond(String sdf){
		return processMetals(sdf,HandleMetal.DELETEBOND);
	}

	public static String substituteMetalBondwithSingleBond(String sdf){
		return processMetals(sdf,HandleMetal.SUBSTITUTEBOND);
	}

	public static String eliminateMetal(String sdf){
		return processMetals(sdf,HandleMetal.DISCONNECTMETAL);
	}

	private static String replaceCharAt(String s, int pos, char c) {
		return s.substring(0,pos) + c + s.substring(pos+1);
	}

	private static int countMetalBonds(String sdf, HandleMetal mode){
		String lines[]=sdf.split("\n");
		int count=0;

		try{
			int atoms=Integer.valueOf(lines[3].substring(0,3).replace(" ","")), bonds=Integer.valueOf(lines[3].substring(4,8).replace(" ",""));
			atoms+=4;
			bonds+=atoms;

			if(mode==HandleMetal.DISCONNECTMETAL){
				pt=new HashSet<Integer>();
				for(int i=4;i<atoms;i++){
					String satoms[]=lines[i].split("\\s++");
					if(satoms[3].equals("Pt")||satoms[4].equals("Pt"))
						pt.add(i-3);
				}
			}

			for(int i=atoms;i<bonds;i++)if(containMetal(lines[i]))count++;
		}catch(Exception e){
			return 0;
		}
		return count;
	}

	static boolean containMetal(String line){
		if(line.charAt(8)=='8')return true;
		if(pt!=null)
			for(int j=0;j<2;j++){
				String mol=line.substring(j*3,(j+1)*3).replaceAll("\\s+","");
				int atom=Integer.valueOf(mol);
				if(pt.contains(atom))return true;
			}
		return false;
	}

	private static String processMetals(String sdf,HandleMetal mode){

		int count=countMetalBonds(sdf,mode);
		if(count==0)return sdf;

		String m="";
		String lines[]=sdf.split("\n");
		try{
			int atoms=Integer.valueOf(lines[3].substring(0,3).replace(" ","")), bonds=Integer.valueOf(lines[3].substring(4,8).replace(" ",""));
			atoms+=4;
			bonds+=atoms;
			for(int i=0;i<lines.length;i++){
				if(i==3&&mode!=HandleMetal.SUBSTITUTEBOND){
					String ss=lines[3].substring(0,3)+String.format("%3d",bonds-count-atoms)+lines[3].substring(6);
					m+=ss+"\n";
					continue;
				}

				if(i>=atoms&i<bonds&&containMetal(lines[i]))
					switch(mode){
					case DELETEBOND: continue;
					case DISCONNECTMETAL: continue; 
					case SUBSTITUTEBOND: lines[i]=replaceCharAt(lines[i],8,'1');break;
					}
				m+=lines[i];
				if(!lines[i].contains("M  END"))m+="\n";
			}
		}catch(Exception e){
			return sdf;
		}

		return m;
	}

	public static String processBond9(String sdf){

		String m="";
		String lines[]=sdf.split("\n");
		try{
			int atoms=Integer.valueOf(lines[3].substring(0,3).replace(" ","")), bonds=Integer.valueOf(lines[3].substring(4,8).replace(" ",""));
			atoms+=4;
			bonds+=atoms;
			for(int i=0;i<lines.length;i++){

				if(i>=atoms&i<bonds&&lines[i].charAt(8)=='9')
					lines[i]=replaceCharAt(lines[i],8,'8');
				m+=lines[i];
				if(!lines[i].contains("M  END"))m+="\n";
			}
		}catch(Exception e){
			return sdf;
		}

		return m;
	}

}

