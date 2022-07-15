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

import java.io.Writer;

import javax.servlet.http.HttpServletRequest;
import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement(name="silicosit-configuration")
public class DescriptorsSilicosItScaffoldConfiguration extends DescriptorsAbstractConfiguration{

	private static final long serialVersionUID = 1L;

	public Boolean oprea;
	public Boolean schuffenhauer;

	public Boolean scaffolds; 
	public Boolean frameworks;

	public DescriptorsSilicosItScaffoldConfiguration() {
		scaffolds=true;
	}

	public void writeScaffoldParameters(Writer bw) throws Exception {
		if (scaffolds!=null && scaffolds) {
			bw.write("RINGS_WITH_LINKERS_2\n");
			bw.write("RINGS_WITH_LINKERS_1\n");
		}
		if (frameworks!=null && frameworks) {
			bw.write("MURCKO_2\n");
			bw.write("MURCKO_1\n");
		}
		if (oprea!=null && oprea) {
			bw.write("OPREA_3\n");
			bw.write("OPREA_2\n");
			bw.write("OPREA_1\n");
		}
		if (schuffenhauer!=null && schuffenhauer) {
			bw.write("SCHUFFENHAUER_5\n");
			bw.write("SCHUFFENHAUER_4\n");
			bw.write("SCHUFFENHAUER_3\n");
			bw.write("SCHUFFENHAUER_2\n");
			bw.write("SCHUFFENHAUER_1\n");
		}
	}


	@Override
	public boolean requires3D() {
		return false;
	}	

	// Note:
	// Descriptor blocks Silicos-It scaffolds:
	// schuffenhauer is for Schuffenhauer scaffolds (hierarchical removing of rings to last 1-5 rings in system)
	// oprea is for Oprea frameworks (keep only topology of rings, 2,3 adding hydrogen donors and acceptors)
	// framworks is for Bemis-Murcko frameworks (2 contracts linkers to unity length)
	// scaffolds is for Bemis-Murcko scaffolds (2 keeps double bonded substituents on rings and linkers)

	@Override
	public String toString(){
		return super.toString() + (oprea != null && oprea ? "oprea" :"") + (schuffenhauer != null && schuffenhauer? "schuffenhauer":"");

	}

	@Override
	DescriptorsSilicosItScaffoldConfiguration setConfiguration(HttpServletRequest request){
		scaffolds = request.getParameter("silicos_scaffolds") != null;
		frameworks = request.getParameter("silicos_frameworks") != null;
		oprea = request.getParameter("silicos_oprea") != null;
		schuffenhauer = request.getParameter("silicos_schuffenhauer") != null;
		return this;
	}

	@Override
	public String getDefaultTypeName() {
		return DescriptorsConfiguration.SilicosItScaffold;
	}
}
