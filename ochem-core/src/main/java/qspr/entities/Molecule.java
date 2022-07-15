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

package qspr.entities;

import java.io.IOException;
import java.net.URLEncoder;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.ThreadScope;
import qspr.dao.Repository;
import qspr.dao.Various;
import qspr.depiction.MoleculeDepiction;
import qspr.workflow.utils.SDFProcessor;

import com.eadmet.utils.OCHEMUtils;


@XmlRootElement(name = "molecule")
public class Molecule
{
	@XmlAttribute
	public Long id;

	@XmlTransient
	public String md5;

	@XmlTransient
	public String pictureMd5;

	@XmlTransient
	private byte[] data;

	@XmlAttribute
	public String methodRestrict;

	@XmlTransient
	public Molecule persistentMolecule;

	@XmlTransient
	public byte[] molImage;

	@XmlTransient
	public Mapping1 mapping1;

	@XmlTransient
	public Mapping2 mapping2;

	@XmlElement(required = true)
	public List<Alert> methodError;

	@XmlElement(name="mapping")
	public Mapping2 getMapping()
	{
		if (!ThreadScope.get().controller.equals("molbrowser"))
			return mapping2;
		else
			return null;
	}

	@XmlElement
	public List<String> searchSynonymList;

	@XmlAttribute
	public String error;

	@XmlAttribute
	public String smile;

	@XmlAttribute
	public Long searchExistingMolId;

	@XmlAttribute
	public String searchCompSign;

	@XmlElement
	public MoleculeName searchedBy;

	public Double molWeight;

	@XmlTransient
	public byte[] getImage() throws IOException 
	{
		if(molImage == null || pictureMd5 == null) {
			if (molWeight == 0) {
				Molecule molecule = Repository.molecule.getEmptyMolecule(); // If we can't calculate even weight; this is erroneous molecule
				molImage = molecule.molImage;
				pictureMd5 =molecule.pictureMd5;
			}
			else
			{
				MoleculeDepiction depiction = MoleculeDepiction.get(OCHEMConfiguration.getCheminfEngine());
				depiction.setMolecule(this.getData());
				if(this.mapping1.inchi1.length() != 14) {
		    		depiction.setError(true);
		    	}
				depiction.setWidth(150);
				depiction.setHeight(150);
				if (error != null) {
					depiction.setColor("F5A9A9");
				}
				molImage = depiction.getImage();
				pictureMd5 = OCHEMUtils.getMD5(molImage);
			}
			Globals.session().saveOrUpdate(this);
		}

		return molImage;
	}


	public void setImage(byte[] image) 
	{
		molImage = image;
	}

	public double getMolWeight() throws IOException
	{
		if (molWeight != null)
			return molWeight;

		if (!isEmptyMolecule() && data != null)
		{
			return molWeight = Various.molecule.getMass(getData());
		}

		return 0.0;
	}

	@XmlTransient
	public byte[] getRawData() {
		return data;
	}


	@XmlTransient
	public String getData()
	{
		if (data == null)
			return null;

		try 
		{
			return new String(OCHEMUtils.MySqlCompatibleUncompress(data));
		} catch (IOException e)
		{
			throw new RuntimeException(e);
		}
	}

	public boolean hasData() {
		return data != null;
	}

	@XmlAttribute(name = "mp2")
	private Integer getMp2()
	{
		return mapping2.id;
	}

	@XmlElement
	public String getEncodedData()
	{
		if (ThreadScope.get().controller.equals("molecule") || ThreadScope.get().controller.equals("ncbisearch"))
		{
			try 
			{
				String ldata = new String(OCHEMUtils.MySqlCompatibleUncompress(data));

				return URLEncoder.encode(ldata,"utf8").replaceAll("\\+", "%20");
			} catch (Exception e)
			{
				throw new RuntimeException(e);
			}
		}
		else
			return null;
	}

	public void setData(String _data)
	{
		try 
		{
			data = _data == null || "NULL".equals(_data.toUpperCase())? null: OCHEMUtils.MySqlCompatibleCompress(_data.replace("$$$$", ""));

			if(data == null)
				throw new Exception("Molecule data are null"); //This is not normal!!!

		} catch (Exception e)
		{
			e.printStackTrace();
			throw new RuntimeException(e);
		}
	}

	@XmlElement(name = "property")
	public Set<CalculatedProperty> calculatedProperties;

	@XmlTransient
	public Set<ExperimentalProperty> experimentalProperties = new HashSet<ExperimentalProperty>();

	@XmlTransient
	public Hashtable<String, Double> fragments = new Hashtable<String, Double>();

	// For loading from SDF file
	@XmlTransient
	public Hashtable<String, String> miscProperties = new Hashtable<String, String>();

	public CalculatedProperty getCalculatedProperty(Property property)
	{
		CalculatedProperty result;

		if (calculatedProperties == null)
			return null;

		Iterator<CalculatedProperty> iCP = calculatedProperties.iterator();
		while (iCP.hasNext())
		{
			result = iCP.next();
			if (result.property == property)
				return result;
		}

		return null;
	}

	public CalculatedProperty addCalculatedProperty(Property property, double value)
	{
		if (calculatedProperties == null)
			calculatedProperties = new HashSet<CalculatedProperty>();
		CalculatedProperty cp = new CalculatedProperty();
		cp.molecule = this;
		cp.property = property;
		cp.value = value;

		calculatedProperties.add(cp);
		return cp;
	}

	public Molecule getPersistentMolecule()
	{
		if (persistentMolecule != null)
			return persistentMolecule;
		else if (id != null)
			return this;
		else
			return null;
	}


	public Molecule()
	{
		calculatedProperties = new HashSet<CalculatedProperty>();
	}

	public static Molecule getStub() {
		Molecule mol = new Molecule();
		mol.mapping2 = new Mapping2();
		mol.mapping2.id = -1;
		mol.mapping2.inchi2 = "none";
		mol.mapping1 = new Mapping1();
		mol.mapping1.id = -1L;
		mol.id = -1L;

		return mol;
	}

	public static Molecule getStub(long depictionId) {
		Molecule mol = Molecule.getStub();
		mol.id = depictionId;
		return mol;
	}

	public static String[] getInChiKeys(String mol){
		String inchiKey = Various.molecule.getInChiKey(mol);

		String[] keys = inchiKey.split("-");

		if (keys.length >= 2) {
			String k[] = {keys[0],keys[1]};
			return k;
		}
		logger.warn("Bad inchi key: " + inchiKey);

		String[] md5 = new String[2];

		md5[0]=md5[1]=OCHEMUtils.getMD5(SDFProcessor.standartize(mol));

		return md5;
	}

	public void updateInchi()
	{
		updateInchi(null);
	}

	/*
	 * if this.getInchi() was called earlier in the code, it doesn't have to be calculated again
	 * in updateInchi.
	 */
	public void updateInchi(String[] inchi)
	{
		try
		{
			String ldata = new String(OCHEMUtils.MySqlCompatibleUncompress(data));
			if (inchi == null) 
				inchi = getInChiKeys(ldata);

			mapping2 = Repository.molecule.getMapping2(inchi[0], inchi[1], ldata);
			mapping1 = mapping2.mapping1;

			if(OCHEMConfiguration.verboseMode > 1)logger.info("New mapping is "+mapping1.id+" "+mapping2.id);

		} catch (Exception e)
		{
			e.printStackTrace();
			logger.error("Something really bad happened with inchies.");
		}

	}

	public boolean isEmptyMolecule()
	{
		return mapping2.isEmpty();
	}

	public String toString()
	{
		return "MolID " + (mapping2 != null ? " M" + mapping2.id : "") + " (Depiction " + id + ")";
	}

	private static Logger logger = LogManager.getLogger(Molecule.class); 


	public String getImageFormat(){
		return mapping1.inchi1.length() == 14? "png:w150,h150,#ffffff" : "png:w150,h150,#F5A9A9";
	}

}
