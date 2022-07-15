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

// Currently only used in ECMolecule 

import java.io.Serializable;
import java.nio.ByteBuffer;
import java.sql.Timestamp;
import java.util.BitSet;
import java.util.Calendar;
import java.util.List;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.FetchType;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.Transient;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import com.eadmet.exceptions.UserFriendlyException;

import qspr.metaserver.configurations.ModelAbstractConfiguration.PredictionScenario;
import qspr.workflow.utils.QSPRConstants;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataRow;

@Entity
@XmlRootElement
public class CachedPrediction implements Serializable
{
	private static final long serialVersionUID = 1L;

	@Id
	@GeneratedValue
	@Column(name="cached_prediction_id")
	@XmlAttribute
	public Integer id;

	@Column(name = "mapping2_id") 
	@XmlAttribute
	public Integer mp2Id;

	@ManyToOne(fetch=FetchType.LAZY)
	@JoinColumn(name="model_id")
	@XmlTransient
	public Model model;

	@Transient
	private AbstractDataRow row;

	@Column(name="data")
	private byte data[]; 

	@Column(columnDefinition = "BINARY(2)", length = 2)
	private byte properties[];

	@Transient
	private BitSet  bs;

	/**
	 * Stores the exact molecule that was cached.
	 * this is useful for consistency testing
	 */
	@Column(name = "molecule_id")
	public Long molId;

	@Column(name = "time")
	@XmlTransient
	public Timestamp cachingTime;

	@Column(columnDefinition = "TINYINT")
	public byte cachetype;
	
	public final static int CACHESIZE = 16; 

	public CachedPrediction()
	{

	}

	public CachedPrediction(Model model, Molecule mol, List<String>columns, AbstractDataRow row, PredictionScenario type)
	{
		this.model = model;
		this.mp2Id = mol.mapping2.id;
		this.molId = mol.id;
		this.cachingTime = new Timestamp(Calendar.getInstance().getTimeInMillis());
		this.row = row;
		bs = new BitSet(CACHESIZE);
		float values[] = new float[columns.size()];
		int nonzero = 0;
		for(int i = 0 ; i < columns.size() ; i++){
			if(columns.get(i).contains(QSPRConstants.INDIVIDUAL_PREDICTIONS))  // so far we do not store Individual predictions; also that is why dtResult is not compact
				continue;
			values[i] = ((Double)row.getValue(i)).floatValue();
			if(values[i] != 0)nonzero = i; // zero elements are stored by default
			if(columns.get(i).startsWith(QSPRConstants.PREDICTION_RESULT_COLUMN)){
				if(i>=CACHESIZE)throw new UserFriendlyException("Cannot cache molecule: increase CACHESIZE="+CACHESIZE+" for columns "+columns.size());
				bs.set(i);
			}
		}
		nonzero++; // at least one element
		properties = bs.toByteArray();
		ByteBuffer buf = ByteBuffer.allocate(4*nonzero);
		for(int i = 0; i < nonzero ; i++)
			buf.putFloat(values[i]);
		this.data =  buf.array();
		this.cachetype = type.getValue();
	}

	public AbstractDataRow getDataRow(int columnSize)
	{		
		if(row == null){
			row = new DataRow();
			ByteBuffer b = ByteBuffer.wrap(data);
			int i;
			for(i=0; i<data.length/4; i++)
				row.setValue(i, (double) b.getFloat());
			for(;i < columnSize; i++)
				row.setValue(i, (double)0); // required, since 0 are not stored!
		}
		row.addAttachment(QSPRConstants.CACHED, true);
		return row;		
	}

	public Double getValue(int propNum)
	{
		if(bs == null)
			bs = BitSet.valueOf(properties);

		for (int i = bs.nextSetBit(0); i >= 0; i = bs.nextSetBit(i+1))
			if(propNum-- == 0)
				return (Double)getDataRow(i).getValue(i);
		throw new UserFriendlyException("Could not find cached property");
	}

}
