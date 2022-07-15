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

import java.util.HashMap;
import java.util.Map;

import javax.persistence.CascadeType;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.FetchType;
import javax.persistence.GeneratedValue;
import javax.persistence.GenerationType;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.Transient;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.criterion.Disjunction;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.MarshallingOption;
import qspr.interfaces.Descriptable;
import qspr.modelling.ModelStatistics;

@Entity
@XmlRootElement(name = "modelMapping")
public class ModelMapping 
{
	private static transient final Logger logger = LogManager.getLogger(ModelMapping.class);

	@Id
	@Column(name = "model_mapping_id")
	@GeneratedValue(strategy = GenerationType.AUTO)
	@XmlAttribute
	public Long id;

	@Column(name = "class")
	@XmlAttribute
	public Long _class;


	@ManyToOne
	@JoinColumn(name = "model_id")
	@XmlTransient
	public Model model;


	@ManyToOne
	@JoinColumn(name = "property_id")
	@XmlElement
	public Property property;

	@ManyToOne
	@JoinColumn(name = "unit_id")
	@XmlElement
	public Unit unit;


	@ManyToOne(cascade={CascadeType.ALL}, fetch = FetchType.LAZY)
	@JoinColumn(name = "stat_original")
	@XmlTransient
	public Attachment<ModelStatistics> statisticsOriginal;

	@ManyToOne(cascade={CascadeType.ALL}, fetch = FetchType.LAZY)
	@JoinColumn(name = "stat_recalculated")
	@XmlTransient
	public Attachment<ModelStatistics> statisticsRecalculated; 

	public void updateDescription()
	{
		Map<String, Object> parameters = new HashMap<String, Object>();
		parameters.putAll(((Descriptable)statisticsRecalculated.getObject()).getParameters());
		logger.info("Property - "+this.property+"\nstatistics - "+parameters);
	}

	@XmlElement(name = "model")
	public Model getModel()
	{
		Model copymodel = new Model();
		copymodel.id = model.id;
		copymodel.template = model.template;
		copymodel.name = model.name;
		copymodel.description = model.description;
		return copymodel;
	}

	@XmlAttribute(name = "trainingSize")
	protected Long getTrainingSize()
	{
		if (!Globals.getMarshallingOption(MarshallingOption.MODEL_SET_SIZE))
			return null;
		Criteria criteria = Globals.session().createCriteria(BasketEntry.class);

		Disjunction excluded = Restrictions.disjunction();

		if (model.microattachment.getObject().excludedBasketEntries.size() > 0)
			excluded.add(Restrictions.in("id", model.microattachment.getObject().excludedBasketEntries));
		else
			excluded.add(Restrictions.sqlRestriction("1 = 0"));

		excluded.add(Restrictions.eq("exclude", Boolean.TRUE));

		criteria.add(Restrictions.not(excluded));

		criteria.add(Restrictions.eq("basket", this.model.trainingSet));
		criteria.createAlias("ep", "e");
		criteria.createAlias("e.property", "p");

		if (this.property.isDirectory)
		{

			criteria.add(Restrictions.eq("p.parent", this.property));
		}
		else
			criteria.add(Restrictions.eq("p.id", this.property.id));

		criteria.setProjection(Projections.rowCount());
		return (Long)criteria.list().get(0);
	}

	@XmlAttribute(name = "validationSize")
	protected Long getValidationSize()
	{
		if(this.model.getValidationSets().size() == 0)
			return 0l;
		
		if (!Globals.getMarshallingOption(MarshallingOption.MODEL_SET_SIZE))
			return null;
		Criteria criteria = Globals.session().createCriteria(BasketEntry.class);

		//
		Disjunction excluded = Restrictions.disjunction();

		if (model.microattachment.getObject().excludedBasketEntries.size() > 0)
			excluded.add(Restrictions.in("id", model.microattachment.getObject().excludedBasketEntries));
		else
			excluded.add(Restrictions.sqlRestriction("1 = 0"));

		excluded.add(Restrictions.eq("exclude", Boolean.TRUE));

		criteria.add(Restrictions.not(excluded));

		criteria.add(Restrictions.eq("basket", this.model.getValidationSets().get(0)));
		criteria.createAlias("ep", "e");
		criteria.createAlias("e.property", "p");
		if (this.property.isDirectory)
		{

			criteria.add(Restrictions.eq("p.parent", this.property));
		}
		else
			criteria.add(Restrictions.eq("p.id", this.property.id));

		criteria.setProjection(Projections.rowCount());
		return (Long)criteria.list().get(0);
	}

	@XmlAttribute(name = "excludedSize")
	protected Long getExcludedSize()
	{
		if (!Globals.getMarshallingOption(MarshallingOption.MODEL_SET_SIZE))
			return null;

		Criteria criteria = Globals.session().createCriteria(BasketEntry.class);

		Disjunction excluded = Restrictions.disjunction();

		if (model.microattachment.getObject().excludedBasketEntries.size() > 0)
			excluded.add(Restrictions.in("id", model.microattachment.getObject().excludedBasketEntries));
		else
			excluded.add(Restrictions.sqlRestriction("1 = 0"));

		excluded.add(Restrictions.eq("exclude", Boolean.TRUE));

		criteria.add(excluded);
		criteria.add(Restrictions.eq("basket", this.model.trainingSet));
		criteria.createAlias("ep", "e");
		criteria.createAlias("e.property", "p");
		if (this.property.isDirectory)
		{
			criteria.add(Restrictions.eq("p.parent", this.property));
		}
		else
			criteria.add(Restrictions.eq("p.id", this.property.id));

		criteria.setProjection(Projections.rowCount());
		return (Long)criteria.list().get(0);
	}

	@XmlTransient
	@Transient
	public Double classificationThreshlod;
	
	@XmlTransient
	public int getIndex()
	{
		return model.modelMappings.indexOf(this);
	}

	public boolean matches(ExperimentalProperty ep)
	{
		return property.equals(ep.property) || property.equals(ep.property.parent);
	}

	@Override
	public boolean equals(Object obj)
	{
		ModelMapping mm = (ModelMapping)obj; 
		if (this.id != null)
			return mm != null && this.id.equals(mm.id);
		else
			return this == obj;
	}

	@Override
	public int hashCode()
	{
		if (this.id != null)
			return this.id.hashCode();
		else
			return super.hashCode();
	}
}
