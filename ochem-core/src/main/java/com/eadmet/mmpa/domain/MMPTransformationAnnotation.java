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

package com.eadmet.mmpa.domain;

import javax.persistence.Column;
import javax.persistence.Embedded;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import qspr.ThreadScope;
import qspr.entities.Basket;
import qspr.entities.Model;
import qspr.entities.Property;

@Entity
@XmlRootElement(name = "mmp-annotated-transformation")
public class MMPTransformationAnnotation
{
	
	@Id
	@GeneratedValue
	@Column(name = "tp_id")
	public Long id;
	
	@ManyToOne
	@JoinColumn(name = "transformation_id")
	@XmlTransient
	public MMPTransformation transformation;
	
	@ManyToOne
	@JoinColumn(name = "property_id")
	public Property property;
	
	@ManyToOne
	@JoinColumn(name = "as_id")
	public MMPAnnotationSet annotationSet;
	
	@Embedded
	public MMPSubsetStatistics subsetStats;
	
	@ManyToOne
	@JoinColumn(name = "basket_id")
	public Basket basket;
	
	@ManyToOne
	@JoinColumn(name = "model_id")
	public Model model;
	
	@XmlElement
	public boolean isIncreasing() {
		if (property.isQualitative())
			return subsetStats.nNP > subsetStats.nPN;
		else
			return subsetStats.deltaMean > 0;
	}

	@XmlElement
	protected MMPTransformation getTransformation() {
		if (ThreadScope.get().requestURI.contains("annotatedTransformationsList"))
			return transformation;
		else
			return null;
	}
}
