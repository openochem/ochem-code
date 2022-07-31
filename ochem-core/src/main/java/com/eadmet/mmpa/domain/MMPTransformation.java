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

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.persistence.CascadeType;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.FetchType;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.OneToMany;
import javax.persistence.Transient;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;
import org.hibernate.sql.JoinType;

import qspr.Globals;
import qspr.MarshallingOption;
import qspr.dao.Various;
import qspr.entities.Mapping2;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.utils.MemoryUtils;

@Entity
@XmlRootElement(name = "mmp-transformation")
@SuppressWarnings("unchecked")
public class MMPTransformation
{
	private static final Logger logger = LogManager.getLogger(MMPTransformation.class);

	@Id
	@Column(name = "transformation_id")
	@GeneratedValue
	public Long id;

	@Column(name = "fragment1_id")
	@XmlTransient
	private long frag1Id;

	@Column(name = "fragment2_id")
	@XmlTransient
	private long frag2Id;

	@Column(name = "pairs_count")
	@XmlTransient
	public long pairsCount;

	@Transient
	public MMPSubsetStatistics statistics;

	/**
	 * An indication that this transformation should be inversed 
	 * (for display and molecule optimisation purposes). 
	 * This flag is generated on the fly depending on the transformation filtering option
	 */
	@XmlAttribute
	@Transient
	public boolean inversed;

	private static Map<String, Long> cache = new HashMap<String, Long>();

	/**
	 * Filtered MMPs for this transformation.
	 * For example, MMPs from a particular basket.
	 */
	public transient List<MMPair> filteredPairs;

	public static MMPTransformation get(long id) {
		return (MMPTransformation) Globals.session().get(MMPTransformation.class, id);
	}

	@XmlElement
	public long getFrag1Id() {
		return inversed ? frag2Id : frag1Id;
	}

	@XmlElement
	public long getFrag2Id() {
		return inversed ? frag1Id : frag2Id;
	}

	public static MMPair getMMPair(long frag1[], Mapping2 mol[]) {
		long transformation = getId(frag1[0],frag1[1]);
		MMPair pair = new MMPair();

		if(transformation > 0){
			pair.molecule1 = mol[0];
			pair.molecule2 = mol[1];
		}else{
			pair.molecule1 = mol[1];
			pair.molecule2 = mol[0];
		}

		pair.transformation = new MMPTransformation();
		pair.transformation.id=Math.abs(transformation);
		return pair;
	}

	private static long getId(long frag1, long frag2) {

		if (cache.isEmpty())
		{
			logger.info("Filling the MMP cache");
			logger.info(MemoryUtils.memorySummary());

			List<Object[]> rows = Globals.session().createCriteria(MMPTransformation.class).setProjection(Projections.projectionList().add(Projections.id()).add(Projections.property("frag1Id")).add(Projections.property("frag2Id"))).list();
			for (Object[] row : rows)
				putCache((Long)row[1],(Long)row[2],(Long)row[0]);

			logger.info("MMP cache filled with " + cache.size() + " elements");
			logger.info(MemoryUtils.memorySummary());
		}

		String cacheString = getCacheString(frag1,frag2);
		Long transformationId = cache.get(cacheString);
		if(transformationId != null) return frag1<frag2?transformationId:-transformationId; // we need in the reversed or direct order

		// now long operation: comparison is first by size and after that according to the IDs

		long id1 = selectMinimalId(frag1, frag2);
		long id2 = id1 == frag1 ? frag2 : frag1;

		MMPTransformation transformation = new MMPTransformation();
		transformation.frag1Id = id1;
		transformation.frag2Id = id2;
		Globals.session().save(transformation);
		putCache(id1,id2,transformation.id);

		return id1 == frag1? transformation.id : -transformation.id;
	}

	static void putCache(Long frag1, Long frag2, Long trans){
		cache.put(getCacheString(frag1,frag2),frag1 < frag2? trans : -trans);
	}

	static String getCacheString(Long frag1, Long frag2){
		return frag1<frag2?Long.toBinaryString(frag1)+Long.toBinaryString(frag2): Long.toBinaryString(frag2)+Long.toBinaryString(frag1);
	}
	
	static double getMass(String smiles) throws Exception {
		try {
			return Various.molecule.getMass(Various.molecule.convertToFormat(smiles, QSPRConstants.SDFNOAROM_NOH));
		} catch (IOException e) {
			e.printStackTrace();
			throw new Exception("Failed to calculate mass.");
		}
	}

	static long selectMinimalId(long frag1, long frag2){
		MMPFragment fra1 = (MMPFragment) Globals.session().get(MMPFragment.class, frag1);
		MMPFragment fra2 = (MMPFragment) Globals.session().get(MMPFragment.class, frag2);

		if(fra1.size != fra2.size)
			return fra1.size < fra2.size ? fra1.id : fra2.id;

		double mw1, mw2;

		try{
			mw1 = getMass(fra1.smiles);
			mw2 = getMass(fra2.smiles);
		}catch(Exception e){
			logger.info("1st attempt failed to calculate MW for " +fra1.smiles + " "+fra2.smiles);
			mw1 = mw2 = 0;
		}

		if(mw1 == mw2  && mw1 == 0)	try{ // second attempt using H
			mw1 = getMass(fra1.smiles.replaceAll("\\[Al\\]", "H"));
			mw2 = getMass(fra2.smiles.replaceAll("\\[Al\\]", "H"));
		}catch(Exception e){
			logger.info("2nd attempt failed to calculate MW for " +fra1.smiles + " "+fra2.smiles);
			mw1 = mw2 = 0;
		}

		if(mw1 == mw2  && mw1 == 0)	try{ // third attempt using nothing
			mw1 = getMass(fra1.smiles.replaceAll("\\[Al\\]", ""));
			mw2 = getMass(fra2.smiles.replaceAll("\\[Al\\]", ""));
		}catch(Exception e){
			logger.info("3rd attempt failed to calculate MW for " +fra1.smiles + " "+fra2.smiles);
			mw1 = mw2 = 0;
		}

		if(fra1.size == 1) // one element is always [Al]
			return mw1 < mw2 ? fra1.id : fra2.id;

		// same size
		if(fra1.carbonChain != fra2.carbonChain){ // carbon to non carbon transformation is going first
			if(fra1.carbonChain) return fra1.id;
			return fra2.id;
		}

		if(mw1 != mw2)
			return mw1 < mw2 ? fra1.id : fra2.id;

		if(fra1.smiles.length() != fra2.smiles.length())
			return fra1.smiles.length()  < fra2.smiles.length() ? fra1.id : fra2.id;

		return fra1.id < fra2.id ? fra1.id : fra2.id;
	}


	@XmlElement(name = "smirks")
	public String getSmirks() {
		long[] frags = new long[]{frag1Id, frag2Id};
		if (inversed)
			frags = new long[]{frag2Id, frag1Id};
		return MMPFragment.getFragment(frags[0]).smiles.replaceAll("\\[Al\\]", "*") + " -> " + MMPFragment.getFragment(frags[1]).smiles.replaceAll("\\[Al\\]", "*");
	}

	@XmlElement(name = "pair")
	public List<MMPair> getPairs() {
		if (!Globals.getMarshallingOption(MarshallingOption.PAIR_TRANSFORMATION))
		{
			if (filteredPairs != null)
				return filteredPairs.subList(0, Math.min(10, filteredPairs.size()));
			else
				return Globals.session().createCriteria(MMPair.class).add(Restrictions.eq("transformation", this)).setMaxResults(10).list();
		}
		else
			return null;
	}

	@XmlElement
	private List<MMPTransformationAnnotation> getAnnotations() {
		if (!Globals.getMarshallingOption(MarshallingOption.TRANSFOMATION_ANNOTATIONS))
			return null;

		MMPAnnotationSet def = new MMPAnnotationSet();
		def.name = "Default";
		def.id = 0L;
		Globals.session().evict(def);		

		List<MMPTransformationAnnotation> l =
				Globals.session().createCriteria(MMPTransformationAnnotation.class)
				.createAlias("annotationSet", "aset", JoinType.LEFT_OUTER_JOIN)
				.add(Restrictions.eq("transformation", this))
				.add(Restrictions.or(
						Restrictions.isNull("annotationSet"),
						Restrictions.eq("aset.user", Globals.userSession().user)
						))
						.addOrder(Order.asc("annotationSet"))
						.list();
		for (MMPTransformationAnnotation transformationAnnotation : l) 
			if (transformationAnnotation.annotationSet == null)
			{
				transformationAnnotation.annotationSet = def;
				Globals.session().evict(transformationAnnotation);
			}

		return l;
	}

	/**
	 * Used for XML only. Returns the number of relevant MMPs (total or filtered)
	 */
	@XmlElement
	private long getPairsCount()
	{
		if (filteredPairs != null)
			return filteredPairs.size();
		else
			return pairsCount;
	}

	public static void clearCache() {
		cache.clear();
	}

	public String toString(){
		return getSmirks();
	}
	
	/**
	 * An ugly ugly workaround around the Hibernate impotency to use several association multiple times.
	 * Can be refactored to HQL later.
	 * 
	 * Sergii - consider re-implementing this local piece using JPA criterias
	 */
	@OneToMany(mappedBy = "transformation", cascade = CascadeType.ALL, fetch = FetchType.LAZY)
	@XmlTransient
	private List<MMPTransformationAnnotation> annotations1 = new ArrayList<MMPTransformationAnnotation>();

	@OneToMany(mappedBy = "transformation", cascade = CascadeType.ALL, fetch = FetchType.LAZY)
	@XmlTransient
	private List<MMPTransformationAnnotation> annotations2 = new ArrayList<MMPTransformationAnnotation>();

	@OneToMany(mappedBy = "transformation", cascade = CascadeType.ALL, fetch = FetchType.LAZY)
	@XmlTransient
	private List<MMPTransformationAnnotation> annotations3 = new ArrayList<MMPTransformationAnnotation>();

	@OneToMany(mappedBy = "transformation", cascade = CascadeType.ALL, fetch = FetchType.LAZY)
	@XmlTransient
	private List<MMPTransformationAnnotation> annotations4 = new ArrayList<MMPTransformationAnnotation>();

	@OneToMany(mappedBy = "transformation", cascade = CascadeType.ALL, fetch = FetchType.LAZY)
	@XmlTransient
	private List<MMPTransformationAnnotation> annotations5 = new ArrayList<MMPTransformationAnnotation>();

	@OneToMany(mappedBy = "transformation", cascade = CascadeType.ALL, fetch = FetchType.LAZY)
	@XmlTransient
	private List<MMPTransformationAnnotation> annotations6 = new ArrayList<MMPTransformationAnnotation>();

}
