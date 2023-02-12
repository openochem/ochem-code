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
import java.io.Serializable;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.TreeMap;

import javax.persistence.CascadeType;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.FetchType;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.JoinTable;
import javax.persistence.ManyToMany;
import javax.persistence.ManyToOne;
import javax.persistence.OneToMany;
import javax.persistence.OrderBy;
import javax.persistence.Transient;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlElementWrapper;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.FlushMode;
import org.hibernate.Hibernate;
import org.hibernate.annotations.Formula;
import org.hibernate.criterion.Disjunction;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.ThreadScope;
import qspr.annotations.Loggable;
import qspr.dao.Repository;
import qspr.dao.Various;
import qspr.metaserver.util.MixtureAttachment;
import qspr.util.AccessChecker;
import qspr.util.CacheMap;
import qspr.util.HashedEntity;
import qspr.util.UserContributedEntity;
import qspr.util.unitconversion.UnitConversion;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.NumericalValueStandardizer;
import com.eadmet.utils.OCHEMUtils;


//@Filters({
//	@Filter(name="userFilter", condition="owner_id = :userId")
//})
@SuppressWarnings("unchecked")
@Entity
@Loggable
@XmlRootElement(name = "exp-property")
public class ExperimentalProperty implements Serializable, UserContributedEntity, HashedEntity<ExperimentalProperty>
{
	private static transient final Logger logger = LogManager.getLogger(ExperimentalProperty.class);
	private static final long serialVersionUID = 1L;

	public static final int NAME_NOT_CHECKED = 0;
	public static final int NAME_CHECKED_I1TRUE_I2TRUE = 1;
	public static final int NAME_CHECKED_I1TRUE_I2FALSE = 2;
	public static final int NAME_NOT_THERE = 3;
	public static final int NAME_CHECKED_I1FALSE_I2FALSE = 4;

	// Statistics
	@XmlTransient
	static public CacheMap<MixtureAttachment> mixtureCache = new CacheMap<MixtureAttachment>();

	@Id
	@GeneratedValue
	@Column(name = "exp_property_id")
	@XmlAttribute
	public Long id;

	@ManyToOne
	@JoinColumn(name = "molecule_id")
	@OrderBy
	public Molecule molecule;

	@ManyToOne
	@JoinColumn(name = "property_id")
	public Property property;

	@ManyToOne
	@JoinColumn(name = "unit_id")
	public Unit unit;

	@Formula("IF(first_entry = exp_property_id, 1, 0)")
	public boolean isPrimary;

	//##	@OneToMany(fetch = FetchType.LAZY, cascade = CascadeType.ALL, mappedBy = "ep")
	//##	@Cascade(org.hibernate.annotations.CascadeType.DELETE_ORPHAN)
	//##	@XmlElementWrapper(name="names")
	//##	public List<NameFact> name = new ArrayList<NameFact>();

	@ManyToMany
	(
			targetEntity = MoleculeName.class,
			cascade={CascadeType.PERSIST, CascadeType.MERGE},
			fetch = FetchType.LAZY
			)
	@JoinTable
	(
			name="ExperimentalPropertyName",
			joinColumns={@JoinColumn(name="exp_property_id")},
			inverseJoinColumns={@JoinColumn(name="molecule_name_id")}
			)
	@XmlTransient
	public List<MoleculeName> moleculenames = new ArrayList<MoleculeName>();

	@XmlTransient
	@OneToMany(mappedBy = "ep", cascade = CascadeType.ALL, fetch = FetchType.LAZY)
	public List<ColoredName> colorednames = new ArrayList<ColoredName>();

	@ManyToOne
	@JoinColumn(name = "article_id")
	@Loggable
	@XmlTransient
	public Article article;

	@ManyToOne
	@JoinColumn(name = "exp_property_id_conn")
	@XmlTransient
	public ExperimentalProperty connectedProperty;

	@ManyToOne
	@JoinColumn(name = "con_set_id")
	@XmlTransient
	public ConditionSet conditions;

	@Column(name = "art_page_num")
	@XmlElement(name = "art-page-num")
	@Loggable(name = "page number")
	public Integer artPageNum;

	@Column(name = "art_table_num")
	@XmlElement(name = "art-table-num")
	@Loggable(name = "table number")
	public String artTableNum;

	@Column(name = "art_line_num")
	@XmlElement(name = "art-line-num")
	@Loggable(name = "line number")
	public Integer artLineNum;

	@Column(name = "art_mol_id")
	@XmlElement(name = "art-mol-id")
	@Loggable(name = "molecule identifier")
	public String artMolId;

	@Column(name = "art_paragraph")
	@Loggable(name = "paragraph")
	public Integer artParagraph;

	@Column
	@XmlElement
	public String other; // is that maybe a comment? // why not, should be

	@Column(name = "error_comment")
	@XmlElement
	public String errorComment;

	public static final int STATUS_ERROR = 0;
	public static final int STATUS_TOVERIFY = 1;
	public static final int STATUS_TOREFERENCE = 2;
	public static final int STATUS_INVALID = 3; // Status for automatically generated errors (duplicates, no obligatory conditions etc..)

	@Column(name = "status")
	public Integer ep_status; // null, 0 or 3 = error in record, 0 is for Trash 1 = referenced by other record (to be verified)

	@Column
	@XmlAttribute
	public Integer rights = Globals.RIGHTS_FREELY_AVAILABLE; //

	@ManyToOne
	@JoinColumn(name = "modifier_id")
	@XmlTransient
	@Loggable(name = "modifier")
	public User owner; // This is actually _modifier_, not an owner. The old name stayed

	@ManyToOne
	@JoinColumn(name = "introducer_id")
	@XmlTransient
	public User introducer;


	@Column(name = "ep_md5", unique = true)
	@XmlAttribute
	@Loggable(exclude = true)
	public String md5;

	@Column
	@Loggable(exclude = true)
	public String hash;


	public double value;

	// A value, converted to the default unit and rounded to 3 significant digits
	@Column(name = "canonical_value")
	public String canonicalValueOld;

	@Column(name = "canonical_value_num")
	public Double canonicalValue;

	// Stands for either accuracy or right bound of interval, dependent on predicate / Midnighter
	@Column(name = "accuracy")
	public Double secondValue;

	@ManyToOne
	@JoinColumn(name = "predicate_id")
	public Predicate predicate;


	@OneToMany(fetch = FetchType.LAZY, mappedBy = "ep")
	@XmlTransient
	public List<BasketEntry> basketEntries;

	///
	@OneToMany(fetch = FetchType.LAZY, mappedBy = "ep")
	@XmlTransient
	public List<BasketEntry> selectedBasketEntries;
	////

	@XmlElement
	@Transient
	public String error;

	@Transient@XmlElement
	public Unit modelUnit;

	@Transient@XmlElement
	public String valueInModelUnits;

	@Transient@XmlElement
	public String predictedInOriginalUnits;

	@Transient
	@XmlElement
	public List<Long> knn; //Experimental k-nearest-neighbours in prediction space functionality for model/outlier interpretation

	@XmlElement
	@Transient
	public String info;

	@XmlAttribute
	@Transient
	public Integer count;

	@Column(name = "time_modified")
	@XmlTransient
	@Loggable(exclude = true)
	public Timestamp time;

	@Column(name = "time_created")
	@XmlTransient
	@Loggable(exclude = true)
	public Timestamp timeCreated;

	@Column
	@XmlTransient
	public Timestamp deleted;

	public Boolean approved = false;


	// "First entry" indicates the earliest published (in the sense of article's publication date) record,
	// that is a partial duplicate to the current one.
	// Calculated in a cron job
	@Column(name = "first_entry")
	@Loggable(exclude = true)
	public Long firstEntry;

	public boolean rejected = false;

	/**
	 * The date when the approval status (request, approval, rejection) has changed
	 */
	@Column(name = "moderation_date")
	@XmlTransient
	public Timestamp moderationDate;

	@Column(name = "external_id")
	@Loggable(exclude = true)
	@XmlElement
	public String externalId;

	// TODO: Sergey DOT PROFILE SNIPPET! TO BE REFGACTORED
	@XmlElement
	@Transient
	public String real;

	@XmlElement
	@Transient
	public String predicted;

	@XmlElement
	@Transient
	private DMValue dmValue;

	@XmlElement
	@Transient
	public DataTable descriptorList;

	@XmlElement
	@Transient
	public String correl;

	@XmlElement
	@Transient
	public List<Float> predictionVector;



	//
	//         DOT PROFILE SNIPPET END
	//
	// TODO: Sergey This piece in only needed in BatchUpload... can be moved to EPProxy but requires time
	//	@XmlElementWrapper(name="excelStrings")
	//	@XmlElement(name="excelString")
	//	@Transient
	//	public List<String> excelStrings = new ArrayList<String>();
	//
	//	@XmlAttribute
	//	@Transient
	//	public Integer index;
	//
	//	@XmlAttribute
	//	@Transient
	//	public Integer state;
	//
	//	@XmlAttribute
	//	@Transient
	//	public Integer status;
	//
	//	@XmlAttribute
	//	@Transient
	//	public Long record_id;
	//
	//	@XmlElement
	//	@Transient
	//	public ExperimentalProperty original_record;
	//
	//

	@Transient
	@XmlElement
	public ExperimentalProperty duplicate;

	@ManyToOne
	@JoinColumn(name = "poption_id")
	public PropertyOption option;

	@OneToMany (fetch = FetchType.LAZY, mappedBy = "ep")
	@XmlTransient
	public List<BatchUploadRow> batchUploadRows;

	@Transient@XmlElement
	public Boolean exclude;

	@XmlAttribute(name = "connected_id")
	public Long getConnectedId()
	{
		if (connectedProperty != null && !ThreadScope.get().controller.equals("molbrowser"))
			return connectedProperty.id;
		else
			return null;
	}

	@XmlElement(name="article")
	public Article getArticle()
	{
		if(!ThreadScope.get().controller.equals("molbrowser"))
			return this.article;
		return null;
	}

	@XmlElement(name="conditions")
	public ConditionSet getConditions()
	{
		if(!ThreadScope.get().controller.equals("molbrowser"))
			return this.conditions;
		return null;
	}

	public List<AbstractMoleculeName> getCheckedNamesIndividual()
	{
		List<AbstractMoleculeName> tempList = new ArrayList<AbstractMoleculeName>();
		determineColorForName();

		if (moleculenames != null)
			for (MoleculeName name : moleculenames)
				tempList.add(name);
		return tempList;
	}

	@XmlElementWrapper(name="moleculenames")
	@XmlElement(name = "moleculename")
	public List<AbstractMoleculeName> getCheckedNames()
	{
		List<AbstractMoleculeName> tempList = new ArrayList<AbstractMoleculeName>();
		if (colorednames != null)
			if (this.id != null && this.id > 0)
			{
				for (ColoredName name : colorednames)
					tempList.add(name);
			} else
				tempList = getCheckedNamesIndividual();

		return tempList;
	}

	public void setValueWithPredicate(String newValue) throws Exception
	{

		try
		{
			Object[] vals = PropertyValue.parseTextIntoPredicateAndValues(newValue);
			predicate = (Predicate)vals[0];
			value = (Double)vals[1];
			secondValue = (Double)vals[2];
		}
		catch (Exception e)
		{
			if (newValue.trim().equals(""))
				throw new UserFriendlyException("No value for property \""+property.getName()+"\" provided");
			else
				throw new UserFriendlyException("Invalid value for property \""+property.getName()+"\": "+newValue);
		}
	}

	public void setOption(PropertyOption option)
	{
		assert property.isQualitative();
		this.option = option;
	}


	@XmlAttribute
	public boolean isSelected()
	{
		//Set<Long> selectionSet = (Set<Long>)ThreadScope.get().localRequest.getSession().getAttribute("selection-list");
		//Set<Long> selectionList = Globals.userSession().selectionList;
		if (id != null)
			return Globals.selectionBasket("1".equals(ThreadScope.get().localRequest.getParameter("trash"))).containsEntry(this.id);
		else
			return false;
		//		if (basketEntries != null)
		//			return basketEntries.size() > 0;
		//			else
		//				return false;
	}

	//@XmlTransient
	public Double getConvertedValue(Unit targetUnit)
	{
		return getConvertedValue(value, targetUnit);
	}


	public Double getConvertedAverageValue(Unit targetUnit)
	{
		if (!predicate.isInterval())
			return getConvertedValue(value,  targetUnit);
		else
			return getConvertedValue((value+secondValue)/2D,  targetUnit);
	}

	public Double getConvertedValue(Double valueToConvert, Unit targetUnit)
	{
		return UnitConversion.convert(valueToConvert, unit, targetUnit, molecule.molWeight);
	}

	public void approve()
	{
		approved = true;
		rejected = false;
		if (!property.approved)
			throw new UserFriendlyException("You are trying to approve a record for a yet unapproved property " + property.getName() + ".\n Please, make the property approved first.");
		if (property.moderator != null && !property.moderator.equals(Globals.userSession().user))
			throw new UserFriendlyException("You are not allowed to approve a record for the property " + property.getName() + ", which is moderated by " + property.moderator.login);
		ThreadScope.get().recordApproved.fire(this);
	}

	public void updateHash()
	{
		// Unique index for avoiding duplicate entries
		// Article + Property + Value + Molecule (ID2) + Disc. conditions /  Midnighter
		if (property.isTextual())
			throw new RuntimeException("Arbitrary text properties are not supported");

		if (property.isDirectory)
			throw new UserFriendlyException("Can't assign a record to a group of properties " + property.getName());

		if (!property.isPublished() && this.isPublished())
		{
			//property.publish();
			this.rights = Globals.RIGHTS_NONE; // the record cannot be made publicly available unless property is publicly available
		}

		String stalled_units ="";

		if (property.isNumeric())
		{
			try
			{
				double val = UnitConversion.convert(value, unit, null, molecule.molWeight, NumericalValueStandardizer.SIGNIFICANT_DIGITS_CANONICAL);
				canonicalValueOld = NumericalValueStandardizer.getSignificantDigitsStr(val, NumericalValueStandardizer.SIGNIFICANT_DIGITS_CANONICAL);
			}
			catch (Exception e)
			{
				canonicalValueOld = NumericalValueStandardizer.getSignificantDigitsStr(value,  NumericalValueStandardizer.SIGNIFICANT_DIGITS_CANONICAL);
				stalled_units = unit + "_";
			}
			canonicalValue = Double.valueOf(canonicalValueOld);

			if (!unit.category.equals(property.unitCategory)) {
				ep_status = ExperimentalProperty.STATUS_ERROR;
				errorComment = " unit.category:" + unit.category + " != property.unitCategory:" + property.unitCategory;
			}

		} else
		{
			canonicalValue = 0D;
			canonicalValueOld = "0";
			if(option == null){
				ep_status = ExperimentalProperty.STATUS_ERROR;
				errorComment = "Condition option is null";
			}
		}

		String oldMD5 = md5;
		boolean markedAsError = ep_status != null && (ep_status == ExperimentalProperty.STATUS_ERROR || ep_status == ExperimentalProperty.STATUS_INVALID);

		if (markedAsError || (deleted != null))
			md5 = null;
		else
		{
			errorComment = null;

			StringBuffer hash = new StringBuffer();

			hash.append("" + article.id + "_" + property.id + "_");

			// Value
			hash.append(stalled_units+(property.isQualitative() ? option.id : canonicalValue) + "_");

			// Molecule
			if (molecule.isEmptyMolecule())
				if (moleculenames != null && moleculenames.size() > 0)
					hash.append(moleculenames.get(0).id);
				else
					if (externalId != null)
						hash.append(externalId);
					else
						hash.append("noname");
			else
				hash.append(molecule.mapping2.id);

			hash.append("_");

			String oblConditionsHash = getObligatoryConditionHash(property);

			if (oblConditionsHash.startsWith(QSPRConstants.ERROR)) {
				md5 = null;
				errorComment = "Some obligatory conditions has not been specified: " + oblConditionsHash.replace(QSPRConstants.ERROR, "");
				ep_status = ExperimentalProperty.STATUS_ERROR;
			}else {
				hash.append(oblConditionsHash);
				this.hash = hash.toString();
				md5 = OCHEMUtils.getMD5(hash.toString());
			}
		}

		if(timeCreated == null || time == null) {
			if(timeCreated == null)timeCreated = new Timestamp(Calendar.getInstance().getTimeInMillis());
			if(time == null)time = timeCreated;
		}

		if(isDummyOrEmpty())firstEntry = id; // dummy or empty record should always point to itself
		else
			if (!("" + md5).equals("" + oldMD5)) // for DUMMY firstEntry does not make sense ...but duplicates should be avoided
			{
				// Something is changed. Recalculate the primary record
				firstEntry = null;

				// Recalculate primary records for the referencing records
				Globals.session().setFlushMode(FlushMode.MANUAL); // Avoid flushing the session with semi-changed record
				if (id != null)
					Globals.session().createQuery("update ExperimentalProperty set firstEntry=null where firstEntry=:id")
					.setParameter("id", id)
					.executeUpdate();
				Globals.session().setFlushMode(FlushMode.AUTO); //

				molecule.mapping1.visible = true;
			}

	}

	/**
	 * Create mixture attachment using specific syntaxes for MIXTURE annotation SMILES1;0.2 SMILES2;0.8 
	 * @param mixture
	 * @return
	 * @throws Exception
	 */
	//TODO add check that molecule == MIXTURE
	public static MixtureAttachment createMixtureAttachment(String mixture) throws Exception {
		if(mixtureCache.containsKey(mixture)) return mixtureCache.get(mixture);

		MixtureAttachment ma = new MixtureAttachment();

		String s[] = mixture.split("\\s+");

		boolean valuesProvided = false;
		for(int i=0; i<s.length; i++) {
			String smiles = s[i], value = "";
			if(smiles.contains(";")) {
				String splits[] = smiles.split(";");
				smiles= splits[0];
				value = splits[1];
				if(Double.valueOf(value) == 0)
					throw new IOException("Molar fraction should not be zero at: " + smiles + " Remove respective component.");
			}
			String inchie = Various.molecule.getInChiKeyNoStereoDeSault(smiles);
			if(i == 0) valuesProvided = value.length() != 0;
			if(valuesProvided != (value.length() != 0))throw new IOException("Syntax error: some SMILES have and some do not have provided molar fractions.");

			if(!valuesProvided) continue;

			double val = Double.valueOf(value);
			if(ma.fractions.containsKey(inchie))
				ma.fractions.put(inchie, val + ma.fractions.get(inchie)); // accumulating for the same structure
			else {
				ma.fractions.put(inchie, val);
				ma.smiles.add(smiles);
			}
		}

		if(valuesProvided){  // check whether they are all zero
			valuesProvided = false;
			double val = -1.;
			for(Double v: ma.fractions.values()) {
				if(val == -1.) val = v;
				if(v != val) valuesProvided = true;
			}
		}

		if(!valuesProvided) ma.fractions = null; // all identical

		mixtureCache.put(mixture,ma);
		return ma;
	}

	/**
	 *  Creates solvent attachment using specific annotation for solvents NAME1; NAME2 (fraction1:fraction2)
	 * @param solvent
	 * @return
	 * @throws Exception
	 */

	public static MixtureAttachment createSolventAttachment(String solvent) throws Exception {
		if(mixtureCache.containsKey(solvent)) return mixtureCache.get(solvent);

		System.out.println(solvent);

		MixtureAttachment ma = new MixtureAttachment();
		ma.fractions = new TreeMap<String,Double>();

		String s[] = solvent.split(";"); // split by name
		double components[] = null;

		if(s.length > 1) {
			if(s[s.length-1].contains("(")){
				String fraction[] = s[s.length-1].split("\\(");
				s[s.length-1] = fraction[0];
				fraction[1] = fraction[1].replace(")", "");
				fraction = fraction[1].split(":");
				if(fraction.length>1) {
					if(fraction.length != s.length)
						throw new IOException("Some components have and some do not have molar fractions: " + solvent);
					components = new double[fraction.length];
					double sum =0;
					for(int i=0;i<fraction.length;i++) 
						sum += components[i]=Double.valueOf(fraction[i]);
					for(int i=0;i<fraction.length;i++) components[i] /= sum; // normalization
				}
			}
		}

		for(int i=0; i<s.length; i++) {
			String name = s[i].trim();

			Molecule mol = 	 Repository.molecule.getBySolventName(name);

			if(mol == null) throw new IOException("Could not find solvent with name: " + name);

			double val = components == null ? 1./s.length: components[i];
			if(ma.fractions.containsKey(mol.mapping1.inchi1))
				throw new IOException("Two identical solvents are not allowed : " + name);
			else {
				ma.fractions.put(mol.mapping1.inchi1, val);
				ma.smiles.add(Various.molecule.convertToFormatFixMetal(mol.getData(),QSPRConstants.SMILESH));
			}

		}

		mixtureCache.put(solvent,ma);
		return ma;
	}

	public String getObligatoryConditionHash(Property property) {

		if(property.obligatoryConditions.size() == 0)return "";

		StringBuffer hash = new StringBuffer();

		// Consider discriminate conditions
		HashSet<Property> obligatoryConditions = new HashSet<Property>();
		obligatoryConditions.addAll(property.obligatoryConditions);
		if (conditions != null)
		{
			// Should have never been an issue... But for two duplicate records with same condition set the hashes where different due to different order of condition values
			Collections.sort(conditions.values, PropertyValue.condNameComp);
			for (PropertyValue conditionEntry : conditions.values)
			{
				if (obligatoryConditions.contains(conditionEntry.property))
				{
					if (conditionEntry.property.isQualitative())try{

						if(conditionEntry.property.isMapable() != null)
							hash.append(createSolventAttachment(conditionEntry.option.name).toString());
						else
							if(conditionEntry.property.getName().equals(QSPRConstants.MIXTURE_CONDITION))
								hash.append(createMixtureAttachment(conditionEntry.option.name).toString());
							else
								hash.append(conditionEntry.option.id);

					}catch(Exception e) {
						return QSPRConstants.ERROR + " failed to resolve: "+ conditionEntry.option.name;
					}

					else if (conditionEntry.property.isNumeric())
						try
					{
							double val = UnitConversion.convert(conditionEntry.value, conditionEntry.unit, null, molecule.molWeight, NumericalValueStandardizer.SIGNIFICANT_DIGITS_CANONICAL);
							hash.append(NumericalValueStandardizer.getSignificantDigitsStr(val, NumericalValueStandardizer.SIGNIFICANT_DIGITS_CANONICAL));
					}
					catch (Exception e)
					{
						hash.append(conditionEntry.unit +"_" + NumericalValueStandardizer.getSignificantDigitsStr(conditionEntry.value, NumericalValueStandardizer.SIGNIFICANT_DIGITS_CANONICAL));
					}
					else
						hash.append(conditionEntry.textualValue.toLowerCase());
					hash.append("_");
					obligatoryConditions.remove(conditionEntry.property);
				}
			}
		}

		if (obligatoryConditions.size() != 0)
			return QSPRConstants.ERROR + obligatoryConditions;

		return hash.toString();
	}

	@XmlAttribute(name="autoEvidence")
	public Long getAutoCheckEvidence()
	{
		return firstEntry == null || (long)firstEntry == (long)id? null:firstEntry;
	}

	public boolean hasConflicts()
	{
		duplicate = null;
		List<ExperimentalProperty> eps = getAllDublicates();
		if (eps != null)
			duplicate = eps.get(0);
		return (duplicate != null);
	}

	public List<ExperimentalProperty> getAllDublicates()
	{
		// Well, the bad thing is
		// this table is expected to be big enough
		// and searches by MD5 before inserting each item could
		// possibly be very slow
		// Well, at the moment let it be, yeah?
		// Midnighter

		if (md5 == null || deleted != null) return null;

		Globals.session().setFlushMode(FlushMode.MANUAL); //Added to avoid session flush before query to database. Session flush may result in dublicate constraint violation.
		Criteria criteria = Globals.session().createCriteria(ExperimentalProperty.class).add(Restrictions.eq("md5", md5));
		if (id != null)
			criteria.add(Restrictions.ne("id", this.id));

		List<ExperimentalProperty> eps = criteria.list(); // in principle maximum one is expected ... but
		Globals.session().setFlushMode(FlushMode.AUTO);

		return eps.size()>0?eps:null;
	}

	public void resolveConflictsAndSave()
	{

		List<ExperimentalProperty> duplicates = getAllDublicates();

		if (duplicates == null)
		{	
			Globals.session().saveOrUpdate(this);
			return;
		}

		duplicates.add(0, this);

		ExperimentalProperty first = this;

		for(ExperimentalProperty dup:duplicates) {
			if((long)dup.id == (long)first.id)continue; // same
			if(first.newIsBetter(dup)) first = dup;
		}

		for(ExperimentalProperty dup:duplicates) {
			if((long)dup.id == (long)first.id) { // record is OK and first; will be updated last
				dup.ep_status = null;
				dup.errorComment = null;
				continue;
			}

			dup.ep_status = ExperimentalProperty.STATUS_INVALID;
			dup.errorComment = "Duplicate of record (RecordID: R"+first.id+", User: "+ first.introducer.login+")";
			dup.md5 = null;
			Globals.session().saveOrUpdate(dup);
		}
		first.updateHash();
		Globals.session().saveOrUpdate(first); // only now, first we need to invalidate the previous md5

	}

	public boolean newIsBetter(ExperimentalProperty newone) {
		if(newone == null || newone.rights == null || newone.article.publicationDate == null) return false;
		if((int)rights != (int)newone.rights)return  (int)newone.rights > (int)rights;
		if(!article.publicationDate.equals(newone.article.publicationDate))
			return article.publicationDate.after(newone.article.publicationDate);
		if(newone.time != null && // TODO to remove once all time are not null
				rights == Globals.RIGHTS_FREELY_AVAILABLE && !time.equals(newone.time))
			return time.after(newone.time);
		return (long)newone.id < (long)id; //less id
	}

	@XmlAttribute(name="visible")
	public String getVisible()
	{
		if (rights == Globals.RIGHTS_FREELY_AVAILABLE)
			return "all";

		if (owner == null)
			return "none";

		if (owner.group == null)
			return (owner.login);

		return (owner.group.name);

	}

	public boolean isPublic()
	{
		return rights == Globals.RIGHTS_FREELY_AVAILABLE;
	}

	public boolean isDummyOrEmpty()
	{
		if((long) property.id == Property.getDummyId())return true;
		if(molecule.isEmptyMolecule())return true;
		return false;
	}

	@XmlElement(name = "owner")
	public String getOwnerStr()
	{
		if(!ThreadScope.get().controller.equals("molbrowser"))
			if (owner != null)
				return owner.login;
			else
				return null;
		return null;
	}

	@XmlElement(name = "introducer")
	public String getIntroducerStr()
	{
		if(!ThreadScope.get().controller.equals("molbrowser"))
			if (introducer != null)
				return introducer.login;
			else
				return null;
		return null;
	}

	@XmlElement(name = "readonly")
	public String getReadonly()
	{
		User curUser = Globals.userSession().user;
		if (owner != null && rights < Globals.RIGHTS_FREELY_AVAILABLE && owner != curUser)
			return (curUser == null) ? "true" :
				(owner.rank >= curUser.rank) ?
						"true" : null;
		return null;
	}

	//NoS 18.10.2011 - Naive method to "merge" duplicates in Batch Upload - all empty non-critical fields in original will be filled with non-empty fields of the argument. Hash should not change in this operation.
	public ExperimentalProperty merge(ExperimentalProperty ep)
	{
		artLineNum = (artLineNum == null) ? ep.artLineNum : artLineNum;
		artMolId = (artMolId == null || artMolId.equals("") || artMolId.startsWith("AUTO_")) ? ep.artMolId : artMolId;
		artPageNum = (artPageNum == null) ? ep.artPageNum : artPageNum;
		artParagraph = (artParagraph == null) ? ep.artParagraph : artParagraph;
		artTableNum = (artTableNum == null || artTableNum.equals("")) ? ep.artTableNum : artTableNum;

		conditions = (conditions == null) ? conditions : conditions.mergeAddOrUpdateWith(ep.conditions); // Non-trivial condition merge procedure...

		connectedProperty = (connectedProperty == null) ? ep.connectedProperty : connectedProperty;
		info = (info == null || info.equals("")) ? ep.info : info;

		for (MoleculeName mn : ep.moleculenames)
			if (!moleculenames.contains(mn))
				moleculenames.add(mn);

		externalId = (externalId == null || externalId.equals("")) ? ep.externalId : externalId;
		other =  (other == null || other.equals("")) ? ep.other : other;
		//		ep_status = (ep_status == null) ? ep.ep_status : ep_status;

		//		//Merge into a basket
		//		if (ep.basketEntries != null)
		//		{
		//			for (BasketEntry be : ep.basketEntries)
		//			{
		//				BasketEntry nbe = new BasketEntry();
		//				nbe.basket = be.basket;
		//				nbe.ep = this;
		//	
		//				if (basketEntries == null)
		//					basketEntries = new ArrayList<BasketEntry>();
		//	
		//				if (!basketEntries.contains(nbe))
		//					basketEntries.add(nbe);
		//			}
		//		}
		//
		updateHash();
		return this;
	}

	public void addName(String name)
	{
		MoleculeName molName = MoleculeName.get(name);
		if (molName != null && !moleculenames.contains(molName))
			moleculenames.add(molName);
	}


	public void determineColorForName()
	{
		if (molecule == null)
			return;

		if (moleculenames != null)
			for (MoleculeName molName : moleculenames)
			{
				ValidatedFact vf = null;
				// default proof
				int foundProof = NAME_NOT_CHECKED;
				String userName =  null;

				Criteria c = Globals.session().createCriteria(ValidatedFact.class)
						.add(Restrictions.eq("validated", ValidatedFact.VALIDATED))
						.add(Restrictions.eq("moleculename", molName));
				// list of validatedfacts with same name
				@SuppressWarnings("rawtypes")
				List results = c.list();

				if (results.size() > 0)
				{
					vf = (ValidatedFact) results.get(0);
					if (vf.mapping == null)
					{
						foundProof = NAME_NOT_THERE;
					} else
					{
						// same structure then it's correctly validated
						if (molecule.mapping2.equals(vf.mapping))
						{
							foundProof = NAME_CHECKED_I1TRUE_I2TRUE;
							if(vf.source == 0)
							{

								User user = (User) Globals.session().get(User.getCurrentClass(), Long.valueOf(vf.sourceid));
								userName = user.login;
							}// same structure then it's correctly validated
						}
						else if (molecule.mapping1.equals(vf.mapping.mapping1)) {
							foundProof = NAME_CHECKED_I1TRUE_I2FALSE;
						}
						// otherwise it's validated as wrong structure-name pair
						else
						{
							foundProof = NAME_CHECKED_I1FALSE_I2FALSE;
						}
					}
				}


				molName.validation = foundProof;
				if(userName != null)
					molName.user = userName;
			}
	}

	public MoleculeName containsName(String name)
	{
		for (MoleculeName molname : moleculenames)
		{
			if (molname.name.equals(name))
				return molname;
		}
		return null;
	}

	public ExperimentalProperty cloneForReference() {
		ExperimentalProperty clone = new ExperimentalProperty();
		clone.property = property;
		clone.value = value;
		clone.secondValue = secondValue;
		clone.predicate = predicate;
		clone.unit = unit;
		clone.molecule = molecule;
		clone.error = error;
		clone.conditions = conditions;
		clone.option = option; //12.01.09 novserj
		//		clone.artLineNum = artLineNum;
		//		clone.artMolId = artMolId;
		//		clone.artPageNum = artPageNum;
		//		clone.artParagraph = artParagraph;
		//		clone.artTableNum = artTableNum;
		clone.moleculenames = moleculenames;
		clone.rights = rights;
		return clone;
	}

	public ExperimentalProperty fullClone() {
		ExperimentalProperty clone = new ExperimentalProperty();
		clone.property = property;
		clone.value = value;
		clone.secondValue = secondValue;
		clone.predicate = predicate;
		clone.unit = unit;
		clone.molecule = molecule;
		clone.error = error;
		clone.conditions = conditions;
		clone.option = option;
		clone.artLineNum = artLineNum;
		clone.artMolId = artMolId;
		clone.artPageNum = artPageNum;
		clone.artParagraph = artParagraph;
		clone.artTableNum = artTableNum;
		clone.moleculenames = moleculenames;
		clone.rights = rights;
		clone.article = article;
		return clone;
	}

	public ExperimentalProperty getReferencedProperty() {
		Criteria criteria = Globals.session().createCriteria(ExperimentalProperty.class)
				.add(Restrictions.isNull("deleted"))
				.add(Restrictions.eq("property", this.property))
				.add(Restrictions.eq("value", this.value))
				.add(Restrictions.eq("article", this.article))
				//			.add(Restrictions.isNotNull("connectedProperty"))
				.createCriteria("molecule")
				.add(Restrictions.eq("mapping2", this.molecule.mapping2));

		@SuppressWarnings("rawtypes")
		List properties = criteria.list();
		if (properties.size() > 0){
			if (properties.size() > 1)
				logger.info("*****************************Warning in getReferencedProperty");
			return (ExperimentalProperty) properties.get(0);
		} else {
			return this;
		}
	}

	@XmlElement
	public String getTime()
	{
		if(!ThreadScope.get().controller.equals("molbrowser"))
			return (time != null)?ThreadScope.get().fullDateFormat.format(time):"";
		return null;
	}

	@XmlElement(name = "time-created")
	public String getTimeCreated()
	{
		if(!ThreadScope.get().controller.equals("molbrowser"))
			return (timeCreated != null)?ThreadScope.get().fullDateFormat.format(timeCreated):"";
		return null;
	}

	public String toString()
	{
		return property.getName()+" "+predicate+" "+value+" ("+id+")";
	}

	@XmlTransient
	public String getStringValue()
	{
		if (property.isQualitative())
			return option.name;
		if (property.isNumeric())
			return "" + value;
		return "";
	}

	@XmlTransient
	public String getPredicatedValue()
	{
		if (predicate == null || predicate.name.equals("="))
			return "" + value;
		else if (predicate.isAccuracy() || predicate.isInterval())
			return value + " " + predicate.name + " " + secondValue;
		else
			return predicate.name + value;
	}

	@XmlElement(name = "molecule-tags")
	public List<Tag> getMoleculeTags()
	{
		if ("epbrowser".equals(ThreadScope.get().controller))
			if (molecule != null && molecule.mapping1 != null)
			{
				Disjunction rights = Restrictions.disjunction();
				rights.add(Restrictions.eq("isPublic", Boolean.TRUE));
				if (Globals.userSession().user != null)
					rights.add(Restrictions.eq("introducer", Globals.userSession().user));

				return Globals.session().createCriteria(Tag.class)
						.createAlias("mapping", "mp")
						.add(Restrictions.eq("mp.id", molecule.mapping1.id))
						.add(Restrictions.eq("showInBrowser", Boolean.TRUE))
						.add(rights).list();
				//return molecule.mapping1.tags;
			}
		return null;
	}

	public static void addAccessRestrictions(Criteria criteria, int requiredPrivileges, User accessingUser, boolean basketJoined)
	{
		boolean showUnapproved = !Globals.isGuestUser() && Globals.userSession().user.getUserSettings().showUnapprovedData();
		addAccessRestrictions(criteria, requiredPrivileges, accessingUser, basketJoined, showUnapproved);
	}


	public static void addAccessRestrictions(Criteria criteria, int requiredPrivileges, User accessingUser, boolean basketJoined, boolean showUnapproved)
	{
		Disjunction rightsRestriction = Restrictions.disjunction();

		// General restriction rules
		AccessChecker.addAccessRestrictions(criteria, requiredPrivileges, accessingUser, rightsRestriction, showUnapproved);

		// Specific restriction rules for ExperimentalProperty
		if (basketJoined)
		{
			if (Globals.userSession().visitedModelsIds != null && !Globals.userSession().visitedModelsIds.isEmpty())
			{
				// Also allow to see records, connected with models, visited by the user
				List<Long>  basketIds = Globals.session().createCriteria(Model.class).createAlias("trainingSet", "ts").add(Restrictions.in("publicId", Globals.userSession().visitedModelsIds)).setProjection(Projections.groupProperty("ts.id")).list();
				List<Long>  basketIds2 = Globals.session().createCriteria(Model.class).createAlias("validationSet", "ts").add(Restrictions.in("publicId", Globals.userSession().visitedModelsIds)).setProjection(Projections.groupProperty("ts.id")).list();
				if (!basketIds2.isEmpty())
					basketIds.addAll(basketIds2);

				if (!basketIds.isEmpty())
					rightsRestriction.add(Restrictions.in("b.id", basketIds));
			}

			// If I have private records in MY basket, show them to me
			if (!Globals.isGuestUser())
				rightsRestriction.add(Restrictions.ge("b.user", Globals.myself()));
		}
	}

	private transient String printableValue;
	@XmlElement(name = "printableValue")
	private String getPrintableValue()
	{
		if (printableValue != null)
			return printableValue;
		// Value, printable on the client

		//Bugfix/patch from 25.01.09 by NoS. This method was throwing null pointer exceptions when marshalling empty records (i.e. in EPBrowserController Action methods)
		if (property == null)
			return null;

		if (molecule == null)
			return null;

		if (property.isQualitative())
			return printableValue = " = " + (option == null?"null":option.name);

		String res;

		if (predicate == null)
			res = "= "+value;
		else
			if (predicate.shortName.equals("+-") || predicate.shortName.equals("-"))
				res = "= "+value+" "+predicate.name+" "+secondValue;
			else
				res = predicate.name+" "+value;

		res += " (in " + unit.getName() + ")";

		if (!property.defaultUnit.equals(unit))
		{
			try
			{
				double val = UnitConversion.convert(value, unit, property.defaultUnit, molecule.molWeight);
				res = res + " = " + formatNum(val) + " (in " + property.defaultUnit.getName() + ")";
			}
			catch (Exception e)
			{
				e.printStackTrace();
				res += " (conversion to " + property.defaultUnit.getName() +" is not possible)";
			}
		}

		return printableValue = res;
	}

	@XmlTransient
	public boolean isDeleted()
	{
		return deleted != null;
	}

	private String formatNum(double x)
	{
		//NoS 11.04.12 - Experimental; change displayed value for consistency of formats
		return NumericalValueStandardizer.getSignificantDigits(x);
		// To treat the weird behaivor of the double type, 204*0.01 = 0.20400000000000001
		//		DecimalFormat df;
		//		if (x >= 0.001)
		//			df = new DecimalFormat("#.####");
		//		else
		//			df = new DecimalFormat("#.####E0");
		//		return df.format(x);
	}

	public void initializeBasicCollections()
	{
		Article a = article;
		while (a != null)
		{
			Hibernate.initialize(a.pdfs);
			Hibernate.initialize(a.authors);
			a = a.parent;
		}
		Property p = property;
		while (p != null)
		{
			Hibernate.initialize(p.options);
			Hibernate.initialize(p.obligatoryConditions);
			p = p.parent;
		}
		if (conditions != null && conditions.values != null)
		{
			for (PropertyValue pv : conditions.values)
			{
				p = pv.property;
				while (p != null)
				{
					Hibernate.initialize(p.options);
					p = p.parent;
				}
			}
		}
		Hibernate.initialize(colorednames);
		Hibernate.initialize(moleculenames);
		Hibernate.initialize(basketEntries);
	}

	public void setDMValue(String dm, String value)
	{
		dmValue = new DMValue();
		dmValue.dm = dm;
		dmValue.value = value;
	}

	static class DMValue
	{
		@XmlAttribute
		public String dm;

		@XmlAttribute
		public String value;
	}

	public ExperimentalProperty()
	{

	}

	public ExperimentalProperty(Molecule m)
	{
		this.molecule = m;
	}

	@Override
	@XmlTransient
	public User getIntroducer() {
		return introducer;
	}

	@Override
	@XmlTransient
	public User getOwner() {
		return owner;
	}

	@Override
	@XmlTransient
	public Integer getRights()
	{
		return rights;
	}

	public ExperimentalProperty getDuplicate()
	{
		return duplicate;
	}

	@XmlTransient
	public boolean isPublished()
	{
		return rights == Globals.RIGHTS_FREELY_AVAILABLE;
	}

	public static void checkIfSetIsPublishable(List<Long> ids) {

		List<Property> props = Globals.session().createQuery("select pr.property from ExperimentalProperty pr where pr.id in (:ids) group by pr.property")
				.setParameterList("ids", ids)
				.list();

		for(Property p:props) {
			if(p.isPublished())continue;
			if(p.rights == Globals.RIGHTS_NONE)
				throw new UserFriendlyException("Property \"" + p.getName() + "\" is not yet public. Publish it before publishing the individual records.");
			if(!p.approved)
				throw new UserFriendlyException("Property \"" + p.getName() + "\" is not yet approved. Contact " + QSPRConstants.INFOEMAIL + " to speed its approval. You can publish records only once the property itself is approved.");
		}

		List<ExperimentalProperty> recs = Globals.session().createQuery("select pr from ExperimentalProperty pr where pr.id in (:ids) and pr.rejected=1 group by pr.rejected")
				.setParameterList("ids", ids)
				.list();

		if(recs.size()>0)
			throw new UserFriendlyException("Record R" +recs.get(0).id+ " was rejected by moderator. It cannot be published until comments of moderator are resolved."); 
	}


	/**
	 *  Publishes records which belong to public models
	 * @param ids
	 */

	public static void publishRecords(List<Long> ids, Long articleId) {
		if(ids ==null || ids.size() == 0)return;

		Globals.session().createSQLQuery("update ExperimentalProperty set rights = 3, modifier_id = " + QSPRConstants.PUBLISHER_ID +", time_modified =NOW() " + 
				" where exp_property_id in (:ids) and rights != 3 and modifier_id != " + QSPRConstants.PUBLISHER_ID).setParameterList("ids", ids)
		.executeUpdate();

		Globals.restartAllTransactions(true);

		Globals.session().createSQLQuery("update ExperimentalProperty join Article using (article_id) set ExperimentalProperty.article_id = " + articleId +
				" where exp_property_id in (:ids) and journal_id = " + QSPRConstants.UNPUBLISHED_JOURNAL).setParameterList("ids", ids)
		.executeUpdate();
	}

}
