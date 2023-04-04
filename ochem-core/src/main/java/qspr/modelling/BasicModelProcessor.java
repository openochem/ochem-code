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

package qspr.modelling;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.dao.Repository;
import qspr.dao.Various;
import qspr.entities.Basket;
import qspr.entities.BasketEntry;
import qspr.entities.ConditionSet;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Property;
import qspr.entities.PropertyOption;
import qspr.entities.PropertyValue;
import qspr.entities.Session;
import qspr.interfaces.ProvidedConditions;
import qspr.metaserver.configurations.MaximalSizeRestriction;
import qspr.metaserver.configurations.MultiLearningAbstractConfiguration;
import qspr.metaserver.util.MixtureAttachment;
import qspr.metaserver.util.ShortCondition;
import qspr.modelling.configurations.CDSConfiguration;
import qspr.modelling.configurations.ExternalCondition;
import qspr.util.ClassCompressor;
import qspr.util.unitconversion.UnitConversion;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.exceptions.UserFriendlyException;

// The abstract implementation of basic functions required for Consensus and CDS models
abstract public class BasicModelProcessor extends ModelProcessor
{
	protected DataTable dtTrainingSetDescriptors = null;
	//protected DataTable dtValidationSetDescriptors = null;
	List<Integer> numOfOptions;

	final static Integer MAX_EXCLUDED = 100000; // maximum number of excluded entries, large number crashes parsing and can be just a mistake
	final static private int TRANSACTION_RESTART = 1000;

	ProvidedConditions getModelConfiguration() {
		return (ProvidedConditions) model.attachment.getObject().configuration;
	}

	@Override
	protected void prepare() throws Exception
	{
		Repository.user.checkEligibility(Repository.basket.countEntries(model.trainingSet.id), QSPRConstants.MODEL_BONUS);
		super.prepare();
		enumeratePropertyOptions();

		ProvidedConditions consConf = getModelConfiguration();
		if(consConf.getConditions() != null)
			enumerateConditions(consConf.getConditions());
	}

	private void enumeratePropertyOptions()
	{
		if (numOfOptions == null || numOfOptions.isEmpty())
		{
			String hash = model.trainingSet.lastModified.toString()+model.trainingSet.id;

			OptionsEnumeration en = null;
			byte[] entry = Session.getHash(hash);
			if(entry != null) en = (OptionsEnumeration)ClassCompressor.byteToObject(entry);

			if (en == null) {
				en = OptionsEnumeration.enumeratePropertyOptions(model.getFilteredSet(QSPRConstants.TRAINING), mapper);
				Session.hashSet(hash,ClassCompressor.objectToByte(en));
			}
			else
				System.out.println("Reusing hashed value");

			if (!model.hasMappedOptions())
				model.attachment.getObject().optionsMapping = en.optionsMapping;
			numOfOptions = en.numOfOptions;

		}
	}

	@Override
	public Serializable getDataForApplier(Basket basket, ConditionSet defaultConditions) throws Exception
	{
		if (basket == null)
			throw new UserFriendlyException("Empty applying set provided");

		DataTable dt = ModelProcessor.basketToSDFTable(basket);

		for(int i=0; i<dt.getRowsSize();i++) // all data for applier are VALIDATION data
			dt.getRow(i).addAttachment(QSPRConstants.VALIDATION, Boolean.TRUE);

		WorkflowNodeData wnd = new WorkflowNodeData(dt);

		ProvidedConditions cond = getModelConfiguration();

		if (cond.hasConditions())
		{
			if(defaultConditions == null)defaultConditions = Repository.property.createDefaultConditions(cond, false); // no default descriptors were provided - restoring them

			// There are some external descriptors, for example conditions. Add them to the input for the teacher
			DataTable dtConditions = new DataTable();
			addConditions(basket, dt, dtConditions, defaultConditions);
			dtConditions.setId(QSPRConstants.CONDITIONS);
			wnd.addPort(dtConditions);
		}

		return wnd;
	}

	@Override
	public Serializable getDataForTeacher() throws Exception
	{
		// Teacher accepts 3 ports (two firsts are obligatory):
		// 1. Training set + Validation set
		// 2. Real (experimental) values for training set
		// 3. Obligatory conditions

		// 1a. Training set
		Basket trainingSet = model.getFilteredSet(QSPRConstants.TRAINING);
		if (trainingSet == null)
			throw new UserFriendlyException("Empty training set provided");

		ProvidedConditions cond = getModelConfiguration();

		String hash =  model.hashOfSets(cond);

		System.out.println("Baskets hash: " + hash);

		DataTable dtMolecules = null;
		byte[] entry = Session.getHash(hash);
		if(entry != null ) {
			dtMolecules = (DataTable) ClassCompressor.byteToObject(entry);			
			System.out.println("Found for Baskets hash: " + hash);
		}else
			System.out.println("Not Found for Baskets hash: " + hash);

		enumeratePropertyOptions(); // usually doing nothing, enumeration is done in prepare

		if(dtMolecules == null) { // not found, clearing previous data and creating new hash set

			dtMolecules = ModelProcessor.basketToSDFTable(trainingSet);

			setStatus("Preparing the training set molecules");

			dtMolecules.id = "mols";

			// 1b. Validation set(s)
			for (int i = 0; i < model.getValidationSets().size(); i++)
			{
				Basket validationSet = model.getFilteredSet(QSPRConstants.VALIDATION + i);
				if (validationSet != null)
					dtMolecules = ModelProcessor.basketToSDFTable(validationSet, dtMolecules);
			}

			// 1c. Excluded records
			Basket excludedSet = model.getFilteredSet(QSPRConstants.EXCLUDED);
			if (excludedSet != null){
				if(excludedSet.getRowsSize() > MAX_EXCLUDED)
					throw new UserFriendlyException("The number of excluded records " + 
							excludedSet.getRowsSize() + " above the allowed limit " + MAX_EXCLUDED +
							". Check prep-processing options for records with ranges."+
							" For Comprehensive Modeling these options can be selected in \"Show advanced options\"."
							);

				int record = dtMolecules.getRowsSize();
				dtMolecules = ModelProcessor.basketToSDFTable(excludedSet, dtMolecules);
				while(record<dtMolecules.getRowsSize())
					dtMolecules.getRow(record++).addAttachment(QSPRConstants.EXCLUDED, "1");
			}

			// 1d. Rearrange records to have additional mixture records at the end

			int training = trainingSet.entries.size() - 1;

			DataTable dtMoleculesNew = dtMolecules.getSlice(0, training);
			// real records
			HashSet<Integer> mols = new HashSet<Integer>(); // molecules for which descriptors will be calculated
			for(int i=training;i<dtMolecules.getRowsSize();i++)
				if(dtMolecules.getRow(i).getAttachment(QSPRConstants.ADDED_MOLECULES) == null) {
					dtMoleculesNew.addRow(dtMolecules.getRow(i));
					mols.add((Integer)dtMolecules.getRow(i).getAttachment(QSPRConstants.MOLECULE_ID_STEREOCHEM));
				}

			for(int i=training;i<dtMolecules.getRowsSize();i++) // adding those for which descriptors are otherwise absent
				if(dtMolecules.getRow(i).getAttachment(QSPRConstants.ADDED_MOLECULES) != null && !mols.contains(dtMolecules.getRow(i).getAttachment(QSPRConstants.MOLECULE_ID_STEREOCHEM))){
					dtMoleculesNew.addRow(dtMolecules.getRow(i));
					mols.add((Integer)dtMolecules.getRow(i).getAttachment(QSPRConstants.MOLECULE_ID_STEREOCHEM));
				}

			dtMolecules = dtMoleculesNew;

			// 1f. Mark non-training set entries as validation data points
			dtMolecules.currentRow = training;
			while (dtMolecules.nextRow())
				dtMolecules.getCurrentRow().addAttachment(QSPRConstants.VALIDATION, Boolean.TRUE);

			dtMolecules.attachment = Integer.valueOf(trainingSet.entries.size());

			// storing for caching
			Session.hashSet(hash,ClassCompressor.objectToByte(dtMolecules));
		}

		setStatus("Enumerating experimental values");

		Object configuration = model.attachment.getObject().configuration;

		if(configuration instanceof CDSConfiguration)
			configuration = ((CDSConfiguration)configuration).modelConfiguration;

		// 1d mark very large molecules as errors
		if(configuration instanceof MaximalSizeRestriction) {
			int maxMoleculeSize = ((MaximalSizeRestriction)configuration).getMaxSize();
			for(int i=0;i<dtMolecules.getRowsSize();i++) { // adding those for which descriptors are otherwise absent
				int size = Various.molecule.getAtomCount((String) dtMolecules.getValue(i, QSPRConstants.SDF_COLUMN));
				if(size <= maxMoleculeSize)continue;
				dtMolecules.getRow(i).setError("Number of atoms " + size + " > " + maxMoleculeSize + " (max atoms)");
			}
		}

		// 2. Real values
		DataTable dtProperties = ModelProcessor.basketToPropertiesTable(model.getFilteredSet(QSPRConstants.TRAINING), null, model, mapper);


		// Missed values handling
		Map<Long,Set<Integer>>  excludedImplicitMoleculesEntries = model.trainingSet.gexExcludedImplicitMolecules();
		if(configuration instanceof MultiLearningAbstractConfiguration 
				&& ((MultiLearningAbstractConfiguration)configuration).getImplicit() != null && excludedImplicitMoleculesEntries != null)

			for(int i = 0; i < dtProperties.getRowsSize() ; i++) {
				Double val[] = Arrays.copyOf(((MultiLearningAbstractConfiguration)configuration).getImplicit(),
						((MultiLearningAbstractConfiguration)configuration).getImplicit().length);

				for(int p = 0;  p<model.modelMappings.size(); p++){
					if(excludedImplicitMoleculesEntries.containsKey(model.modelMappings.get(p).property.id) &&
							excludedImplicitMoleculesEntries.get(model.modelMappings.get(p).property.id).
							contains(trainingSet.entries.get(i).ep.molecule.mapping2.id))
						val[p] = null; 
				}
				dtMolecules.getRow(i).addAttachment(QSPRConstants.IMPLICIT_ATTACHMENT, val);  // for both values and molecules
				dtProperties.getRow(i).addAttachment(QSPRConstants.IMPLICIT_ATTACHMENT, val); 
			}

		if( dtProperties.getRowsNoErrorsSize() != dtProperties.getRowsSize())
			for(int i = 0; i < dtProperties.getRowsSize() ; i++) 
				if(dtProperties.getRow(i).isError())
					dtMolecules.getRow(i).setError(dtProperties.getRow(i).detailedStatus);

		for (int i = 0; i < dtMolecules.getRowsSize() - trainingSet.entries.size(); i++)
			dtProperties.addStubRow();

		dtProperties.id = "exp-values";

		WorkflowNodeData wnd = new WorkflowNodeData(dtMolecules).addPort(dtProperties);

		// 3. Conditions


		if (cond.hasConditions())
		{
			setStatus("Enumerating conditions");

			trainingSet = model.getFilteredSet(QSPRConstants.TRAINING);
			DataTable dtConditions = new DataTable();
			ConditionSet defaultConditions = Repository.property.createDefaultConditions(cond, true);
			addConditions(trainingSet, dtMolecules,  dtConditions, defaultConditions);
			for (int i = 0; i < model.getValidationSets().size(); i++)
			{
				Basket validationSet = model.getFilteredSet(QSPRConstants.VALIDATION + i); 
				if (validationSet != null)
					addConditions(validationSet, dtMolecules, dtConditions, defaultConditions);
			}
			if (model.getFilteredSet(QSPRConstants.EXCLUDED) != null)
				addConditions(model.getFilteredSet(QSPRConstants.EXCLUDED), dtMolecules, dtConditions, defaultConditions);

			dtConditions.setId("conditions");

			while(dtConditions.getRowsSize() < dtMolecules.getRowsSize())
				dtConditions.addRow(); // adding to have the same number of conditions

			wnd.addPort(dtConditions);
		}

		return wnd;
	}

	private void addConditions(Basket basket, DataTable data, DataTable dtConditions, ConditionSet defaultConditions)
	{
		int i = 0;

		if(defaultConditions == null) throw new UserFriendlyException("Default conditions have not been provided!");

		ProvidedConditions conf = getModelConfiguration(); 

		List<ExternalCondition> eDescs = new ArrayList<ExternalCondition>();

		for (ShortCondition eDesc : conf.getConditions())
			eDescs.add(eDesc instanceof ExternalCondition ?  (ExternalCondition) eDesc : new ExternalCondition(eDesc));

		for (BasketEntry be : basket.entries) 
		{
			if(i++ % TRANSACTION_RESTART == 0)
				Globals.restartAllTransactions(true); // empirical to avoid automatic closing of session

			ConditionSet set = be.ep.conditions != null 
					? (ConditionSet) Globals.session().get(ConditionSet.class, be.ep.conditions.id) 
							: defaultConditions;

			dtConditions.addRow();
			for (ExternalCondition eDesc : eDescs) 
			{
				Property condition = eDesc.getProperty();

				PropertyValue pv = set.getValue(condition); // first try conditions of the given record
				if(pv == null) pv = defaultConditions.getValue(condition); 

				if(pv == null) throw new UserFriendlyException("Default condition for " + condition + " has not been provided" );

				double value  = 0;

				if (condition.isQualitative())
				{
					if(pv.option == null || model == null || 16888l == pv.option.id)
						System.out.println("Strange! pv.option= " + (pv.option == null?" null": pv.option.name) + " model = " + model);

					/*
					System.out.println("pv"+pv);
					System.out.println("pv.option"+pv.option);
					System.out.println("pv.option.id"+pv.option.id);
					System.out.println("model"+model);
					 */

					Long val = model.getMappedOption(pv.option.id);
					if(ProvidedConditions.isReplaceableCondition(condition.getName()) != null) {
						try{
							MixtureAttachment at = ExperimentalProperty.createSolventAttachment(pv.option.name);
							dtConditions.getCurrentRow().addAttachment(QSPRConstants.SOLVENT_ATTACHMENT, at); // do not need to be added
							// we also add SMILES for later use

							String smiles = 
									Various.molecule.convertToFormatFixMetal(data.getSDF(dtConditions.currentRow),QSPRConstants.SMILESH)  + "." + at.smiles();

							data.getRow(dtConditions.currentRow).addAttachment(QSPRConstants.SMILES_ATTACHMENT,smiles);

						}catch(Exception e) {
							data.getRow(dtConditions.currentRow).setError(e.getMessage());
						}
					}else
						if(val == null)
							data.getRow(dtConditions.currentRow).setError("Option " + pv.option.name +" for condition " + condition.getName() + " was not used during model training ");
						else
							value = (double) val;

					dtConditions.setColumnAttachment(condition.getName(), QSPRConstants.IS_QUALITATIVE_COLUMN, true);
					dtConditions.getCurrentRow().addAttachment(condition.getName().toLowerCase(), pv.option.name.toLowerCase());
				}
				else
					try{
						value = pv.value;
						if (!pv.unit.id.equals(eDesc.unitId))
							// Convert to the default unit, specified for this condition
							value = UnitConversion.convert(pv.value, pv.unit, eDesc.getUnit(), be.ep.molecule.molWeight);


						if(Math.abs(value) > QSPRConstants.MAXVALUE)
							data.getRow(dtConditions.currentRow).setError(("Absolute value of condition " + condition.getName() + " in " + eDesc.getUnit() + " unit is too large " + value + " and likely to impact performance of machine learning methods. Use another unit (e.g. log(" + eDesc.getUnit() +") scale or check your data"));

						if(value != 0 && Math.abs(value) < 1./QSPRConstants.MAXVALUE)
							data.getRow(dtConditions.currentRow).setError(("Absolute value of condition  " + condition.getName() + " in " + eDesc.getUnit() + " unit is too small " + value + " and likely to impact performance of machine learning methods. Use another unit (e.g. log(" + eDesc.getUnit() +") scale or check your data"));

					}catch(UserFriendlyException e) {
						throw new UserFriendlyException(e.getMessage() +  " for record R" +  pv.id);
					}

				dtConditions.setValue(condition.getName(), value);
				dtConditions.setColumnAttachment(condition.getName(), QSPRConstants.IS_CONDITION_COLUMN, true);
			}
		}

		if(data.getRowsNoErrorsSize() == 0 && data.getRowsSize() > 10)
			throw new UserFriendlyException("All records failed, example error: " + data.getRow(0).detailedStatus);


		Globals.restartAllTransactions(true);

	}

	@SuppressWarnings("unchecked")
	void enumerateConditions(List<ShortCondition> externalDescriptors) {

		Map<Long, Long> optionsMapping = model.attachment.getObject().optionsMapping;

		setStatus("Enumerating options of qualitative conditions...");
		List<PropertyOption> options = Globals.session().createCriteria(ExperimentalProperty.class)
				.add(Restrictions.in("id", model.getFilteredSet(QSPRConstants.TRAINING).getIdentifiers()))
				.createAlias("conditions", "c")
				.createAlias("c.values", "v")
				.add(Restrictions.isNotNull("v.option"))
				.setProjection(Projections.groupProperty("v.option")).list();

		for (ShortCondition eDesc : externalDescriptors)
		{
			eDesc.calculateMergings();
			// Consequently enumerate all the appearing options for this condition (eDesc)
			long num = 0;
			Property condition = (Property) Globals.session().get(Property.class, eDesc.id);
			if (!condition.isQualitative())
				continue;
			for (PropertyOption option : options) 
				if (option.property.equals(condition))
				{
					long sameOptionId = eDesc.getAnalogueOption(option.id);
					if (optionsMapping.containsKey(sameOptionId))
						optionsMapping.put(option.id, optionsMapping.get(sameOptionId));
					else
					{
						optionsMapping.put(sameOptionId, num);
						optionsMapping.put(option.id, num);
						num++;
					}
					setStatus(condition.getName() + "->" + option.name + " is mapped to " + optionsMapping.get(option.id) + " -- " +condition.id +":" + option.id);
				}

			if (eDesc.optionsMergings != null)
			{
				// Merge some of the options together
				for (ShortCondition.OptionsMerging om : eDesc.optionsMergings) {
					Long num1 = optionsMapping.get(om.id1);
					Long num2 = optionsMapping.get(om.id2);
					if (num1 == null || num2 == null)
						continue;
					long min = Math.min(num1, num2);
					if (num1 != min)
						optionsMapping.put(om.id1, min);
					if (num2 != min)
						optionsMapping.put(om.id2, min);
				}
			}
		}

	}

}


