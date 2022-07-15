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

package qspr.metaserver.cs;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

import com.eadmet.exceptions.CriticalException;

import qspr.dao.ChemInfEngine;
import qspr.dao.Various;
import qspr.metaserver.Event;
import qspr.metaserver.EventListener;
import qspr.metaserver.configurations.StandartizationOptions;
import qspr.metaserver.util.MoleculeStandartizer;
import qspr.workflow.utils.QSPRConstants
;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;

public class MolStandartizationServer extends WorkflowNodeServer
{
	private Map<String, MoleculeStandartizer> stCache = null; // only one should be used

	public MolStandartizationServer()
	{
		supportedTaskType = QSPRConstants.MolStandartizer;
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
	}

	@Override
	public WorkflowNodeData calculate(WorkflowNodeData task, Serializable configuration) throws Exception
	{

		DataTable dtMolecules = task.ports.get(0);
		StandartizationOptions stOptions = (StandartizationOptions) configuration;

		stCache = new HashMap<String, MoleculeStandartizer>();

		return new WorkflowNodeData(standartizeSDFMolecules(dtMolecules, stOptions));
	}

	private DataTable standartizeSDFMolecules(DataTable datatable, StandartizationOptions standartization) throws Exception
	{
		datatable.compressStrings = true;

		// Standartization of molecules
		// we can have different standartisers and different options for them
		// we create a list of them and 

		EventListener<String> statusListener = new EventListener<String>()
		{
			@SuppressWarnings("rawtypes")
			@Override
			public void onEvent(Event event, String arg)
			{
				setStatus(arg);
			}
		};

		if (standartization != null && standartization.standartizationRequired()){

			standartization.synchronise();

			// We initialize 
			Various.molecule = Various.getCheminfImpl(standartization.getDefault());

			if (standartization.desaltWith != null){
				getStandartizer(standartization.desaltWith).setDesalt();
				getStandartizer(standartization.desaltWith).statusChange.addListener(statusListener);
			}

			if (standartization.neutralizeWith != null){
				getStandartizer(standartization.neutralizeWith).setNeutralize();
				getStandartizer(standartization.neutralizeWith).statusChange.addListener(statusListener);
			}

			if (standartization.standardizeWith != null){
				getStandartizer(standartization.standardizeWith).setStandardize();
				getStandartizer(standartization.standardizeWith).statusChange.addListener(statusListener);
			}

			if (standartization.cleanStructureWith != null){
				getStandartizer(standartization.cleanStructureWith).setCleanStructure();
				getStandartizer(standartization.cleanStructureWith).statusChange.addListener(statusListener);
			}

			if (standartization.addExplicitHydrogensWith != null){
				getStandartizer(standartization.addExplicitHydrogensWith).addExplicitHydrogens = true;
				getStandartizer(standartization.addExplicitHydrogensWith).statusChange.addListener(statusListener);
			}

			if (standartization.dearomatizeWith != null){
				getStandartizer(standartization.dearomatizeWith).setDearomatize();
				getStandartizer(standartization.dearomatizeWith).statusChange.addListener(statusListener);
			}

			// doing different standartizations
			for (MoleculeStandartizer m : stCache.values()){
				if (standartization.outFormat != null)
					m.outFormat = standartization.outFormat;
				datatable = m.doStandartizationTable(datatable, this);
			}

			out.print("Standartization is finished");

		}

		if(Various.molecule != null && Various.molecule.engine != ChemInfEngine.CHEMAXON)
			datatable.setInfEngine(Various.molecule.engine);

		return datatable;
	}

	private MoleculeStandartizer getStandartizer(ChemInfEngine engine) throws Exception
	{
		String key = engine.toString();

		if (stCache.containsKey(key))
			return stCache.get(key);
		else
		{
			MoleculeStandartizer ms = MoleculeStandartizer.getInstance(engine,workingDirectory);
			stCache.put(key, ms);
			if(stCache.size() >1) throw new CriticalException("Only one Standartizer should be used");
			return ms;
		}
	}

}
