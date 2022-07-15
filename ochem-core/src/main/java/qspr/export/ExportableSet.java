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

package qspr.export;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlElementWrapper;
import javax.xml.bind.annotation.XmlRootElement;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.entities.Model;
import qspr.entities.ModelMapping;
import qspr.entities.Property;
import qspr.entities.ShuffleKey;
import qspr.entities.Unit;
import qspr.export.ExportableMolecule.ExportableValue;
import qspr.export.ExportableMolecule.Prediction;
import qspr.util.unitconversion.UnitConversion;
import qspr.workflow.datatypes.DataTable;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.MemoryUtils;

/**
 * The universal class for any data exported from OCHEM
 * Basically, its a set of ExportableMolecules with some meta information.
 * @author midnighter
 *
 */

@XmlRootElement(name = "exportable-data")
public class ExportableSet
{
	private static transient final Logger logger = LogManager.getLogger(ExportableSet.class);

	@XmlElement
	public List<ExportableColumn> selectedColumns = new ArrayList<ExportableColumn>(); // Columns selected by user in a dialog

	//Do we _really_ need a whole model here?
	public transient List<Model> models = new ArrayList<Model>();
	public transient Map<ModelMapping, List<String>> dmNames = new HashMap<ModelMapping, List<String>>();
	public transient List<Property> conditions = new ArrayList<Property>();
	public List<Property> properties = new ArrayList<Property>();
	public transient List<ExportableMolecule> exportableMolecules = new ArrayList<ExportableMolecule>();
	private DataTable dtDescriptors; // TODO Refactoring - how about descriptors for different models?
	public ShuffleKey shuffleKey;
	public Map<Property, Unit> unitsOfExport = new HashMap<Property, Unit>();
	public Map<String, Object> supplementaryData = new HashMap<String, Object>(); // Sheet name -> Object (typically String, but can be more universal, e.g. another Map)

	// Sheets are Excel-specific, but what the hack (maybe, in other formats we can also use this information)
	protected transient List<String> sheets = new ArrayList<String>();
	private int currentSheet = 0;

	// Information for UI
	@XmlElement
	private double freeBonuses;

	@XmlElement
	private double earnedBonuses;

	@XmlElement
	private List<UIColumn> availableColumns = new ArrayList<UIColumn>(); // All available columns for this export data
	@XmlElementWrapper(name="available-shuffle-keys")@XmlElement(name="shuffle-key")
	private List<ShuffleKey> availableShuffleKeys = new ArrayList<ShuffleKey>();

	public void selectAll()
	{
		selectedColumns.clear();
		for (UIColumn uiColumn : availableColumns)
			selectedColumns.add(uiColumn.column);
	}

	public void configure(ExportableSetConfiguration config)
	{
		// Selected columns
		if (config.selectAll)
			selectAll();
		else if (!config.useDefaults)
		{
			selectedColumns.clear();
			for (UIColumn uiColumn : availableColumns)
			{
				if (config.potentialColumnNames.contains((uiColumn.column.name())))
					selectedColumns.add(uiColumn.column);
			}
		}

		// Shuffle key for descriptors
		if (config.useShuffleKey)
		{
			if ("new".equals(config.shuffleKey))
				shuffleKey = new ShuffleKey();
			else
				shuffleKey = (ShuffleKey) Globals.session().get(ShuffleKey.class, Long.valueOf(config.shuffleKey));
		}

		// Units of export
		for (Map.Entry<Long,Long> entry : config.propertyToUnitMap.entrySet())
		{
			Property p = (Property) Globals.session().get(Property.class, entry.getKey());
			Unit u = (Unit) Globals.session().get(Unit.class, entry.getValue());
			unitsOfExport.put(p, u);
		}
	}

	public void convertUnits()
	{
		if (unitsOfExport.isEmpty())
			return;

		logger.info("Converting units..." + unitsOfExport);
		for (ExportableMolecule eMol : exportableMolecules)
		{
			// Convert predictions
			for (Model m : models)
				for (ModelMapping mm : m.modelMappings)
					if (eMol.predictedValues.keySet().contains(mm.id))
						if (unitsOfExport.containsKey(mm.property))
						{
							Prediction prediction = eMol.predictedValues.get(mm.id);
							Object value = (Double) prediction.predictedValue;
							if (value != null && value instanceof Double)
								prediction.predictedValue = UnitConversion.convert((Double)value, mm.unit, unitsOfExport.get(mm.property), eMol.molWeight);
						}

			// Convert experimental values
			for (Property p : eMol.expValuesConverted.keySet())
			{
				Unit u = unitsOfExport.get(p);
				if(u == null){
					if(p.parent != null)u = unitsOfExport.get(p.parent);
				}

				ExportableValue v = eMol.expValuesConverted.get(p);

				if (u == null || v == null)
					continue;

				if (v.predicate != null && (v.predicate.equals("-") || v.predicate.equals("\u2014")) 
						&& v.value != null && v.secondValue != null){ // use average value
					v.value = (v.value + v.secondValue)/2;
				}

				v.predicate = null; // we do not export predicates, intervals, etc. for converted units
				v.secondValue = null;

				try {
					if (v.value != null)
						v.value =  UnitConversion.convert(v.value, v.unit, u, eMol.molWeight);
				}catch(UserFriendlyException e) {
					eMol.error = e.getMessage();
				}

				v.unit = u;
			}
		}
		logger.info("Units have been converted");
	}

	public void addMolecule(ExportableMolecule m)
	{
		exportableMolecules.add(m);
		m.parent = this;
		m.sheet = currentSheet;
		if (exportableMolecules.size() % 5000 == 0)
			logger.info("["+Globals.now()+"] Preparing the data to export: Prepared " + exportableMolecules.size() + " records");

		if (MemoryUtils.getCurrentMemoryUsedFraction() > 0.95)
			throw new UserFriendlyException("Unfortunately, the dataset you're exporting is too large and can not be handled now. Try exporting a smaller set or try again later.");
	}


	public void addModel(Model model)
	{
		if (!models.contains(model))
		{
			models.add(model);
			if (shuffleKey != null && shuffleKey.name == null)
				shuffleKey.name = model.name;
		}
	}

	public void setDescriptors(DataTable dtDescriptors)
	{
		this.dtDescriptors = dtDescriptors;

		// We are creating a new shuffle key
		if (shuffleKey != null && shuffleKey.id == null)
			shuffleKey.createKeyFrom(dtDescriptors);
	}

	public DataTable getDescriptors()
	{
		return dtDescriptors;
	}

	@SuppressWarnings("unchecked")
	public ExportableSet()
	{
		// Presume, all the columns are available. 
		// In each particular place of export, remove unavailable columns manually
		for (ExportableColumn c : ExportableColumn.values())
			availableColumns.add(new UIColumn(c));
		freeBonuses = 0;
		earnedBonuses = 0;
		// Load available shuffle keys
		Criteria cShuffleKeys = Globals.session().createCriteria(ShuffleKey.class);
		if (Globals.isGuestUser())
			cShuffleKeys.add(Restrictions.eq("session", Globals.userSession()));
		else
			cShuffleKeys.createCriteria("session").add(Restrictions.eq("user", Globals.userSession().user));
		cShuffleKeys.addOrder(Order.desc("id"));
		availableShuffleKeys = cShuffleKeys.list();
		uncheckColumn(ExportableColumn.COMPRESSED);
		uncheckColumn(ExportableColumn.INCHI_KEY);
		uncheckColumn(ExportableColumn.COMMENTS);
		uncheckColumn(ExportableColumn.MODIFIER);
		uncheckColumn(ExportableColumn.INTRODUCER);
		uncheckColumn(ExportableColumn.RECORDID);
		uncheckColumn(ExportableColumn.MOLECULEID);
	}

	public void removeColumn(ExportableColumn column)
	{
		Iterator<UIColumn> i = availableColumns.iterator();
		while (i.hasNext())
			if (i.next().column == column)
				i.remove();
	}

	public void clearColumns()
	{
		availableColumns.clear();
	}

	public void addColumn(ExportableColumn column)
	{
		availableColumns.add(new UIColumn(column));
	}

	public void uncheckColumn(ExportableColumn uncheckedColumn)
	{
		for (UIColumn uiCol : availableColumns)
			if (uiCol.column == uncheckedColumn)
				uiCol.checked = false;
	}

	public void useSheet(String title)
	{
		if (!sheets.contains(title))
			sheets.add(title);

		currentSheet = sheets.indexOf(title);
	}
}

class UIColumn
{
	public ExportableColumn column;
	public boolean checked = true;
	public boolean restricted = false;
	public String title;

	public UIColumn(ExportableColumn c)
	{
		this.column = c;
		this.title = c.toString();
		this.restricted = c.isRestricted();
	}

	public UIColumn()
	{

	}
}
