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

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.Globals;
import qspr.entities.Model;
import qspr.entities.ModelMapping;
import qspr.entities.Property;
import qspr.entities.PropertyValue;
import qspr.export.ExportableMolecule.ExportableValue;
import qspr.export.ExportableMolecule.Prediction;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.MemoryUtils;
import com.eadmet.utils.NumericalValueStandardizer;

/**
 * An abstract writer for exportable data (molecules, predictions, experimental records, etc)
 * 
 * @author midnighter
 *
 */

abstract public class ExportWriter
{
	private static transient final Logger logger = LogManager.getLogger(ExportWriter.class);

	abstract public void writeSupplementaryData(String key, Object value);
	@SuppressWarnings("rawtypes")
	abstract public void writeRow(List row);

	abstract public void initialize() throws IOException;

	abstract public String getHttpContentType();

	protected ExportableSet data;
	public OutputStream os;
	protected List<String> columns = new ArrayList<String>();
	protected String fileName = "exported_data";
	protected String directory = "";
	protected ExportableMolecule currentMolecule;

	protected String format;
	//	private int restrictedAccessCount = 0;
	private boolean restrictedAccessExport = false;
	//	private int restrictedAccessClearance = 10;
	public Double exportCost = 0D;

	public final static String MISSED_VALUE = "-";

	public boolean skipDescriptor(Object desc){
		return false;
	}

	public void write() throws Exception
	{
		//		restrictedAccessClearance = AccessChecker.getWeeklyExportLimit(Globals.userSession().user) - ExportAction.getWeeklyExport(Globals.userSession());
		//		logger.info("Restricted records clearance: " + restrictedAccessClearance + " restricted records allowed");

		data.convertUnits();

		if (data.shuffleKey != null && data.getDescriptors() != null)
		{
			if (data.shuffleKey.name == null)
				data.shuffleKey.name = data.getDescriptors().getRowsSize() + " descriptors";
			data.shuffleKey.setTargetDescriptors(data.getDescriptors());
			data.shuffleKey.save();
			fileName += "_ShuffleKey_" + data.shuffleKey.id;
		}

		boolean processed[] = null;

		if(data.selectedColumns.contains(ExportableColumn.COMPRESSED)) {
			processed = new boolean[data.exportableMolecules.size()];
			data.selectedColumns.remove(ExportableColumn.COMPRESSED);
			data.selectedColumns.remove(ExportableColumn.INTRODUCER);
			data.selectedColumns.remove(ExportableColumn.MODIFIER);
			data.selectedColumns.remove(ExportableColumn.N);
			data.selectedColumns.remove(ExportableColumn.EXTERNAL_ID);
			data.selectedColumns.remove(ExportableColumn.APPLICABILITY_DOMAIN);
		}

		// Define the column names
		createColumns();
		initialize();

		// Export the data
		for (int i = 0; i < data.exportableMolecules.size(); i++)
			if(processed==null || !processed[i])
			{
				currentMolecule = data.exportableMolecules.get(i);
				List<Object> row = processMolecule(currentMolecule);
				int mol = currentMolecule.moleculeId == null ? -1 : currentMolecule.moleculeId;

				for (int j = i + 1; j < data.exportableMolecules.size(); j++) 
					if(processed !=null && !processed[j])
					{
						if(mol == (int)data.exportableMolecules.get(j).moleculeId) {
							processed[j] = true;
							List<Object> rownew =  processMolecule(data.exportableMolecules.get(j));
							for(int k = 0; k < row.size(); k++) {
								if(row.get(k) ==  null)
									row.set(k, rownew.get(k));
							}
						}
					}

				if(processed !=null)currentMolecule.sheet = 0; // just use one sheet for xls

				writeRow(row);

				if (i % 100 == 0)
					logger.info("["+Globals.now()+"] Writing the exported data: Processed "+i+" records out of " + data.exportableMolecules.size());

				if (MemoryUtils.getCurrentMemoryUsedFraction() > 0.95)
					throw new UserFriendlyException("Unfortunately, the dataset you're exporting is too large and can not be handled now. Try exporting a smaller set or try again later.");

				if (i % 500 == 0)
					Globals.restartAllTransactions(true);
			}

		// Supplementary data
		for (String key : data.supplementaryData.keySet()) 
			writeSupplementaryData(key, data.supplementaryData.get(key));

		flush();

		// After everything is finished with no errors, register the fact of the data export
	}

	private List<Object> processMolecule(ExportableMolecule eMol)
	{
		// Create a row based on the selected columns
		List<Object> row = new ArrayList<Object>();
		for (ExportableColumn column : data.selectedColumns)
		{
			switch (column)
			{
			case ARTICLE:
				if (eMol.article != null)
				{
					row.add("A" + eMol.article.id);
					row.add(eMol.article.pmid);
					row.add(eMol.artPageNum);
					row.add(eMol.artTableNum);
				}
				else
				{
					row.add(MISSED_VALUE);
					row.add(MISSED_VALUE);
					row.add(MISSED_VALUE);
					row.add(MISSED_VALUE);
				}
				break;
			case SMILES:
				row.add(eMol.smiles);
				break;
			case APPLICABILITY_DOMAIN:
				row.add(eMol.ousideOfAD?"FALSE":"TRUE");
				break;
			case COMMENTS:
				row.add(eMol.comments);
				break;
			case INCHI_KEY:
				row.add(eMol.inchiKey);
				break;
			case NAMES:
				row.add(eMol.name);
				row.add(eMol.name2);
				break;
			case CASRN:
				row.add(eMol.casRN);
				break;
			case ERROR:
				String error = eMol.error;
				if (error == null && eMol.descriptors != null && eMol.descriptors.isError())
					error = eMol.descriptors.detailedStatus;
				row.add(error);
				break;
			case PREDICTED_VALUE:
				for (Model model : data.models)
					for (ModelMapping mm : model.modelMappings)
					{
						Prediction pred = eMol.getPrediction(mm);
						row.add(pred == null ? null : pred.predictedValue);
						if (mm.property.isQualitative())
							row.add(pred == null ? null : pred.numericPredictedValue);
					}
				break;
			case EXP_VALUE:
				for (int p = 0; p < data.properties.size(); p++)
				{
					Property prop = data.properties.get(p);
					ExportableValue v = eMol.expValues.get(prop);
					if (v != null)
					{
						row.add(v.getPredicatedValue());
						if (prop.isNumeric())
							row.add(v.unit == null ? MISSED_VALUE : v.unit.getName());
					}
					else
					{
						row.add(null);
						if (prop.isNumeric()) // units for numeric values
							row.add(null);
					}
				}
				break;
			case EXP_VALUE_CONVERTED:
				for (int p = 0; p < data.properties.size(); p++)
				{
					Property prop = data.properties.get(p);
					if (prop.isNumeric()){

						ExportableValue v = eMol.expValuesConverted.get(prop);
						if (v != null)
						{
							row.add(v.getPredicatedValue());
							if (prop.isNumeric())
								row.add(v.unit == null ? MISSED_VALUE : v.unit.getName());
						}
						else
						{
							row.add(null);
							if (prop.isNumeric()) // units for numeric values
								row.add(null);
						}
					}
				}
				break;					
			case DESCRIPTORS:
				if (data.getDescriptors() == null)
					break;
				int numDescriptors = data.shuffleKey != null ? data.shuffleKey.numDescriptors : data.getDescriptors().getColumnsSize();
				for (int p = 0; p < numDescriptors; p++)
					if (eMol.descriptors == null || eMol.descriptors.isError())
						row.add(MISSED_VALUE);
					else{

						Object desc = data.shuffleKey == null 
								?eMol.descriptors.getValue(p)
										:data.shuffleKey.getDescriptor(eMol.descriptors, p);

								if(skipDescriptor(desc))
									row.add(null);
								else{
									if (desc != null && desc instanceof Double)
										desc = NumericalValueStandardizer.getSignificantDigitsDouble((Double) desc,  NumericalValueStandardizer.SIGNIFICANT_DIGITS);
									row.add(desc == null ? "" : "" + desc);
								}
					}
				break;
			case RECORDID:
				row.add(eMol.recordId != null ? "R" + eMol.recordId : "");
				break;
			case INTRODUCER:
				row.add(eMol.introducer);
				break;
			case MODIFIER:
				row.add(eMol.modifier);
				break;
			case CONDITIONS:
				for (Property condition : data.conditions)
				{
					PropertyValue pv = eMol.conditions.get(condition);
					row.add(pv == null ? "" :  pv.getStringValue());
					if (condition.isNumeric())
						row.add(pv == null ? MISSED_VALUE : pv.unit.getName());
				}
				break;
			case DM_VALUE:
				for (ModelMapping mm : data.dmNames.keySet())
					for (String dm : data.dmNames.get(mm))
					{
						Prediction pred = eMol.getPrediction(mm);
						if (pred != null && pred.dms != null && pred.dms.get(dm) != null)
							row.add(NumericalValueStandardizer.getSignificantDigitsDouble(pred.dms.get(dm),NumericalValueStandardizer.SIGNIFICANT_DIGITS));
						else
							row.add(null);
					}
				break;
			case MOLECULEID:
				if (eMol.moleculeId != null)
					row.add("M" + eMol.moleculeId);
				else
					row.add(null);
				break;
			case ACCURACY:
				for (Model model : data.models)
					for (ModelMapping mm : model.modelMappings)
						if (data.dmNames.get(mm) != null)
						{
							Prediction pred = eMol.getPrediction(mm);
							row.add(pred == null || pred.accuracy == null? null : NumericalValueStandardizer.getSignificantDigitsDouble(pred.accuracy,NumericalValueStandardizer.SIGNIFICANT_DIGITS));
						}
				break;
			case EXTERNAL_ID:
				row.add(eMol.inhouseRecordId);
				break;
			case N:
				row.add(eMol.molNum);
				break;
			default:
				row.add(MISSED_VALUE);
			}

		}

		if (!data.selectedColumns.contains(ExportableColumn.SMILES))
			currentMolecule.sdf = "1\n2\n3\n 0  0  0  0  0  0            999 V2000\nM  END"; // empty molecule

		if (!eMol.ownRecord && !eMol.freelyAvailable)
			if (restrictedAccessExport)
			{
				//				restrictedAccessCount++;
				exportCost += getRecordCost(eMol);
			}

		return row;
	}

	protected void createColumns()
	{
		restrictedAccessExport = false;
		columns.clear();
		
		boolean predictedExist = data.selectedColumns.contains(ExportableColumn.PREDICTED_VALUE);
			
		
		for (ExportableColumn column : data.selectedColumns)
		{
			if (column.isRestricted())
				restrictedAccessExport = true;
			switch (column)
			{
			case ARTICLE:
				// names should be the same as those used for batch upload 
				columns.add("ARTICLEID");
				columns.add("PUBMEDID");
				columns.add("PAGE");
				columns.add("TABLE");
				break;
			case NAMES:
				columns.add("NAME");
				columns.add("NAME");
				break;
			case EXP_VALUE:
				for (Property property : data.properties)
				{
					columns.add(property.getName() + (predictedExist? " {measured}":""));
					if (property.isNumeric())
						columns.add("UNIT {"+property.getName()+"}");
				}
				break;
			case PREDICTED_VALUE:
				for (Model model : data.models)
					for (ModelMapping mm : model.modelMappings)
					{
						String unitSuffix = "";
						if (mm.property.isNumeric())
							unitSuffix = " in " + (data.unitsOfExport.containsKey(mm.property) ? data.unitsOfExport.get(mm.property) : mm.unit);
						columns.add(mm.property + " {predicted by ochem.eu/model/" + model.publicId + unitSuffix + "}");
						if (mm.property.isQualitative())
							columns.add("Numeric prediction for " + mm.property + " {predicted by ochem.eu/model/" + model.publicId + unitSuffix + "}");
					}
				break;
			case EXP_VALUE_CONVERTED:
				for (Property property : data.properties)
					if (property.isNumeric())
					{
						columns.add(property.getName() + " {measured, converted}");
						columns.add("UNIT {"+property.getName()+"}");
					}
				break;					
			case DESCRIPTORS:
			case DESCRIPTORSNAMES:
				if (data.getDescriptors() == null)
					break;
				int numDescriptors = data.shuffleKey != null ? data.shuffleKey.numDescriptors : data.getDescriptors().getColumnsSize();
				for (int i = 0; i < numDescriptors; i++)
					if (data.shuffleKey != null)
						columns.add("SHUFFLED_DESC_" + (i + 1));
					else
						columns.add(data.getDescriptors().getColumn(i));
				break;
			case CONDITIONS:
				for (Property condition : data.conditions)
				{
					columns.add(condition.getName());
					if (condition.isNumeric())
						columns.add("UNIT {"+condition.getName()+"}");
				}
				break;
			case DM_VALUE:
				for (ModelMapping mm : data.dmNames.keySet())
					for (String dm : data.dmNames.get(mm))
						columns.add(dm + "{" + mm.property + " by " + mm.model.name + "}");
				break;
			case ACCURACY:
				for (Model model : data.models)
					for (ModelMapping mm : model.modelMappings)
						if (data.dmNames.get(mm) != null)
							columns.add((mm.property.isQualitative() ? "Estimated accuracy" : "Estimated RMSE") + "{"+mm.property +" predicted by "+model.name+"}");
				break;
			case EXTERNAL_ID:
				columns.add(QSPRConstants.EXTERNALID);
				break;
			default:
				columns.add(column.name());
			}
		}

		if ("sdf".equals(format))
			restrictedAccessExport = true;
	}

	public void flush() throws IOException
	{
		os.flush();
		os.close();
	}

	public void setFileName(String filename)
	{
		fileName = filename;
	}

	public String getFileName()
	{
		return fileName;
	}

	public static ExportWriter createWriter(String format, ExportableSet data, OutputStream os)
	{
		ExportWriter writer = ExportWriter.getWriterByFormat(format);

		writer.format = format;
		writer.data = data;
		writer.os = os;

		return writer;
	}

	public static ExportWriter createWriter(String format, ExportableSet data, String directory, String fileNameWithoutExtension) throws FileNotFoundException
	{
		ExportWriter writer = ExportWriter.getWriterByFormat(format);

		writer.format = format;
		writer.data = data;

		writer.setFileName(fileNameWithoutExtension + "." + writer.getFileExtension());
		writer.os = new BufferedOutputStream(new FileOutputStream(new File(directory + "/" + writer.fileName)));

		return writer;
	}

	protected static ExportWriter getWriterByFormat(String format)
	{
		if (QSPRConstants.EXCEL.equals(format))
			return new ExcelExportWriter();
		else if ("csv".equals(format))
			return new CSVExportWriter();
		else if ("sdf".equals(format))
			return new SDFExportWriter();
		else if ("r".equals(format))
			return new RExportWriter();
		else
			throw new RuntimeException("Unsupported format: " + format);
	}

	public String getFileExtension()
	{
		return format; // can be overrriden. For example, "sdf" will have "sdf.gz" extension
	}

	public String getFullFilePath()
	{
		return directory + "/" + fileName;
	}

	@SuppressWarnings({ "rawtypes", "unchecked" })
	public void writeRow(Object[] row)
	{
		List list = new ArrayList<Object>();
		for (Object object : row)
			list.add(object);
		writeRow(list);
	}

	private Double getRecordCost(ExportableMolecule eMol)
	{
		return 1D;
	}

}
