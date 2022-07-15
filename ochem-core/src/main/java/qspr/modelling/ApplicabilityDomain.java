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

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import qspr.Globals;
import qspr.entities.Model;
import qspr.entities.ModelMapping;
import qspr.metaserver.configurations.ADConfiguration;
import qspr.modelling.applier.ModelApplierTaskProcessor;
import qspr.modelling.configurations.CDSModelData;
import qspr.workflow.datatypes.DataTable;

@XmlRootElement(name = "ad")
public class ApplicabilityDomain 
{
	@XmlElement
	public boolean showNegativeValues = false;

	@XmlElement
	public boolean useStandardizedResiduals = false;

	@XmlElement(name = "model")
	public ModelMapping modelMapping;

	@XmlTransient
	private DataTable dtNewPredictions;

	@XmlAttribute(name = "name")
	public String dmName;

	//@XmlElementWrapper(name = "training-set")
	//@XmlElement(name = "point")
	//protected List<ADPoint> trainingSet = new ArrayList<ADPoint>();

	//@XmlElementWrapper(name = "test-set")
	//@XmlElement(name = "point")
	//protected List<ADPoint> testSet = new ArrayList<ADPoint>();

	@XmlElement(name = "ad-configuration")
	public ADConfiguration adConfigurationTraining;

	//@XmlElement(name = "ad-configuration-test")
	//public ADConfiguration adConfigurationTest;

	@XmlTransient
	public SetStatistics ssTraining;
	//@XmlTransient
	//protected SetStatistics errorSource;

	//	private void addPoint(double dm, double error, long epId, boolean isTrainingSet)
	//	{
	//		ADPoint point = new ADPoint();
	//		point.dm = dm;
	//		point.error = error;
	//		if (!showNegativeValues)
	//			point.error = Math.abs(point.error);
	//		point.epId = (int)epId;
	//		
	//		if (isTrainingSet)
	//			trainingSet.add(point);
	//		else
	//			testSet.add(point);
	//	}

	@XmlElement
	public Double dmThreshold;

	public int countOutsideAD = 0;

	@XmlElement
	private ModelStatistics ms;

	public double getPredictedError(double dm)
	{
		if (adConfigurationTraining == null || adConfigurationTraining.intervals == null)
			// No AD estimated for this model
			return 0; 

		return adConfigurationTraining.getAverageAccuracy(dm);
	}

	public Double getPredictedError(int i)
	{
		Double dm = getPredictedDM(i);
		if (dm == null)
			return null;
		return getPredictedError(dm);
	}

	public Double getPredictedDM(int i)
	{
		return (Double)dtNewPredictions.getValue(i, getDMColumnName());
	}

	private String getDMColumnName()
	{
		return (modelMapping.model.modelMappings.size() == 1 ?  "DM" : "DM" + modelMapping.getIndex()) + ":" + dmName;
	}

	public void setModel(ModelMapping sourceModel, String _name, ADConfiguration adConfiguration) throws Exception
	{
		this.modelMapping = sourceModel;
		ssTraining = ModelStatistics.get(sourceModel).sets.get(0);

		if (_name == null)
			dmName = ssTraining.distancesToModel.get(0);
		else
			dmName = _name;

		//adConfiguration = errorSource.adConfigurations.get(dmName);//targetModel.getAttachment().adConfiguration;
		int dmIndex = ssTraining.distancesToModel.indexOf(dmName);
		this.adConfigurationTraining = ssTraining.getADConfiguration(dmName, ssTraining, modelMapping, adConfiguration);
		dmThreshold = adConfigurationTraining.threshold;

		int descriptors = 1;

		if(sourceModel.model.calcDescriptors != null)
			descriptors = sourceModel.model.getCalculatedDescriptors().getColumnsSize() + 1;
		else
			if(sourceModel.model.attachment != null && sourceModel.model.attachment.getObject().configuration instanceof CDSModelData) {
				CDSModelData cds = (CDSModelData)sourceModel.model.attachment.getObject().configuration;
				if(cds.selectionConfiguration != null)descriptors = cds.selectionConfiguration.descriptorAsStrings().size()+1;
			}

		if (useStandardizedResiduals)
			dmThreshold = 3.0 * descriptors / (ssTraining.points.size());

		boolean considerPredicates = Globals.considerPredicates();
		for (SetStatistics ss : ModelStatistics.get(sourceModel).sets)
		{
			//SetStatistics dmSourceValidation = ((ModelStatistics)targetModel.statisticsRecalculated.getObject()).getSetStatistics(QSPRConstants.VALIDATION);

			if (!ss.distancesToModel.contains(dmName))
				continue;

			if (ss.getRowsSize() > 0)
				ss.adConfiguration = ss.getADConfiguration(dmName, ss, modelMapping, adConfiguration);

			//if (!targetModel.property.isQualitative())
			for (int i = 0; i < ss.points.size(); i++) 
			{
				PointStatistics ps = ss.points.get(i);
				if (ps.error == null)
				{
					double residual = ps.getRealValue(considerPredicates) - ps.predicted;
					double dmValue = ps.distancesToModel.get(dmIndex);
					if (!showNegativeValues)
						residual = Math.abs(residual);
					if (useStandardizedResiduals)
					{
						double variance = ssTraining.rmse * Math.sqrt(ssTraining.points.size() / (ssTraining.points.size() - descriptors - 1));
						residual = residual / (variance * Math.sqrt(1 - Math.min(1, dmValue)));
					}
					ps.setADData(dmValue, residual);
				}
			}

			if (useStandardizedResiduals)
				break; // only the training set is shown
		}
	}

	public void setTargetModelResults(ModelApplierTaskProcessor modelTask) throws Exception
	{
		countOutsideAD = 0;
		this.dtNewPredictions = modelTask.wndResult.ports.get(0);
		String dmColumn = getDMColumnName();
		if (!dtNewPredictions.containsColumn(dmColumn))
			return;


		ModelStatistics ms = ModelStatistics.get(modelMapping);
		SetStatistics ssPredicted = new SetStatistics().setId("predicted");
		ms.sets.add(ssPredicted);

		dtNewPredictions.reset();
		while (dtNewPredictions.nextRow())
		{
			if (!"error".equals(dtNewPredictions.getCurrentRow().status))
			{
				// Predicted error for test set
				Double dm = (Double)dtNewPredictions.getValue(dmColumn);
				PointStatistics ps = new PointStatistics();
				ps.setADData(dm, getPredictedError(dm));
				ssPredicted.points.add(ps);

				if (dmThreshold != null && dm > dmThreshold)
					countOutsideAD++;
				//addPoint(dm, new Double(getPredictedError(dm)), -1, false);
			}
		}

		this.ms = ms;

		if (this.ms.getRowsSize() > 20000)
			this.ms = null;

		Collections.sort(ssPredicted.points, new PointStatistics.DMComparator());
	}

	public static boolean hasDM(Model model) throws Exception
	{
		if(model.modelMappings == null || model.modelMappings.isEmpty())
			return false;

		ModelStatistics ms = ModelStatistics.get(model.modelMappings.get(0));
		if (ms == null) 
			return false;
		return ms.sets.get(0).distancesToModel.size() > 0;
	}

	public static String getDmName(Model model) throws Exception
	{
		return ModelStatistics.get(model.modelMappings.get(0)).sets.get(0).distancesToModel.get(0);
	}

}

//class ADPoint implements Comparable<ADPoint>
//{
//	@XmlAttribute
//	double error;
//	
//	@XmlAttribute
//	double dm;
//	
//	@XmlAttribute
//	int epId;
//	
//	@XmlAttribute
//	int molId;
//	
//	@XmlAttribute
//	Boolean dtDeleted;
//
//	public int compareTo(ADPoint o) 
//	{
//		return o.dm > dm ? -1 : o.dm == dm ? 0 : 1;
//	}
//}

class Histogram
{
	public List<HistogramBar> accuracy = new ArrayList<HistogramBar>();

	class HistogramBar
	{
		public double accuracy;
		public int frequency;
	}
}


