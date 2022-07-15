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

package qspr.business.toxalert;

import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.dao.Various;
import qspr.entities.Basket;
import qspr.entities.BasketEntry;
import qspr.entities.PendingTask.TaskType;
import qspr.entities.SubstructureAlert;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.metaserver.configurations.StandartizationOptions;
import qspr.metaserver.configurations.DescriptorsStructuralAlertsConfiguration;
import qspr.metaserver.protocol.Task;
import qspr.metaserver.protocol.Task.TaskPriority;
import qspr.metaserver.util.ExtendedSMART;
import qspr.workflow.utils.QSPRConstants;
import qspr.modelling.AbstractTaskProcessor;
import qspr.modelling.CompoundsProvider;
import qspr.modelling.ModelApplierAttachment;
import qspr.util.Operation;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.NodeConfiguration;
import qspr.workflow.datatypes.NodesConfiguration;
import qspr.workflow.datatypes.WorkflowConfiguration;
import qspr.workflow.datatypes.WorkflowNodeData;

/**
 * Screens chemicals compounds against a set of structural alerts.
 * Posts a task to the Metaserver and registers a pending task on OCHEM side.
 * 
 * @author midnighter
 *
 */
public class ScreeningProcessor extends AbstractTaskProcessor 
{
	private static transient final Logger logger = LogManager.getLogger(ScreeningProcessor.class);

	public CompoundsProvider compoundsProvider = new CompoundsProvider();
	//public Basket basket;
	public Index<Integer, Long> alertsByCompounds = new Index<Integer, Long>();
	public Index<Long, Integer> compoundsByAlerts = new Index<Long, Integer>();
	public Index<Long, Integer> compoundsByEndpoints = new Index<Long, Integer>();
	public Index<Long, Integer> compoundsByPublications = new Index<Long, Integer>();
	public List<Integer> orderedCompounds = new ArrayList<Integer>();
	public List<SubstructureAlert> alerts;
	public AlertsFilter alertsFilter;

	public StandartizationOptions standartization = new StandartizationOptions();

	public ScreeningProcessor() {
		setOperationID(Operation.generateID());
		taskClass = TaskType.TOXALERT_SCREENING;
	}

	public static List<String> getSMARTs(List<SubstructureAlert> alerts)
	{
		List<String> alertPatterns = new ArrayList<String>();
		Iterator<SubstructureAlert> iAlerts = alerts.iterator();
		while (iAlerts.hasNext())
		{
			SubstructureAlert alert = iAlerts.next();
			if ((ExtendedSMART.create(alert.getFullSMARTS(), OCHEMConfiguration.getCheminfEngine())).invalid)
			{
				iAlerts.remove();
				logger.info("WARNING: Not a valid SMART - " + alert.smart);
			}
			else
				alertPatterns.add(alert.smart);
		} 

		return alertPatterns;
	}

	public int getTaskPriority()
	{
		int priority = compoundsProvider.getBasket().entries.size() > 10000 ? TaskPriority.NORMAL : TaskPriority.HIGH;
		return Math.min(priority, defaultTaskPriority);
	}

	@Override
	public void onTaskPosted() {
		compoundsProvider = null;
		System.gc();
	}

	public void onTaskReceived(Task task) throws IOException, ClassNotFoundException
	{
		compoundsProvider = new CompoundsProvider();
		compoundsProvider.setBasket(((ScreeningAttachment)pTask.attachment.getObject()).getWorkData(false));

		setStatus("Screening results received - building alert indices");
		logger.info("Screening results received - building alert indices");
		DataTable dtMatches = WorkflowNodeData.fromTask(task).ports.get(0);
		int[] bits = new int[(alerts.size() - 1)/32 + 1];
		for (int i = 0; i < compoundsProvider.basket.getRowsSize(); i++)
		{
			Integer mappingId2 = null;

			if (!dtMatches.getRow(i).isError()) {

				mappingId2 = compoundsProvider.basket.entries.get(i).ep.molecule.mapping2.id;

				for (int k = 0; k < dtMatches.getColumnsSize(); k++)
					bits[k] = Float.floatToIntBits(((Double)dtMatches.getValue(i, k)).floatValue());
				for (int k = 0; k < alerts.size(); k++)
				{
					if (bits[k / 32] == (bits[k / 32] | (1 << (k % 32))))
					{
						int cmpIndex = alertsByCompounds.addValue(mappingId2, alerts.get(k).id);
						compoundsByAlerts.addValue(alerts.get(k).id, cmpIndex);
						compoundsByEndpoints.addValue(alerts.get(k).property.id, cmpIndex);
						compoundsByPublications.addValue(alerts.get(k).article.id, cmpIndex);
					}
				}
			}

			orderedCompounds.add(mappingId2);

		}

		// Release memory
		compoundsProvider = null;
	}

	private ScreeningAttachment attachment;

	@Override
	public ScreeningAttachment getAttachment()
	{
		if (attachment != null)
			return attachment;

		attachment = new ScreeningAttachment();
		attachment.setWorkData(compoundsProvider.getBasket());
		attachment.setAlerts(alerts);

		return attachment;
	}

	@Override
	protected String getTaskType()
	{
		return QSPRConstants.Workflow;
	}

	@SuppressWarnings("unchecked")
	private DescriptorsConfiguration getAlertsConfiguration() throws Exception
	{
		if (alerts == null)
		{
			Criteria c = Globals.session().createCriteria(SubstructureAlert.class);
			alertsFilter.filterCriteria(c);
			alerts = c.list();
		}

		setStatus("Preparing SMARTS");
		DescriptorsStructuralAlertsConfiguration conf = new DescriptorsStructuralAlertsConfiguration();
		conf.compactMode = true;
		conf.alertPatterns = getSMARTs(alerts);

		DescriptorsConfiguration descConf = new DescriptorsConfiguration();
		descConf.addDescriptorType(DescriptorsConfiguration.StructuralAlerts, conf).skipCache();

		return descConf;
	}

	@Override
	protected Serializable getTaskConfiguration() throws Exception
	{
		NodesConfiguration nodesConfiguration = new NodesConfiguration();
		nodesConfiguration.nodes.add(new NodeConfiguration("mol-standardizer", standartization));
		nodesConfiguration.nodes.add(new NodeConfiguration("mol-optimiser", null).setSkipNode(true));
		nodesConfiguration.nodes.add(new NodeConfiguration("descriptors-processor",  getAlertsConfiguration()));
		nodesConfiguration.nodes.add(new NodeConfiguration("descriptors-selector", null).setSkipNode(true));
		return new WorkflowConfiguration("CDS", nodesConfiguration);
	}

	@Override
	protected Serializable getTaskData() throws Exception
	{
		compoundsProvider.status.addListener(statusTracker);

		Basket basket = compoundsProvider.getBasket();

		Globals.session().evict(basket);
		for (BasketEntry be : basket.entries)
			if (Globals.session().contains(be))
				Globals.session().evict(be);

		DataTable dtMolecules = new DataTable();
		dtMolecules.addColumn(QSPRConstants.SDF_COLUMN);
		dtMolecules.setInfEngine(Various.molecule.engine); 
		Iterator<BasketEntry> iEntries = basket.entries.iterator();
		while (iEntries.hasNext())
		{
			BasketEntry entry = iEntries.next();
			dtMolecules.addRow();
			dtMolecules.setValue(entry.ep.molecule.error != null? "": entry.ep.molecule.getData());
		}

		if(dtMolecules.getRowsNoErrorsSize() == 0)
			throw new IOException("No molecules were provided. Please, review your data.");

		return new WorkflowNodeData(dtMolecules);
	}

	@Override
	protected void restoreFromAttachment(Serializable uncastedAttachment)
	{
		ScreeningAttachment attachment = (ScreeningAttachment) uncastedAttachment;
		setStatus("Unpacking attachment");
		logger.info("Alerts...");
		alerts = attachment.getAlerts();
		logger.info("Molecules...");
		compoundsProvider = new CompoundsProvider();
		compoundsProvider.setBasket(attachment.getWorkData(false));
	}

	public static void main(String[] args)
	{
		int v = Float.floatToIntBits(Float.intBitsToFloat(1 << 31));
		System.out.println(v == (v | (1 << 31)));
	}
}

class ScreeningAttachment extends ModelApplierAttachment
{
	private static final long serialVersionUID = 1L;
	Long[] alertIDs;

	public void setAlerts(List<SubstructureAlert> alerts)
	{
		alertIDs = new Long[alerts.size()];
		for (int i = 0; i < alerts.size(); i++)
			alertIDs[i] = alerts.get(i).id;
	}

	public List<SubstructureAlert> getAlerts()
	{
		Globals.session().createCriteria(SubstructureAlert.class).add(Restrictions.in("id", alertIDs)).list();
		List<SubstructureAlert> alerts = new ArrayList<SubstructureAlert>();
		for (long alertID : alertIDs)
			alerts.add(SubstructureAlert.getByID(alertID));

		return alerts;
	}
}


