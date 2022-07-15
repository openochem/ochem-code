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

package qspr.controllers;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.SessionVariable;
import qspr.dao.MetalBondParserSdf;
import qspr.dao.Repository;
import qspr.dao.Various;
import qspr.entities.Alert;
import qspr.entities.Basket;
import qspr.entities.Molecule;
import qspr.entities.PendingTask;
import qspr.export.ExportableColumn;
import qspr.export.ExportableSet;
import qspr.export.ExportableSetConfiguration;
import qspr.frontend.DescriptorsRow;
import qspr.frontend.DescriptorsSummary;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;
import qspr.metaserver.configurations.CorinaConfiguration;
import qspr.modelling.CompoundsProvider;
import qspr.modelling.DataPreprocessingParser;
import qspr.modelling.configurations.CDSConfiguration;
import qspr.modelling.configurators.DescriptorsConfigurator;
import qspr.util.AccessChecker;
import qspr.util.ExportThread;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;

import com.eadmet.business.DescriptorsCalculatorProcessor;
import com.eadmet.useractions.EventFactory;
import com.eadmet.utils.NumericalValueStandardizer;

@Controller
public class DescriptorsCalculatorController extends BrowserWrapper
{
	public DescriptorsCalculatorController()
	{
		sessionRequired = true;
	}

	public ModelAndView show(HttpServletRequest request, HttpServletResponse response)
	{
		AccessChecker.requestRegisteredUserPrivileges();
		return new WebModel().setTemplate("descriptorscalculator/dc-choose-compounds").getModelAndView();
	}

	public ModelAndView submitCompounds(HttpServletRequest req, HttpServletResponse res) throws Exception {
		AccessChecker.requestRegisteredUserPrivileges();
		DescriptorsCalculatorProcessor processor = new DescriptorsCalculatorProcessor();

		CompoundsProvider provider = new CompoundsProvider();
		final Basket basket = provider.parseUI().getBasket();
		processor.setDescription = provider.setDescription;

		if (basket == null)
			return new WebModel(new Alert(provider.error)).setTemplate("descriptorscalculator/choose-compounds").getModelAndView();

		processor.basket = basket;
		Globals.setSessionAttribute(SessionVariable.DESCRIPTORS_CALCULATOR_PROCESSOR, processor);

		return redirect("descriptorscalculator/dataPreprocessing.do");
	}

	public ModelAndView dataPreprocessing(HttpServletRequest request, HttpServletResponse response)
	{
		AccessChecker.requestRegisteredUserPrivileges();
		return new WebModel().setTemplate("descriptorscalculator/dc-data-preprocessing").getModelAndView();
	}

	public ModelAndView dataPreprocessingSubmit(HttpServletRequest request, HttpServletResponse response) throws NumberFormatException, Exception
	{
		AccessChecker.requestRegisteredUserPrivileges();
		DescriptorsCalculatorProcessor processor = getProcessor(request);
		DataPreprocessingParser.parseStandartizationUI(request, processor.structureStandartisation, null, null);
		return redirect("descriptorscalculator/selectDescriptors.do");
	}

	public ModelAndView selectDescriptors(HttpServletRequest request, HttpServletResponse response)
	{
		AccessChecker.requestRegisteredUserPrivileges();
		return new DescriptorsConfigurator().descriptorBlocks().setTemplate("descriptorscalculator/dc-descriptors").getModelAndView();
	}

	public ModelAndView descriptorsSubmit(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		AccessChecker.requestRegisteredUserPrivileges();
		DescriptorsCalculatorProcessor processor = getProcessor(request);
		DescriptorsConfigurator descConfigurator = new DescriptorsConfigurator();
		CDSConfiguration configuration = descConfigurator.configureDescriptors(request, new CDSConfiguration());
		processor.descConfig = configuration.descriptors;

		if (configuration.descriptors.requires3D())
		{
			processor.structureOptimisation = new CorinaConfiguration();
			return redirect("descriptorscalculator/selectOptimisation.do");
		}
		else
			return start(processor);
	}

	private ModelAndView start(DescriptorsCalculatorProcessor processor)
	{
		processor.start();
		processor.operation.successURL = "descriptorscalculator/results.do";
		Globals.setSessionAttribute(SessionVariable.DESCRIPTORS_CALCULATOR_PROCESSOR, processor);

		EventFactory.document("Descriptor calculation", null, String.format("has started calculation of descriptors for %d compounds", processor.basket.getRowsSize()));

		return redirect("longoperations/operationWaitingScreen.do?operation-id="
				+ processor.operation.operationId);
	}

	public ModelAndView selectOptimisation(HttpServletRequest request, HttpServletResponse response)
	{
		AccessChecker.requestRegisteredUserPrivileges();
		return new DescriptorsConfigurator().descriptorBlocks().setTemplate("descriptorscalculator/dc-structure-optimisation").addParam("recommendation", "1").getModelAndView();
	}

	public ModelAndView optimisationSubmit(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		AccessChecker.requestRegisteredUserPrivileges();
		DescriptorsCalculatorProcessor processor = getProcessor(request);
		processor.structureOptimisation = DescriptorsConfigurator.configureStructureOptimisation(request);
		return start(processor);
	}


	private DescriptorsRow formRow(DescriptorsCalculatorProcessor processor, DataTable descriptors, int index) throws IOException
	{
		DescriptorsRow dt = new DescriptorsRow();
		dt.id = index;
		AbstractDataRow row = descriptors.getRow(index);
		if (row.isError())
			dt.error = row.detailedStatus;

		Long molid = processor.basket.entries.get(index).ep.molecule.id;
		Molecule mol = Repository.molecule.getMolecule(molid);
		if (mol != null)
		{
			dt.name = Various.molecule.convertToCanonicalName(MetalBondParserSdf.eliminateMetalBond(mol.getData()));
			dt.mol = mol;
		}
		else
		{
			dt.mol = Molecule.getStub();
			dt.name = "Invalid molecule";
		}

		for (int i = 0; i < Math.min(descriptors.getColumnsSize(), 150); i++)
			if ((Double)row.getValue(i) != null)
				dt.descriptors.add(NumericalValueStandardizer.getSignificantDigits((Double)row.getValue(i)));
			else
				dt.descriptors.add("-");

		return dt;
	}

	private boolean matchesFilters(AbstractDataRow row)
	{
		if (assertParam("errors_only"))
			return (row.isError());
		else
			return true;
	}

	@Override
	public ModelAndView list(HttpServletRequest req, HttpServletResponse res) throws Exception 
	{
		AccessChecker.requestRegisteredUserPrivileges();
		DescriptorsCalculatorProcessor processor = getProcessor(req);
		DataTable descriptors = processor.dtDescriptors;
		int pagenum = getPageNum();
		int pagesize = getPageSize(15);
		List<Integer> index = new ArrayList<Integer>();

		for (AbstractDataRow row : descriptors)
			if (matchesFilters(row))
				index.add(descriptors.currentRow);

		List<Integer> subindex = index.subList(pagesize * (pagenum - 1), Math.min(pagesize * pagenum - 1, index.size()));
		List<Object> result = new ArrayList<Object>();

		for (Integer i : subindex) 
			result.add(formRow(processor, descriptors, i));

		WebList list = new WebList();
		list.list = result;
		list.size = index.size();
		list.pageNum = pagenum;
		list.pageSize = pagesize;

		return new WebModel(list).getModelAndView();
	}

	public ModelAndView results(HttpServletRequest request, HttpServletResponse response) throws NumberFormatException, Exception
	{
		AccessChecker.requestRegisteredUserPrivileges();
		if (request.getParameter("task") != null)
			Globals.setSessionAttribute(SessionVariable.DESCRIPTORS_CALCULATOR_PROCESSOR, null);
		DescriptorsCalculatorProcessor processor = getProcessor(request);
		DataTable descriptors = processor.dtDescriptors;
		DescriptorsSummary ds = new DescriptorsSummary();
		ds.columns.addAll(descriptors.getColumns().subList(0, Math.min(descriptors.getColumnsSize(), 150)));
		for (AbstractDataRow row : descriptors)
		{
			if (row.isError())
				ds.error++;
			else
				ds.valid++;
			ds.total++;
		}
		return new WebModel(ds).setTemplate("descriptorscalculator/dc-results").getModelAndView();
	}

	public ModelAndView download(HttpServletRequest request, HttpServletResponse response) throws NumberFormatException, Exception
	{
		AccessChecker.requestRegisteredUserPrivileges();
		final ExportableSet eData = new ExportableSet();
		eData.removeColumn(ExportableColumn.APPLICABILITY_DOMAIN); 
		eData.removeColumn(ExportableColumn.PREDICTED_VALUE); 
		eData.removeColumn(ExportableColumn.EXP_VALUE);
		eData.removeColumn(ExportableColumn.EXP_VALUE_CONVERTED); 
		eData.removeColumn(ExportableColumn.DM_VALUE);
		eData.removeColumn(ExportableColumn.ACCURACY);
		eData.removeColumn(ExportableColumn.DESCRIPTORSNAMES);

		DescriptorsCalculatorProcessor processor = getProcessor(request);

		if (assertParam("submit"))
		{
			final String format = getParam("format");
			ExportThread eThread = processor.getExportThread(format, ExportableSetConfiguration.configureFromDialog(request));
			eThread.start();
			Thread.sleep(100);
			return redirect("longoperations/operationWaitingScreen.do?operation-id="+eThread.operationID);
		}
		else
		{
			eData.selectedColumns.add(ExportableColumn.DESCRIPTORS);
			return new WebModel(eData).setTemplate("export").getModelAndView();
		}
	}

	private DescriptorsCalculatorProcessor getProcessor(HttpServletRequest req) throws NumberFormatException, Exception
	{
		DescriptorsCalculatorProcessor processor = (DescriptorsCalculatorProcessor) Globals.getSessionAttribute(SessionVariable.DESCRIPTORS_CALCULATOR_PROCESSOR);;
		if (req.getParameter("task") != null)
		{
			Long taskId = Long.valueOf(req.getParameter("task"));
			if (processor == null || processor.pTask == null || processor.pTask.id == null || !processor.pTask.id.equals(taskId))
			{
				processor = new DescriptorsCalculatorProcessor();
				processor.restoreFromPendingTask(PendingTask.getById(taskId));
			}
			Globals.setSessionAttribute(SessionVariable.DESCRIPTORS_CALCULATOR_PROCESSOR, processor);
		}

		return processor;

	}
}