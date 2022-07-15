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
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.SessionVariable;
import qspr.business.SetComparisonProcessor;
import qspr.business.WebFilters;
import qspr.entities.Alert;
import qspr.entities.Basket;
import qspr.entities.Mapping2Filter;
import qspr.entities.Molecule;
import qspr.entities.PendingTask;
import qspr.export.CSVExportWriter;
import qspr.frontend.BrowserModel;
import qspr.frontend.SetCompareResult;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;
import qspr.modelling.CompoundsProvider;
import qspr.modelling.DataPreprocessingParser;
import qspr.modelling.configurations.CDSConfiguration;
import qspr.modelling.configurators.DescriptorsConfigurator;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;

import com.eadmet.useractions.EventFactory;

/**
 * A controller that supports the "SetCompare" : a comparison of two molecule sets based on their features (scaffolds, alerts, etc)
 * @author midnighter
 *
 */
@Controller
public class SetComparisonController extends BrowserWrapper
{
	public ModelAndView select(HttpServletRequest request, HttpServletResponse response) 
	{
		Globals.setSessionAttribute(SessionVariable.SETCOMPARE_PROCESSOR, new SetComparisonProcessor());
		return new WebModel().setTemplate("setcomparison/sc-select").getModelAndView();
	}

	public ModelAndView selectSubmit(HttpServletRequest request, HttpServletResponse response)
	{
		CompoundsProvider provider1 = new CompoundsProvider();
		provider1.scope = "set1";
		Basket basket1 = provider1.parseUI().getBasket();

		CompoundsProvider provider2 = new CompoundsProvider();
		provider2.scope = "set2";
		Basket basket2 = provider2.parseUI().getBasket();

		SetComparisonProcessor processor = getProcessor();
		processor.setDescription = provider1.setDescription + " vs. " + provider2.setDescription;
		processor.basket1 = basket1;
		processor.basket2 = basket2;

		return redirect("setcomparison/standardize.do");
	}

	public ModelAndView standardize(HttpServletRequest request, HttpServletResponse response) 
	{
		return new WebModel().setTemplate("setcomparison/sc-standardize").getModelAndView();
	}

	public ModelAndView exportResults(HttpServletRequest request, HttpServletResponse response) throws IOException 
	{
		SetComparisonProcessor processor = getProcessor();
		DataTable dtResult = processor.wndResult.ports.get(0);
		List<Object> results = new ArrayList<Object>();

		int pageNum = getPageNum();
		int pagesize = getPageSize(15);
		for (int i = (pageNum - 1) * pagesize; i < Math.min(dtResult.getRowsSize(), pageNum * pagesize); i++)
			results.add(new SetCompareResult(dtResult.getRow(i)).setId(i));

		OutputStream os = response.getOutputStream();
		response.setContentType("application/csv");
		response.setHeader("Content-Disposition", "attachment; filename=setcompare-results.csv");

		CSVExportWriter eWriter = new CSVExportWriter();
		eWriter.os = os;
		eWriter.pw = new PrintWriter(os);
		eWriter.initialize();
		eWriter.writeRow(new String[]{
				"Descriptor ID",
				"Alert ID",
				"In Set1",
				"In Set1 (%)",
				"In Set2",
				"In Set2 (%)",
				"pValue"
		});

		int n1 = processor.basket1.entries.size() - processor.getErrorIndices(0).size();
		int n2 = processor.basket2.entries.size() - processor.getErrorIndices(1).size();

		for (AbstractDataRow row : dtResult)
		{
			SetCompareResult rowResult = new SetCompareResult(row);
			eWriter.writeRow(new Object[]{
					rowResult.displayName, 
					rowResult.alertId != null ? "TA" + rowResult.alertId : "",
							rowResult.inSet1, 100.0 * rowResult.inSet1 / n1,
							rowResult.inSet2, 100.0 * rowResult.inSet2 / n2,
							rowResult.pValue1
			});
		}

		eWriter.flush();

		return null;
	}

	public ModelAndView standardizeSubmit(HttpServletRequest request, HttpServletResponse response)
	{
		SetComparisonProcessor processor = getProcessor();
		DataPreprocessingParser.parseStandartizationUI(request, processor.standartization, null, null);
		processor.standartization.addExplicitHydrogensWith = processor.standartization.getDefault();
		processor.standartization.dearomatizeWith = processor.standartization.getDefault();
		return redirect("setcomparison/descriptors.do");
	}

	public ModelAndView descriptors(HttpServletRequest request, HttpServletResponse response)
	{
		DescriptorsConfigurator descConfigurator = new DescriptorsConfigurator();
		return descConfigurator.descriptorBlocks().setTemplate("setcomparison/sc-descriptors").getModelAndView();
	}

	public ModelAndView descriptorsSubmit(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		SetComparisonProcessor processor = getProcessor();
		DescriptorsConfigurator descConfigurator = new DescriptorsConfigurator();
		CDSConfiguration configuration = descConfigurator.configureDescriptors(request, new CDSConfiguration());
		processor.descriptors = configuration.descriptors;

		processor.start();
		processor.operation.successURL = "setcomparison/results.do";
		Globals.setSessionAttribute(SessionVariable.SETCOMPARE_PROCESSOR, processor);

		EventFactory.document("SetCompare", null, "is running the SetCompare utility using descriptors: " + configuration.descriptors.types);

		return redirect("longoperations/operationWaitingScreen.do?operation-id="
				+ processor.operation.operationId);
	}

	public ModelAndView results(HttpServletRequest request, HttpServletResponse response) throws Exception 
	{
		SetComparisonProcessor processor = getProcessor();
		if (assertParam("task"))
		{
			processor = new SetComparisonProcessor();
			processor.restoreFromPendingTask(PendingTask.getById(getLongParam("task")));
			Globals.setSessionAttribute(SessionVariable.SETCOMPARE_PROCESSOR, processor);
		}
		WebModel wm = new WebModel();
		wm.addParam("set1", (processor.basket1.entries.size() - processor.getErrorIndices(0).size()) + " unique molecules");
		wm.addParam("set2", (processor.basket2.entries.size() - processor.getErrorIndices(1).size()) + " unique molecules");

		wm.addParam("errors1", "" + processor.getErrorIndices(0).size());
		wm.addParam("errors2", "" + processor.getErrorIndices(1).size());

		return wm.setTemplate("setcomparison/sc-results").getModelAndView();
	}

	public ModelAndView molecules(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		return new WebModel().setTemplate("setcomparison/sc-molecules").getModelAndView();
	}

	private List<Integer> getMoleculeIndices()
	{
		SetComparisonProcessor processor = getProcessor();

		String descName = getParam("name");
		int set = getIntParam("set") - 1;
		List<Integer> indices;
		if (assertParam("errors"))
			indices = processor.getErrorIndices(set);
		else
			indices = processor.getMoleculeIndices(descName, set);

		return indices;
	}

	public ModelAndView listMolecules(HttpServletRequest request, HttpServletResponse response) throws Exception 
	{
		SetComparisonProcessor processor = getProcessor();

		int set = getIntParam("set") - 1;
		List<Integer> indices = getMoleculeIndices();
		Basket basket = (set == 0) ? processor.basket1 : processor.basket2;

		WebList wl = new WebList();
		wl.pageNum = getPageNum();
		wl.pageSize = getPageSize(100);
		wl.size = indices.size();
		wl.list = new ArrayList<Object>();

		for (int i = (wl.pageNum - 1) * wl.pageSize; i < Math.min(wl.size , wl.pageNum * wl.pageSize); i++)
		{
			int index = indices.get(i);

			Molecule mol = basket.entries.get(index).ep.molecule;
			if (assertParam("errors"))
				mol.error = processor.getErrorMessage(set, index);
			wl.list.add(mol);
		}

		WebFilters filters = formFilters(request);
		return new BrowserModel().setFilters(filters).setObject(wl).getModelAndView();
	}

	public SetComparisonController()
	{
		sessionRequired = true;
	}

	/**
	 * Get the currently active processor saved in the user session
	 */
	private SetComparisonProcessor getProcessor()
	{
		SetComparisonProcessor processor = (SetComparisonProcessor) Globals.getSessionAttribute(SessionVariable.SETCOMPARE_PROCESSOR);
		if (processor == null)
			Globals.setSessionAttribute(SessionVariable.SETCOMPARE_PROCESSOR, processor = new SetComparisonProcessor());

		return processor;
	}

	@Override
	public ModelAndView list(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		SetComparisonProcessor processor = getProcessor();
		DataTable dtResult = processor.wndResult.ports.get(0);
		List<Object> results = new ArrayList<Object>();

		int pageNum = getPageNum();
		int pagesize = getPageSize(15);
		for (int i = (pageNum - 1) * pagesize; i < Math.min(dtResult.getRowsSize(), pageNum * pagesize); i++)
			results.add(new SetCompareResult(dtResult.getRow(i)).setId(i));

		WebList list = new WebList();
		list.list = results;
		list.size = dtResult.getRowsSize();
		list.pageNum = pageNum;
		list.pageSize = pagesize;

		WebFilters filters = formFilters(req);
		return new BrowserModel().setFilters(filters).setObject(list).getModelAndView();
	}

	public ModelAndView openEPBrowser(HttpServletRequest request, HttpServletResponse response)
	{
		SetComparisonProcessor processor = getProcessor();
		List<Integer> molIndices = getMoleculeIndices();

		int set = getIntParam("set") - 1;

		Basket basket = (set == 0) ? processor.basket1 : processor.basket2;

		List<Integer> mp2IDs = new ArrayList<Integer>();
		for (Integer index : molIndices)
			mp2IDs.add(basket.entries.get(index).ep.molecule.mapping2.id);

		int filterID = Mapping2Filter.generateFilterID();
		Mapping2Filter.addCompoundToFilter(filterID, mp2IDs);

		return new WebModel(new Alert("epbrowser/show.do?mol-filter=" + filterID)).getModelAndView();
	}
}
