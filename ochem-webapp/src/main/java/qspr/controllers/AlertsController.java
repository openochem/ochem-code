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

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.OutputStream;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.HibernateException;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.SessionVariable;
import qspr.ThreadScope;
import qspr.business.WebFilters;
import qspr.business.toxalert.AlertsFilter;
import qspr.business.toxalert.ScreeningProcessor;
import qspr.entities.Alert;
import qspr.entities.AlertSubstitutionVariable;
import qspr.entities.Article;
import qspr.entities.Basket;
import qspr.entities.BasketEntry;
import qspr.entities.Mapping2Filter;
import qspr.entities.PendingTask;
import qspr.entities.Property;
import qspr.entities.SubstructureAlert;
import qspr.entities.Tag;
import qspr.export.ExportableColumn;
import qspr.export.ExportableSet;
import qspr.export.ExportableSetConfiguration;
import qspr.frontend.AvailableAlertsData;
import qspr.frontend.BrowserModel;
import qspr.frontend.LabeledValue;
import qspr.frontend.TreeNode;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;
import qspr.metaserver.Event;
import qspr.metaserver.EventListener;
import qspr.metaserver.util.ExtendedSMART;
import qspr.modelling.CompoundsProvider;
import qspr.modelling.DataPreprocessingParser;
import qspr.util.AccessChecker;
import qspr.util.ExportThread;
import qspr.util.RWriter;
import qspr.workflow.datatypes.DataRow;
import qspr.OCHEMConfiguration;

import com.eadmet.business.AlertsService;
import com.eadmet.business.AlertsUploadReport;
import com.eadmet.business.SmartsMatcher;
import com.eadmet.business.TagsService;
import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.useractions.EventFactory;
import com.eadmet.utils.MAILERConstants;
import com.eadmet.utils.mailer.Email;
import com.eadmet.utils.mailer.Mailer;


/**
 * A controller servicing the ToxAlerts front-end
 * @author midnighter
 *
 */
@Controller
@SuppressWarnings({"unchecked","rawtypes"})
public class AlertsController extends BrowserWrapper 
{
	private static transient final Logger logger = LogManager.getLogger(AlertsController.class);

	public AlertsController()
	{
		sessionRequired = true;
	}

	public ModelAndView show(HttpServletRequest req, HttpServletResponse res) 
	{
		WebModel wm = new WebModel();
		wm.addObject(new AvailableAlertsData());
		if (isModerator())
			wm.addParam("moderator", "true");
		return wm.setTemplate("browsers/sa-browser").getModelAndView();
	}

	public ModelAndView autocomplete(HttpServletRequest request, HttpServletResponse response) 
			throws Exception
	{
		Property functionalGroup = Property.getByName("Functional groups");
		List<SubstructureAlert> alerts = Globals.session().createCriteria(SubstructureAlert.class)
				.add(Restrictions.like("name", getParam("term") + "%"))
				.add(Restrictions.eq("property", functionalGroup))
				.setMaxResults(10).list();

		List<LabeledValue> labels = new ArrayList<LabeledValue>();

		for (SubstructureAlert alert : alerts)
			labels.add(new LabeledValue(alert.name, alert.id));

		return new WebModel().addObjects(labels).getModelAndView();
	}

	@Override
	public ModelAndView list(HttpServletRequest req, HttpServletResponse res)
			throws Exception {
		Criteria criteria = Globals.session().createCriteria(
				SubstructureAlert.class);

		getAlertsFilter().filterCriteria(criteria);

		WebList list = new WebList();
		list.useEntity(SubstructureAlert.class);
		list.loadDistinctFromCriteria(criteria, getPageNum(), getPageSize(100));
		WebFilters filters = this.formFilters(req);
		return new BrowserModel().setSelectionSize(Globals.userSession().selectedAlerts.size()).setFilters(filters).setObject(list)
				.getModelAndView();
	}

	/**
	 * Export all the alerts into a CSV file
	 * @param req
	 * @param res
	 * @return
	 * @throws Exception
	 */
	public ModelAndView export(HttpServletRequest req, HttpServletResponse res)
			throws Exception 
	{
		if (!isModerator())
			throw new UserFriendlyException("Only ToxAlerts moderators are allowed to access this functionality");

		OutputStream os = res.getOutputStream();
		res.setContentType("application/csv");
		res.setHeader("Content-Disposition", "attachment; filename=alerts.csv");

		AlertsService.exportAlerts(os);

		return null;
	}

	public ModelAndView exportResults(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		ExportableSet eData = new ExportableSet();
		eData.clearColumns();
		eData.addColumn(ExportableColumn.SMILES);
		eData.addColumn(ExportableColumn.MOLECULEID);
		eData.addColumn(ExportableColumn.DESCRIPTORS);

		if (assertParam("submit"))
		{
			final ScreeningProcessor processor = getProcessor();
			ExportThread eThread = AlertsService.getResultExporter(getParam("format"), ExportableSetConfiguration.configureFromDialog(req), processor);
			eThread.start();

			return redirect("longoperations/operationWaitingScreen.do?operation-id=" + eThread.operationID);
		}
		else
			return new WebModel(eData).setTemplate("export").getModelAndView();
	}

	private AlertsFilter getAlertsFilter() throws Exception
	{
		AlertsFilter filter = new AlertsFilter();
		filter.parseUI(ThreadScope.get().getHttpServletRequest());
		return filter;
	}

	public ModelAndView action(HttpServletRequest req, HttpServletResponse res) throws HibernateException, Exception 
	{
		String action = getParam("action");
		AlertsFilter filter = getAlertsFilter();
		if (!action.contains("select") && filter.alertId == null)
			filter.selectedOnly = true;

		if (action.contains("approve") && !isModerator())
			throw new UserFriendlyException("Moderator privileges are required for this action");

		if (action.contains("select"))
		{
			List<Long> alertIds = filter.filterCriteria().setProjection(Projections.property("id")).list();
			if ("addselect".equals(action))
				Globals.userSession().selectedAlerts.addAll(alertIds);
			else if ("removeselect".equals(action))
				Globals.userSession().selectedAlerts.removeAll(alertIds);
			else if ("toggleselect".equals(action))
				for (Long alertId : alertIds) 
				{
					if (Globals.userSession().selectedAlerts.contains(alertId))
						Globals.userSession().selectedAlerts.remove(alertId);
					else
						Globals.userSession().selectedAlerts.add(alertId);
				}
		}
		else 
		{
			List<SubstructureAlert> alerts = filter.filterCriteria().list();
			for (SubstructureAlert alert : alerts) 
			{
				if (!action.contains("select"))
					AccessChecker.requestModificationPermission(alert);

				if ("delete".equals(action))
					Globals.session().delete(alert);
				else if ("approve".equals(action) && isModerator())
					alert.approved = true;
				else if ("disapprove".equals(action))
					alert.approved = false;
			}
		}
		return new WebModel().setTemplate("toxalerts/sa-upload").getModelAndView();
	}

	public ModelAndView contents(HttpServletRequest req, HttpServletResponse res) throws HibernateException, Exception 
	{
		Globals.session().createSQLQuery("update SubstructureAlert set number=sa_id where number is null").executeUpdate();

		AlertsFilter filter = new AlertsFilter();
		filter.includeFolders = true;
		filter.endpointId = Property.getByName("Functional groups").id;

		List<TreeNode> nodes = new ArrayList<TreeNode>(); 
		Criteria c = filter.filterCriteria();
		c.addOrder(Order.asc("number"));
		List<SubstructureAlert> alerts = c.list();

		Map<SubstructureAlert, TreeNode> nodesMap = new HashMap<SubstructureAlert, TreeNode>();

		for (SubstructureAlert alert : alerts)
		{
			TreeNode node = new TreeNode();
			node.object = alert;
			nodes.add(node);
			nodesMap.put(alert, node);
		}

		for (TreeNode node : nodes)
		{
			for (SubstructureAlert alert : ((SubstructureAlert)node.object).parents)
			{
				node.parents.add(nodesMap.get(alert));
				nodesMap.get(alert).children.add(node);
			}
			node.object = ((SubstructureAlert) node.object).getSimpleAlert();
		}

		TreeNode root = new TreeNode();

		SubstructureAlert rootAlert = new SubstructureAlert();
		rootAlert.name = "Contents";
		rootAlert.id = -1L;
		rootAlert.folder = true;
		root.object = rootAlert;
		for (TreeNode node : nodes)
		{
			SubstructureAlert alert = (SubstructureAlert) node.object;
			if (node.parents.isEmpty() && (!node.children.isEmpty() || alert.folder))
				root.children.add(node);
		}



		WebModel wm = new WebModel(root);
		return wm.setTemplate("toxalerts/sa-contents").getModelAndView();
	} 

	public ModelAndView screen(HttpServletRequest req, HttpServletResponse res) 
	{
		WebModel wm = new WebModel();
		wm.addObject(new AvailableAlertsData());
		return wm.setTemplate("toxalerts/sa-screen").getModelAndView();
	}

	public ModelAndView screenSubmit(HttpServletRequest req, HttpServletResponse res) throws Exception {

		CompoundsProvider provider = new CompoundsProvider();
		provider.parseUI();
		if (provider.error != null)
			return new WebModel(new Alert(provider.error)).setTemplate("toxalerts/sa-screen").getModelAndView();

		ScreeningProcessor processor = new ScreeningProcessor();

		processor.setDescription = provider.setDescription;
		processor.alertsFilter = getAlertsFilter();
		processor.compoundsProvider = provider;

		DataPreprocessingParser.parseStandartizationUI(req, processor.standartization, null, null);
		processor.standartization.addExplicitHydrogensWith = processor.standartization.getDefault();
		processor.standartization.dearomatizeWith = processor.standartization.getDefault();

		processor.start();
		processor.operation.successURL = "alerts/screenResults.do";
		Globals.setSessionAttribute(SessionVariable.ALERT_SCREENING_PROCESSOR, processor);

		provider.basketLoaded.addListener(new EventListener<Basket>()
		{
			@Override
			public void onEvent(Event event, Basket basket)
			{
				EventFactory.document("ToxAlerts screening", null, String.format("has started alerts screening for %d compounds", basket.getRowsSize()));
			}
		});

		return redirect("longoperations/operationWaitingScreen.do?operation-id="
				+ processor.operation.operationId);
	}

	public ModelAndView image(HttpServletRequest req, HttpServletResponse res)
			throws FileNotFoundException, IOException {
		SubstructureAlert alert = SubstructureAlert.getByID(getLongParam("id"));
		if (alert.image != null) {
			res.setHeader("Content-type", "image/png");
			res.getOutputStream().write(alert.image);
		} else
			res.sendRedirect("/img/no-image.jpg");

		return null;
	}

	public ModelAndView uploadImage(HttpServletRequest req, HttpServletResponse res) throws Exception,
	FileNotFoundException, IOException {
		SubstructureAlert alert = SubstructureAlert.getByID(getLongParam("id"));
		File f = Globals.getUploadedFile();
		if (f != null) {
			logger.info("Updating image for alert " + alert.id);
			byte[] image = new byte[(int) f.length()];
			FileInputStream ff = new FileInputStream(f);
			ff.read(image);
			alert.image = image;
			Globals.session().saveOrUpdate(alert);
			ff.close();
		}
		return null;
	}

	public ModelAndView match(HttpServletRequest req,
			HttpServletResponse res) throws Exception 
	{
		if (assertParam("smiles"))
		{
			SmartsMatcher matcher = SmartsMatcher.get(OCHEMConfiguration.getCheminfEngine());
			return new WebModel(new Alert(matcher.match(getParam("smiles"), getParam("smarts")) ? "1" : "0")).getModelAndView();
		}

		return new WebModel().setTemplate("smart-matcher").getModelAndView();
	}


	public ModelAndView screenResults(HttpServletRequest req,
			HttpServletResponse res) throws Exception {
		ScreeningProcessor processor = getProcessor();

		if (assertParam("task"))
		{
			processor = new ScreeningProcessor();
			processor.restoreFromPendingTask(PendingTask.getById(getLongParam("task")));
			Globals.setSessionAttribute(SessionVariable.ALERT_SCREENING_PROCESSOR, processor);
		}

		// Count compounds by alerts
		List<Object> objects = new ArrayList<Object>();
		for (int i = 0; i < processor.compoundsByAlerts.index.size(); i++) {
			SubstructureAlert alert = SubstructureAlert.getByID(processor.compoundsByAlerts.index.get(i));
			alert.countOccurences = processor.compoundsByAlerts.values.get(i).size();
			objects.add(alert);
		}

		for (int i = 0; i < processor.compoundsByEndpoints.index.size(); i++) {
			Property endpoint = Property.getById(processor.compoundsByEndpoints.index.get(i));
			endpoint.count = 1L * processor.compoundsByEndpoints.values.get(i).size();
			objects.add(endpoint);
		}

		for (int i = 0; i < processor.compoundsByPublications.index.size(); i++) {
			Article article = Article.getById(processor.compoundsByPublications.index.get(i));
			article.count = 1L * processor.compoundsByPublications.values.get(i).size();
			objects.add(article);
		}

		return new WebModel().addObjects(objects)
				.setTemplate("toxalerts/sa-screen-results").getModelAndView();
	}

	// Get indexes of the screened compounds filtered by an alert, by a publication, by an endpoint or by nothing
	private List<Integer> getScreenResultCompoundIndex(final ScreeningProcessor processor)
	{
		if (assertParam("alert")) 
		{
			if (getParam("alert").startsWith("e"))
				return processor.compoundsByEndpoints.values
						.get(processor.compoundsByEndpoints.index
								.indexOf(Long.valueOf(getParam("alert")
										.substring(1))));
			else if (getParam("alert").startsWith("a"))
				return processor.compoundsByPublications.values
						.get(processor.compoundsByPublications.index
								.indexOf(Long.valueOf(getParam("alert")
										.substring(1))));
			else
				return processor.compoundsByAlerts.values
						.get(processor.compoundsByAlerts.index
								.indexOf(getLongParam("alert")));
		}
		else
			return new AbstractList<Integer>() {
			@Override
			public Integer get(int index) {
				return index;
			}

			@Override
			public int size() {
				return processor.alertsByCompounds.index.size();
			}
		};
	}

	public ModelAndView screenlist(HttpServletRequest req,
			HttpServletResponse res) 
	{
		ScreeningProcessor processor = getProcessor();
		int pagenum = getPageNum();
		int pagesize = getPageSize(15);

		ArrayList<Object> result = new ArrayList<Object>();
		WebList list = new WebList();
		list.list = result;
		list.pageNum = pagenum;
		list.pageSize = pagesize;

		List<Integer> compoundIndices = getScreenResultCompoundIndex(processor);
		for (int i = (pagenum - 1) * pagesize; i < Math.min(
				compoundIndices.size(), pagenum * pagesize); i++) {
			int cmpIndex = compoundIndices.get(i);
			DataRow row = new DataRow();
			row.addValue(processor.alertsByCompounds.index.get(cmpIndex));

			for (Long alertID : processor.alertsByCompounds.values
					.get(cmpIndex))
				row.addValue(SubstructureAlert.getByID(alertID));
			result.add(row);
		}
		list.size = compoundIndices.size();

		return new WebModel(list).getModelAndView();
	}

	// An insider function to export screening results into R
	public ModelAndView exportResultsInR(HttpServletRequest request, HttpServletResponse response) throws IOException 
	{
		ScreeningProcessor processor = getProcessor();
		response.setContentType("application/r");
		response.setHeader("Content-Disposition", "attachment; filename=\"toxalerts-results-"+processor.setDescription+".R\"");
		RWriter writer = new RWriter(response.getOutputStream());

		for (BasketEntry be : processor.compoundsProvider.getBasket().entries) {
			writer.addValue("compoundIDs", be.ep.molecule.mapping2.id);
		}

		// Compounds with any alerts
		for (int i = 0; i < processor.alertsByCompounds.index.size(); i++)
			writer.addValue("compoundsWithAlerts", processor.alertsByCompounds.index.get(i));

		// Compounds by alerts
		for (int i = 0; i < processor.compoundsByAlerts.index.size(); i++)
			for (int k = 0; k < processor.compoundsByAlerts.values.get(i).size(); k++)
			{
				writer.addValue("alerts$id", processor.compoundsByAlerts.index.get(i));
				writer.addValue("alerts$compound", processor.compoundsByAlerts.values.get(i).get(k));
			}

		// Compounds by endpoints
		for (int i = 0; i < processor.compoundsByEndpoints.index.size(); i++)
			for (int k = 0; k < processor.compoundsByEndpoints.values.get(i).size(); k++)
			{
				writer.addValue("endpoints$id", processor.compoundsByEndpoints.index.get(i));
				writer.addValue("endpoints$compound", processor.compoundsByEndpoints.values.get(i).get(k));
			}

		// Compounds by publications
		for (int i = 0; i < processor.compoundsByPublications.index.size(); i++)
			for (int k = 0; k < processor.compoundsByPublications.values.get(i).size(); k++)
			{
				writer.addValue("publications$id", processor.compoundsByPublications.index.get(i));
				writer.addValue("publications$compound", processor.compoundsByPublications.values.get(i).get(k));
			}

		writer.writeAndClose();

		return null;
	}

	public ModelAndView edit(HttpServletRequest req, HttpServletResponse res) {
		SubstructureAlert alert = SubstructureAlert.getByID(getLongParam("id"));
		return new WebModel(alert).setTemplate("toxalerts/sa-edit").getModelAndView();
	}

	public ModelAndView contentsAction(HttpServletRequest req, HttpServletResponse res) throws Exception {
		SubstructureAlert folder = SubstructureAlert.getByID(getLongParam("id"));
		String action = getParam("action");
		if ("delete".equals(action))
			Globals.session().delete(folder);
		else if ("addchild".equals(action))
		{	
			SubstructureAlert child = new SubstructureAlert();
			child.name = getParam("name");
			child.folder = true;
			if (folder != null)
				child.parents.add(folder);
			child.approved = true;
			child.introducer = child.owner = Globals.userSession().user;
			child.article = Article.getByPmId(12345);
			child.property = Property.getByName("Functional groups");
			child.updateHash();
			Globals.session().save(child);	
			return new WebModel(child).getModelAndView();
		}
		else if ("rename".equals(action)){
			folder.name = getParam("name");
			folder.updateHash();
			Globals.session().saveOrUpdate(folder);	
		}
		else if ("move".equals(action))
		{
			folder.parents.clear(); // FIXME: Consider multiple parents case
			if (getLongParam("parent") > 0)
				folder.parents.add(SubstructureAlert.getByID(getLongParam("parent")));
			if (assertParam("replaced-node"))
			{
				SubstructureAlert replacedNode = SubstructureAlert.getByID(getLongParam("replaced-node"));
				folder.number = replacedNode.number;
				Globals.session().createSQLQuery("update SubstructureAlert set number=number+1 where number >= " + replacedNode.number).executeUpdate();
			}
			Globals.session().saveOrUpdate(folder);	 
		}
		return new WebModel().getModelAndView();
	}

	public ModelAndView addFolder(HttpServletRequest req, HttpServletResponse res) throws Exception {

		return new WebModel().getModelAndView();
	}

	public ModelAndView deleteFolder(HttpServletRequest req, HttpServletResponse res) throws Exception {
		SubstructureAlert folder = SubstructureAlert.getByID(getLongParam("id"));
		Globals.session().delete(folder);

		return new WebModel().getModelAndView();
	}

	public ModelAndView save(HttpServletRequest req, HttpServletResponse res) {
		SubstructureAlert alert = SubstructureAlert.getByID(getLongParam("id"));
		AccessChecker.requestModificationPermission(alert);
		if (assertParam("smarts"))
			alert.smart = getParam("smarts");
		else
			throw new UserFriendlyException("SMARTS cannot be empty");
		alert.name = getParam("name");
		alert.description = getParam("description");
		alert.smartsDescription = getParam("smarts-description");
		alert.comment = getParam("comment");
		alert.property = Property.getById(getLongParam("property"));
		alert.article = Article.getById(getLongParam("article"));

		alert.parents.clear();
		if (assertParam("parents"))
		{
			String[] parentIDs = getParam("parents").split(",");
			for (String parentID : parentIDs)
				alert.parents.add(SubstructureAlert.getByID(Long.valueOf(parentID)));
		}


		if (ExtendedSMART.create(alert.getFullSMARTS(), OCHEMConfiguration.getCheminfEngine()).invalid)
			throw new UserFriendlyException("SMART " + alert.smart
					+ " is not a valid SMARTS string");

		alert.updateHash();
		Globals.session().saveOrUpdate(alert);

		return new WebModel(alert).setTemplate("toxalerts/sa-edit").getModelAndView();
	}

	public ModelAndView home(HttpServletRequest req, HttpServletResponse res) {
		return new WebModel().setTemplate("toxalerts/sa-home").getModelAndView();
	} 

	public ModelAndView upload(HttpServletRequest req, HttpServletResponse res) {
		return new WebModel().setTemplate("toxalerts/sa-upload").getModelAndView();
	}

	public ModelAndView variables(HttpServletRequest req, HttpServletResponse res) 
	{
		List<AlertSubstitutionVariable> variables = Globals.session().createCriteria(AlertSubstitutionVariable.class).list();
		return new WebModel()
				.setList(variables)
				.setTemplate("toxalerts/sa-variables")
				.getModelAndView();
	}

	// S.A. upload dublicates many features of the normal batch upload. This
	// part can be refactored and some code can be re-used
	public ModelAndView uploadSubmit(HttpServletRequest req, HttpServletResponse res) throws NumberFormatException, Exception {

		AccessChecker.requestRegisteredUserPrivileges();

		File f = null;

		try
		{
			f = Globals.getUploadedFile();
		} catch (Exception e)
		{
			return redirect("alerts/upload.do");
		}

		AlertsUploadReport report = AlertsService.uploadAlerts(new FileInputStream(f), f.getName(), assertParam("private-upload"), null, null);

		String userLogin = (null == Globals.userSession().user) ? "Guest" : Globals.userSession().user.login;
		notifyModerators("New alerts have been uploaded into ToxAlerts!", "User " + userLogin + " has uploaded new alerts and they are awaiting your approval.\nDetails: " + report.toString() + ".\n\nThank you for your kind attention.");
		EventFactory.document("ToxAlerts upload", null, String.format("has uploaded %d structural alerts into the ToxAlerts database", report.successes));

		return new WebModel(new Alert(report.toString())).setTemplate("toxalerts/sa-upload-results")
				.getModelAndView();
	}

	public ModelAndView openPropertyBrowser(HttpServletRequest req, HttpServletResponse res) 
	{
		ScreeningProcessor processor = getProcessor();
		int filterID = Mapping2Filter.generateFilterID();
		List<Integer> compoundIndices = getScreenResultCompoundIndex(processor);
		for (int i = 0; i < compoundIndices.size(); i++)
			Mapping2Filter.addCompoundToFilter(filterID,
					processor.alertsByCompounds.index.get(compoundIndices
							.get(i)));

		return redirect("epbrowser/show.do?mol-filter=" + filterID);
	}

	// Tag the filtered molecules from screening results
	public ModelAndView addTag(HttpServletRequest req, HttpServletResponse res)
	{
		Tag tag = Tag.getByID(getLongParam("tag"));
		ScreeningProcessor processor = getProcessor();
		List<Integer> cmpIndex = getScreenResultCompoundIndex(processor);

		List<Integer> mp2IDs = new ArrayList<Integer>();
		for (int i = 0; i < cmpIndex.size(); i++)
			mp2IDs.add(processor.alertsByCompounds.index.get(cmpIndex.get(i)));

		TagsService.addTagByMP2(tag, mp2IDs);

		return new WebModel().getModelAndView();
	}

	private boolean isModerator()
	{
		if (Globals.isSuperUser())
			return true;

		List<String> toxalertsModerators = new ArrayList<String>();
		toxalertsModerators.add(MAILERConstants.ADMIN);
		return Globals.userSession().user != null && toxalertsModerators.contains(Globals.userSession().user.login);
	}

	private void notifyModerators(String subject, String message)
	{
		Mailer.postMailSafely(new Email(MAILERConstants.EMAIL_ADMIN, subject, "Dear ToxAlerts moderator,\n\nwe got news.\n" + message + "\n\nSincerely yours,\nToxAlerts server"));
	}

	protected ScreeningProcessor getProcessor()
	{
		return (ScreeningProcessor) Globals.getSessionAttribute(SessionVariable.ALERT_SCREENING_PROCESSOR);
	}

}




