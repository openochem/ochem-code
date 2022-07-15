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

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.SessionVariable;
import qspr.business.PendingTaskFilter;
import qspr.business.PendingTaskPeer;
import qspr.business.Privileges;
import qspr.business.SetComparisonProcessor;
import qspr.entities.Alert;
import qspr.entities.Article;
import qspr.entities.Attachment;
import qspr.entities.Attachment.AttachmentType;
import qspr.entities.AttachmentSource;
import qspr.entities.Model;
import qspr.entities.PendingTask;
import qspr.entities.PendingTask.TaskType;
import qspr.entities.ReadyModelAttachment;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;
import qspr.metaserver.CalculationClient;
import qspr.metaserver.protocol.Task;
import qspr.metaserver.protocol.Task.DebugLevel;
import qspr.modelling.ModelFactory;
import qspr.modelling.ModelProcessor;
import qspr.modelling.TaskPriorityManager;
import qspr.modelling.applier.ModelApplier;
import qspr.util.AccessChecker;
import qspr.util.WrapperThread;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.business.DescriptorsCalculatorProcessor;
import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.mailer.Mailer;

@Controller
public class PendingTasksController extends BrowserWrapper
{
	public PendingTasksController()
	{
		sessionRequired = true;
	}

	public ModelAndView fetch(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		return redirect("modelconfigurator/configure.do?page=save&model=" + getLongParam("id") + "&render-mode=popup");
	}

	public ModelAndView profile(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		PendingTask pTask = PendingTask.getById(getLongParam("id"));
		if (!pTask.getPrivileges(request).canView)
			throw new UserFriendlyException("You are not authorized to view this task");

		return new WebModel(pTask).setTemplate("pending-task-profile").getModelAndView();
	}

	public ModelAndView editProfile(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		PendingTask pTask = PendingTask.getById(getLongParam("id"));
		if (!pTask.getPrivileges(request).canEdit)
			throw new UserFriendlyException("You are not authorized to edit this task published by " + pTask.session.user.login);
		pTask.name = getParam("name");

		pTask.description = getParam("description");
		pTask.setDescription = getParam("set-description");

		if (assertParam("article"))
			pTask.article = Article.getById(getLongParam("article"));

		Globals.session().saveOrUpdate(pTask);

		return new WebModel().getModelAndView();
	}

	// This function will replace the old "fetch" after debugging
	public ModelAndView fetchnew(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		String suffix = "?task=" + getLongParam("id") + "&render-mode=popup";
		PendingTask pTask = (PendingTask) Globals.session().get(PendingTask.class, getLongParam("id"));
		if (!pTask.getPrivileges(request).canView)
			throw new UserFriendlyException("You cannot access this task");

		if (pTask.type == TaskType.MODEL_TRAINING)
			return redirect("modelconfigurator/configure.do?page=save&model=" + pTask.model.id + "&render-mode=popup");
		else if (pTask.type == TaskType.MODEL_APPLICATION)
			return redirect("modelapplier/results.do" + suffix);
		else if (pTask.type == TaskType.DESCRIPTOR_CALCULATION)
			return redirect("descriptorscalculator/results.do" + suffix);
		else if (pTask.type == TaskType.TOXALERT_SCREENING)
			return redirect("alerts/screenResults.do" + suffix);
		else if (pTask.type == TaskType.SET_COMPARISON)
			return redirect("setcomparison/results.do" + suffix);
		else
			throw new UserFriendlyException("Unsupported task type: " + pTask.type);
	}

	public ModelAndView recalculateAllFailed(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		BatchRecalculator recalculator = new BatchRecalculator();
		recalculator.filter = getTaskFilter();
		recalculator.start();

		return new WebModel().getModelAndView();
	}

	public ModelAndView deleteAll(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		BatchDeleter batcher = new BatchDeleter();
		batcher.filter = getTaskFilter();
		batcher.start();

		return new WebModel().getModelAndView();
	}

	public ModelAndView fetchReady(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		BatchFetcher fetcher = new BatchFetcher();
		fetcher.filter = getTaskFilter();
		fetcher.start();

		return new WebModel().getModelAndView();
	}

	public ModelAndView action(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		String action = getParam("action");
		PendingTask pTask = (PendingTask) Globals.session().get(PendingTask.class, getLongParam("id"));
		if (pTask == null)
			throw new UserFriendlyException("This item has already been deleted");

		if (!pTask.getPrivileges(request).canEdit)
			throw new UserFriendlyException("Not sufficient privileges for this action" + (pTask.getPrivileges(request).canView?" model was created by "+
					pTask.session.user.login:""));

		if ("delete".equals(action))
		{
			PendingTaskPeer.deletePendingTask(pTask,false);
		}
		else if ("recalculate".equals(action))
		{
			recalculateTask(pTask, assertParam("debug"));
			Globals.restartAllTransactions(true);
		}
		else if (Task.KILL.equals(action))
		{
			if (pTask.taskId == 0)
				WrapperThread.getPendingTaskThread(pTask.id).interrupt();
			else
				getClient().killTask(pTask.taskId);
		}
		else if ("set_priority".equals(action))
		{
			int newPriority = getIntParam("priority");
			newPriority = TaskPriorityManager.getNewTaskPriority(newPriority, pTask, Globals.userSession());
			if (Task.isAliveStatus(pTask.status))
			{
				if (pTask.taskId != null)
					getClient().setTaskPriority(pTask.taskId, newPriority);
				else
					ModelProcessor.getPendingTaskThread(pTask.id).pTask.setPriority(newPriority);
			}
			pTask.setPriority(newPriority);
			Globals.session().saveOrUpdate(pTask);
		}
		else if ("publish".equals(action))
		{
			pTask.published = true;
			if (pTask.taskId > 0) {
				if(pTask.article == null) throw new UserFriendlyException("Before publication click edit button and select an article in which this task will be published." +
						" \n\nOnce it is published only the OCHEM administrators can unpublish it (after some time). Do not publish tasks unless they are important to illustrate your study.");
				Task task = pTask.retrieveTask(false);
				task.setLocalTask();
				task.setResult(task.getResult()); // will be storing results in the local database
				pTask.readyTask = new Attachment<Task>(task, AttachmentType.SERIALIZABLE, AttachmentSource.ReadyPendingTask);
				pTask.attachment.publish();
				Mailer.notifyAdmins("A new task is published", "A task " + pTask.name +  " has been published by user " + pTask.session.user.login);
			}
			Globals.session().saveOrUpdate(pTask);
		}

		else if ("unpublish".equals(action))
		{
			try{
				if(pTask.article != null) PendingTask.retrieveTask(pTask.taskId, false); // if we can still retrieve task; we allow to unpublish it
				pTask.published = false;
				Mailer.notifyAdmins("A task is unpublished", "A task " + pTask.name +  " has been un-published by user " + pTask.session.user.login);
				Globals.session().saveOrUpdate(pTask);
			}catch(Exception e){
				throw new UserFriendlyException("This task can be only deleted by the administrator.");
			};
		}
		else
			throw new UserFriendlyException("Unknown action requested - " + action);

		return new WebModel().getModelAndView();
	}

	private void recalculateTask(PendingTask pTask, boolean debug) throws IOException, ClassNotFoundException, Exception
	{
		if (WrapperThread.getPendingTaskThread(pTask.id) != null)
			throw new UserFriendlyException("Recalculation of this task has already been initiated");

		pTask.readyTask = null; // it is not yet ready task anyway...

		if (pTask.type == TaskType.MODEL_TRAINING)
		{
			logger.info("Recalculating model " + pTask.model.publicId);
			ModelProcessor proc = ModelFactory.getProcessor(pTask.model.template);
			proc.taskDebugLevel = debug ? DebugLevel.ERRORS : DebugLevel.NONE;
			proc.defaultTaskPriority = pTask.getPriority();
			proc.model = pTask.model;
			proc.pTask = pTask;
			pTask.taskId = 0;
			pTask.status = "init";
			pTask.detailedStatus = ("Initializing recalculation...");
			Globals.restartAllTransactions(true);
			proc.updateSessionTime = false;
			proc.start();
			proc.exitAfterPosting();
		}
		else if (pTask.type == TaskType.MODEL_APPLICATION)
		{
			pTask.taskId = 0;
			pTask.status = "init";
			pTask.detailedStatus = "Initializing recalculation...";
			Globals.restartAllTransactions(true);
			ModelApplier applier = new ModelApplier(pTask);
			applier.taskDebugLevel = debug ? DebugLevel.ERRORS : DebugLevel.NONE;
			applier.defaultTaskPriority = pTask.getPriority();
			applier.start();
		}
		else if (pTask.type == TaskType.TOXALERT_SCREENING)
		{
			throw new UserFriendlyException("Resubmission of this kind of tasks is not yet supported!"); // Need to think how to load molecules properly
			//ScreeningProcessor proc = new ScreeningProcessor(pTask);
			//proc.start();
		}
		else if (pTask.type == TaskType.SET_COMPARISON)
		{
			pTask.status = "init";
			pTask.taskId = 0;
			Globals.restartAllTransactions(true);
			SetComparisonProcessor processor = new SetComparisonProcessor();
			processor.pTask = pTask;
			processor.start();
		}	
		else if (pTask.type == TaskType.DESCRIPTOR_CALCULATION)
		{
			pTask.status = "init";
			pTask.taskId = 0;
			Globals.restartAllTransactions(true);
			DescriptorsCalculatorProcessor processor = new DescriptorsCalculatorProcessor();
			processor.pTask = pTask;
			processor.taskDebugLevel = debug ? DebugLevel.ALL : DebugLevel.NONE;
			processor.start();
		}
	}

	public ModelAndView kill(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		Model model = (Model) Globals.session().get(Model.class, getLongParam("id"));

		if (model.taskId == null)
			throw new UserFriendlyException("There is no task to kill!");

		getClient().killTask(model.taskId);

		return new WebModel().setTemplate("model/recalculate").getModelAndView();
	}

	public ModelAndView recalculate(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		ModelApplier mData = new ModelApplier();
		Globals.setSessionAttribute(SessionVariable.MODEL_APPLIER, mData);
		Model model = (Model) Globals.session().get(Model.class, getLongParam("id"));
		if(model.published)throw new UserFriendlyException("This model has been alraedy published and can't be recalculated!");

		Privileges privileges = model.getPrivileges(request);

		if (!privileges.canEdit)
			throw new UserFriendlyException("You dont have sufficient privileges for this action");

		if (!PendingTask.getByModel(model, TaskType.MODEL_TRAINING).isEmpty())
			throw new UserFriendlyException("Model is already being calculated");

		mData.teacher = ModelFactory.getProcessor(model.template);
		mData.teacher.model = model;
		//model.taskId = 0;
		//model.status = Task.INIT;
		//Globals.session().saveOrUpdate(model);
		Globals.restartAllTransactions(true);
		//mData.teacher.threadId = model.id;
		mData.teacher.start();

		return new WebModel().setTemplate("model/recalculate").getModelAndView();
	}

	public ModelAndView recalculatestatus(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		ModelApplier mbData = (ModelApplier) Globals.getSessionAttribute(SessionVariable.MODEL_APPLIER);
		String status = "";
		if (mbData == null)
			status = "Modeller is not initialized";
		else
			status = mbData.teacher.getStatus();

		if ("Finished".equals(mbData.teacher.getStatus()))
			Globals.session().saveOrUpdate(mbData.teacher.model);

		return new WebModel(new Alert(status)).getModelAndView();
	}

	private PendingTaskFilter getTaskFilter()
	{
		PendingTaskFilter filter = new PendingTaskFilter();
		if (assertParam("set"))
			filter.setID = getLongParam("set");

		if (assertParam("model"))
			filter.modelID = getLongParam("model");

		if (assertParam("toBeDeleted"))
			filter.toBeDeleted = true;

		filter.status = getParam("status");
		if (assertParam("task-type"))
			filter.taskType = TaskType.valueOf(getParam("task-type"));

		if (assertParam("failures"))
			filter.status = "error,kill,killed";

		filter.taskName = getParam("task-name");
		filter.publishedOnly = assertParam("published");
		filter.ownTasksOnly = !assertParam("published") || assertParam("own-tasks-only");
		filter.showGroupModels = assertParam("group");
		filter.unpublishedTasksOfOtherUsers = assertParam("other-users") || assertParam("failures");
		if (assertParam("article-id"))
			filter.articleID = Long.valueOf(getParam("article-id").replaceAll("A", "").trim());

		if (assertParam("published"))
			filter.status = "ready";

		return filter;
	}

	@Override
	public ModelAndView list(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		PendingTaskFilter filter = getTaskFilter();
		if (filter.status != null)
			PendingTaskPeer.updateTaskStatuses(Globals.userSession().user, null, 200);

		WebList wl = new WebList();
		wl.useEntity(PendingTask.class);
		wl.loadDistinctFromCriteria(filter.createCriteria(), getPageNum(), getPageSize(15));

		// Distinguish running threads and tasks
		ArrayList<Integer> taskIdsToRequest = new ArrayList<Integer>();
		for (Object obj : wl.list)
		{
			PendingTask pt = (PendingTask) obj;

			if (pt.taskId > 0)
			{
				// Reguest task status only if its not published
				if (!pt.published || pt.readyTask == null)
					taskIdsToRequest.add(pt.taskId);
			}
			else if (pt.taskId == -1)
			{
				// A dummy ready task
				pt.status = "ready";
			}
			else
			{
				// A running task not yet posted to MetaServer
				WrapperThread thread = WrapperThread.getPendingTaskThread(pt.id);
				if (thread == null)
					pt.status = "error";
				else
					pt.detailedStatus = thread.getStatus();
			}
		}

		// Update the statuses
		if (!assertParam("status") && !taskIdsToRequest.isEmpty())
		{
			List<Task> tasks = getClient().getTaskStatuses(taskIdsToRequest);
			for (Object obj : wl.list)
			{
				PendingTask pt = (PendingTask) obj;
				for (int i = 0; i < taskIdsToRequest.size(); i++)
				{
					if (taskIdsToRequest.get(i).equals(pt.taskId))
					{
						pt.update(tasks.get(i));
					}
				}
			}
		}

		return new WebModel(wl).addObject(filter).getModelAndView();
	}

	private CalculationClient getClient()
	{
		CalculationClient client = new CalculationClient("Pending-Tasks");
		client.setDeepSleepTime(0);

		return client;
	}

	public ModelAndView tasks(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		AccessChecker.requestRegisteredUserPrivileges();
		//if (true)
		//throw new UserFriendlyException("The modeling section is under maintenance now and will be available again in about a day. Sorry for the inconvenience.");
		PendingTaskFilter filter = getTaskFilter();
		return new WebModel().addObject(filter)
				.setTemplate("browsers/pending-tasks-browser")
				.getModelAndView();
	}

	public ModelAndView published(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		return new WebModel()
				.setTemplate("browsers/published-tasks-browser")
				.getModelAndView();
	}

	class BatchRecalculator extends WrapperThread
	{
		public PendingTaskFilter filter;

		@SuppressWarnings("unchecked")
		@Override
		public void wrapped() throws Exception 
		{
			String statuses[] ={Task.KILL,Task.ERROR};
			for(String status:statuses) {
				filter.status = status;
				setStatus("fetching the tasks...");
				List<PendingTask> pTasks = filter.createCriteria().list();
				int i = 0;
				for (PendingTask pTask : pTasks) 
				{
					setStatus("Posting task " + (++i) + " out of " + pTasks.size());
					pTask = (PendingTask) Globals.session().get(PendingTask.class, pTask.id);
					recalculateTask(pTask, false);
				}
			}
			setStatus("Finished");
		}
	}

	class BatchDeleter extends WrapperThread
	{
		public PendingTaskFilter filter;

		@SuppressWarnings("unchecked")
		@Override
		public void wrapped() throws Exception 
		{
			setStatus("Deleting the tasks...");
			List<Long> pTaskIDs = filter.createCriteria().setProjection(Projections.groupProperty("id")).list();
			int i = 0;
			int total = pTaskIDs.size();
			while (!pTaskIDs.isEmpty())
			{
				List<Long> batchIDs = pTaskIDs.subList(0, Math.min(100, pTaskIDs.size()));
				List<PendingTask> pTasks = (List<PendingTask>) Globals.session().createCriteria(PendingTask.class).add(Restrictions.in("id", batchIDs)).list();
				for (PendingTask pTask : pTasks) 
				{
					if (cancelRequested)
						break;
					setStatus("Deleting task " + (++i) + " out of " + total);
					pTask = (PendingTask) Globals.session().get(PendingTask.class, pTask.id);
					if (pTask != null)
						PendingTaskPeer.deletePendingTask(pTask, true);
				}

				batchIDs.clear();

				Globals.restartAllTransactions(true);
			}

			setStatus("Finished");
		}
	}

	class BatchFetcher extends WrapperThread
	{
		public PendingTaskFilter filter;

		@SuppressWarnings("unchecked")
		@Override
		public void wrapped() throws Exception 
		{
			filter.status = "ready";
			filter.taskType = TaskType.MODEL_TRAINING;
			setStatus("Fetching the tasks...");
			List<PendingTask> pTasks = filter.createCriteria().list();
			int i = 0;
			for (PendingTask pTask : pTasks) 
			{
				pTask = PendingTask.getById(pTask.id); // refetch from DB
				if (pTask.model.taskId != null && !pTask.model.isStatisticsCalculated)
				{
					setStatus("Fetching task " + (++i));
					ModelProcessor processor = ModelFactory.getProcessor(pTask.model);
					processor.onTaskReceived(getClient().getTask(pTask.taskId));
					processor.saveModel();
					if (processor.model.readyModelAttachment.getDataLength() > QSPRConstants.MAXMODELSIZE)
					{
						processor.model.readyModelAttachment.setObject(new ReadyModelAttachment(), AttachmentType.MARSHALABLE);
						logger.info("Warning: the model is too big and will not be saved.");
					}
					Globals.restartAllTransactions(true);
				}
			}

			logger.info("" + i + " ready models have been fetched");

			setStatus("Finished");
		}
	}

	private static Logger logger = LogManager.getLogger(PendingTask.class);
}


