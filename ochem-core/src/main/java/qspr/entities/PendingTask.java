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

package qspr.entities;

import java.io.IOException;
import java.io.Serializable;
import java.sql.Timestamp;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.List;

import javax.persistence.CascadeType;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.EnumType;
import javax.persistence.Enumerated;
import javax.persistence.FetchType;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.Transient;
import javax.servlet.http.HttpServletRequest;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.MarshallingOption;
import qspr.business.Privileges;
import qspr.entities.Attachment.AttachmentType;
import qspr.metaserver.CalculationClient;
import qspr.metaserver.protocol.Command;
import qspr.metaserver.protocol.Task;
import qspr.util.AccessChecker;

import com.eadmet.exceptions.UserFriendlyException;

// An entity for all kinds of pending tasks (model training, descriptor calculation, model application, DM calculation, etc.)
// Midnighter

@Entity
@XmlRootElement(name = "pending-task")
public class PendingTask
{
	@Id
	@Column(name = "pending_task_id")
	@GeneratedValue
	public Long id;

	@ManyToOne
	@JoinColumn(name = "session_id")
	public Session session; // who submitted this task

	@Column(name = "task_id")
	public Integer taskId; // the task ID on metaserver
	public String status;

	@Column(name = "detailed_status")
	public String detailedStatus;

	@Column(name = "time_posted")
	@XmlTransient
	public Timestamp timePosted;

	@Column(name = "time_assigned")
	@XmlTransient
	public Timestamp timeAssigned;

	@Column(name = "time_completed")
	@XmlTransient
	public Timestamp timeCompleted;

	@ManyToOne
	@JoinColumn(name = "model_id")
	public Model model; // A model that is being trained or applied (otherwise empty)

	@Column@XmlElement
	private Integer priority;

	@Column@XmlElement
	public boolean published;

	/**
	 * Contains "published" pending task, i.e. when attachment is tosred in the local datatabase
	 */

	@XmlTransient
	@ManyToOne(cascade={CascadeType.ALL}, fetch = FetchType.EAGER)
	@JoinColumn(name = "ready_task")
	public Attachment<Task> readyTask; 

	@ManyToOne(fetch = FetchType.EAGER)
	@JoinColumn(name = "article_id")
	@XmlTransient
	public Article article;

	@Column(name="calc_server_id")
	@XmlElement
	public String calcServerId;

	@Column
	@XmlElement
	public String name;

	/**
	 * A detailed textual description of the task 
	 */
	@Column
	@XmlElement
	public String description;

	/**
	 * A textual description of the used compound set
	 */
	@Column(name="set_description")
	@XmlElement
	public String setDescription;

	@XmlTransient
	@ManyToOne(cascade={CascadeType.ALL}, fetch = FetchType.EAGER)
	@JoinColumn(name = "attachment_id")
	public Attachment<Serializable> attachment; // A task-specific data necessary to recalculate/fetch the task

	@Enumerated(EnumType.STRING)
	@Column
	public TaskType type;

	public enum TaskType {
		MODEL_TRAINING, 
		MODEL_APPLICATION,
		DESCRIPTOR_CALCULATION,
		TOXALERT_SCREENING,
		SET_COMPARISON,
	};

	public String toString(){
		return (taskId != null? "taskId="+taskId:"") + (name !=null?name:"");
	}

	private void setPostedTask(Integer taskId)
	{
		timePosted = new Timestamp(Calendar.getInstance().getTimeInMillis());
		this.taskId = taskId;
		if (taskId == null)
			throw new UserFriendlyException("The task has not been posted and cannot be stored in pending tasks.");
	}

	@XmlElement
	protected Article getArticle()
	{
		if (!Globals.getMarshallingOption(MarshallingOption.ARTICLE_PENDING_TASKS))
			return article;
		return null;
	}

	@XmlElement(name = "toBeDeleted")
	private Boolean isToBeDeleted() {
		if (!published && model != null && (model.deleteAfterDays != null && !model.published))
			return true;
		return null;
	}

	public void update(Task task)
	{
		if(Task.ERROR.equals(status) && status.equals(task.status)){ // We will keep previous detailed status of error even if task is deleted from the metaserver
			String stat = task.getDetailedStatus();
			if(stat != null && !stat.startsWith(Command.UNKNOWN_TASK))
				detailedStatus = task.getDetailedStatus();
		}else{
			status = task.status;
			detailedStatus = task.getDetailedStatus();
		}

		if (task.timeAssigned != null)
			timeAssigned = task.timeAssigned;
		if (task.timeCompleted != null)
			timeCompleted = task.timeCompleted;
		calcServerId = task.calcServerId;
	}

	public Task retrieveTask(boolean tolerateMetaserverDown) throws ClassNotFoundException, IOException 
	{
		// Did we store the received task locally?
		if (readyTask != null && readyTask.getObject() != null)
			return readyTask.getObject();

		Task task = retrieveTask(taskId,tolerateMetaserverDown);

		if (task != null)
			status = task.status;
		return task;
	}

	public static Task retrieveTask(Integer taskId, boolean tolerateMetaserverDown) throws ClassNotFoundException, IOException{
		CalculationClient client = new CalculationClient("Pending-Tasks");
		if(tolerateMetaserverDown)
			client.setTolerateMetaserverDown();

		return client.getTask(taskId);
	}

	public PendingTask()
	{

	}

	/**
	 * Creates Task without attachment and without indicating task ID, to be indicated later 
	 * @param type
	 */

	public PendingTask(TaskType type)
	{
		this.type = type;
		session = Globals.userSession();
		setPostedTask(0);
	}

	public PendingTask(TaskType type, Integer taskId, Serializable attachedObject)
	{
		this.type = type;
		session = Globals.userSession();
		setPostedTask(taskId);
		if (attachedObject != null)
			setAttachment(attachedObject);
	}

	public void setAttachment(Serializable attachedObject)
	{
		attachment = new Attachment<Serializable>(attachedObject, AttachmentType.SERIALIZABLE, AttachmentSource.PendingTask);
	}

	public PendingTask(TaskType type, Integer taskId)
	{
		this(type, taskId, null);
	}

	public PendingTask setModel(Model model)
	{
		this.model = model;
		return this;
	}

	public PendingTask setPriority(int priority)
	{
		priority = Math.min(priority, AccessChecker.getMaximumTaskPriority(Globals.getCurrentUser()));
		this.priority = priority;
		return this;
	}

	public PendingTask setDescription(String desc)
	{
		this.setDescription = desc;
		return this;
	}

	public Integer getPriority()
	{
		return priority;
	}

	public static PendingTask getByTaskId(int taskId)
	{
		return (PendingTask) Globals.session().createCriteria(PendingTask.class)
				.add(Restrictions.eq("taskId", taskId))
				.uniqueResult();
	}

	public static PendingTask getById(long pTaskID)
	{
		return (PendingTask) Globals.session().get(PendingTask.class, pTaskID);
	}

	@SuppressWarnings("unchecked")
	public static List<PendingTask> getByModel(Model model, TaskType type)
	{
		return Globals.session().createCriteria(PendingTask.class)
				.add(Restrictions.eq("model", model))
				.add(Restrictions.eq("type", type))
				.list();
	}

	public void updatePostedTime()
	{
		timePosted = new Timestamp(Calendar.getInstance().getTimeInMillis());
	}

	public boolean isActive()
	{
		return Task.isAliveStatus(status);
	}

	public boolean isReady()
	{
		return Task.READY.equals(status);
	}


	@Transient
	@XmlAttribute(name="timePosted")
	protected String getTimePosted()
	{
		DateFormat df = new SimpleDateFormat ("yyyy-MM-dd HH:mm:ss");
		if (timePosted != null)
			return df.format(timePosted);
		else
			return null;
	}

	public Privileges getPrivileges(HttpServletRequest request)
	{
		Privileges privileges = new Privileges("pending task");
		User user = Globals.userSession().user;

		if(type == TaskType.MODEL_TRAINING && model != null)return model.getPrivileges(request);
		
		// Only the owner of the a model can edit it
		privileges.canEdit = session == Globals.userSession()
				|| (session.user != null && session.user.equals(user));
		
		privileges.canView = privileges.canEdit || published || Globals.isOCHEMDeveloper() 
				|| AccessChecker.isFromGroup(session.user, Globals.userSession().user);

		return privileges;
	}

	public String getDetailedStatus() {
		return detailedStatus;
	}
}
