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

package qspr.metaserver.protocol;

import java.io.IOException;
import java.io.Serializable;
import java.sql.Timestamp;
import java.util.Calendar;
import java.util.regex.Pattern;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.GenerationType;
import javax.persistence.Id;
import javax.persistence.Transient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.annotations.GenericGenerator;

import qspr.exceptions.CalculationException;
import qspr.metaserver.MemoryEstimate;
import qspr.metaserver.transport.DataReference;
import qspr.metaserver.transport.DataReferenceFactory;
import qspr.metaserver.transport.DataReferencer;
import qspr.util.ClassCompressor;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.MAILERConstants;
import com.eadmet.utils.OCHEMUtils;

/**
 * A low-level representation of the metaserver-managed calculation task.
 * The tasks are transferred between the metaserver and its clients in a Java-serializable format
 * 
 * The class contains the comprehensive summary of a task:
 * data, configuration, task type, priority, status, etc.
 * 
 * @author midnighter
 */
@Entity
public class Task implements Serializable
{
	public static final long serialVersionUID = 2;

	static transient Logger logger = LogManager.getLogger(Task.class);

	public enum REFERENCE_TYPE {CFG, DATA, RESULT};

	//	public enum STATUS {INIT,ASSIGNED,KILL};
	public static final String INIT="init", ASSIGNED="assigned", KILL="kill", ERROR="error", READY = "ready", STOP = "stop", FETCHED = "fetched", KILLED="killed";

	public static final String Workflow = "Workflow";
	public static final String TASK_DATABASE = "metaserver";
	public static final int REFERENCE_LENGTH = 24;
	public static final int MAX_STATUS_LENGTH = 10240;
	public static final String TEST_USER_PREFIX = "Test_";
	public static final String EMAIL_OCHEM = MAILERConstants.EMAIL_OCHEM;

	@Column(name = "task_id", updatable = false, nullable = false)
	@GeneratedValue(strategy = GenerationType.AUTO, generator = "native")
	@GenericGenerator(name = "native", strategy = "native")
	@Id
	public Integer id;

	@Column(name = "reference_id")
	public String referenceId;

	@Column(name = "parent_task_id")
	public Integer parentTaskId;

	@Column
	public String status=null; // TODO Midnighter: put an enumerated type here

	@Column
	public String client;

	@Column
	public Integer datarows;

	@Column
	public Integer cols;

	@Column
	public Integer nonzero;

	@Column(name = "detailed_status")
	protected String detailedStatus;

	@Column(name="task_type")
	public String taskType;

	@Column(name = "calc_server_id")
	public String calcServerId;

	@Column
	public Timestamp time;

	@Column(name = "time_assigned")
	public Timestamp timeAssigned;

	@Column(name = "time_completed")
	public Timestamp timeCompleted;

	@Column
	public double priority = TaskPriority.NORMAL;

	@Column(name = "peak_memory_usage")
	public Integer peakMemoryUsage;

	@Column(name = "task_name")
	public String taskName;

	/**
	 * These columns contain data or references to data
	 */

	@Column
	protected byte[] configuration;

	@Column
	protected byte[] data;

	@Column
	protected byte[] result;

	/*
	 * These three variables are used as proxy for local Task to avoid compression/decompression and decrease memory usage for Results
	 */
	@Transient
	private Serializable localResult, localData, localConfig;

	@Column(name = "last_access")
	public Timestamp lastAccess;

	@Column(name = "scheduled_kill")
	public boolean scheduledKill;

	// In mbytes
	@Column(name = "min_required_memory")
	public Integer minRequiredMemory;

	// Number of times task was resubmitted
	@Column(name = "resubmitted")
	public Integer resubmitted;

	/**
	 * How many clients are interested in this task?
	 * 0 - the task in non-cachable
	 */
	@Column(name = "ref_count")
	public int referenceCount;

	/**
	 * MD5 key of the task input (data + configuration).
	 * Used for the task-level caching
	 */
	@Column(name = "task_md5")
	public String md5;

	@Column(name = "preferred_server")
	public String preferredServer;

	/**
	 * Indicates that this is a debug task and extra logging information should be returned by the calculation servers.
	 * 0 - no debug
	 * 1 - debug errors
	 * 2 - debug all tasks
	 */
	@Column
	public int debug = 0;

	@Transient
	static DataReferencer referencer;

	/**
	 * The user who originated the calculation
	 */
	@Column
	protected String user;

	/**
	 * A flag to disable the task-level cache for this particular task.
	 * Its a serialisable, but not persistent flag.
	 */
	@Transient
	public boolean disableTaskLevelCache;

	@Transient 
	private boolean isLocalTask = false;

	@Transient 
	int tmpSize = 0, cfgSize = 0;

	public Task()
	{
	}

	public Task(String type, Serializable configuration, Serializable data, boolean skipSerialization) throws IOException
	{
		long time = Calendar.getInstance().getTimeInMillis();

		isLocalTask = skipSerialization;	
		this.taskType = type;
		this.setConfiguration(configuration);
		this.setData(data);

		if(isLocalTask)logger.info("Task serialization was not performed");
		else
			logger.info("Task serialization took " + (Calendar.getInstance().getTimeInMillis() - time)/1000 + " sec.");
	}

	public int numberCrashes(){
		return resubmitted == null? 0: resubmitted.intValue();
	}

	public Task(String type, Serializable configuration, Serializable data) throws IOException
	{
		long time = Calendar.getInstance().getTimeInMillis();

		isLocalTask = true; // by default we create a local task...
		this.taskType = type;
		this.setConfiguration(configuration);
		this.setData(data);

		logger.info("Task serialization took " + (Calendar.getInstance().getTimeInMillis() - time)/1000 + " sec.");
	}

	public static Task fromCommand(Command command)
	{
		return (Task) command.data;
	}

	public String getDetailedStatus()
	{
		return detailedStatus;
	}

	public boolean shouldKeepLogs()
	{
		return debug > DebugLevel.NONE && ("error".equals(status) || debug >= DebugLevel.ALL);
	}

	public boolean hasMultipleReferences()
	{
		return referenceCount > 1;
	}

	public void correctTaskPriority(double initialPriority)
	{
		priority = Math.floor(initialPriority);
		long minutes = (System.currentTimeMillis()/1000 - 1551076798)/60; // from time when I started 
		priority += 1./(Math.floor(Math.log(datarows != null && datarows > 10? datarows: 10)) + 1.); // smallest tasks first -- largest impact
		priority += 1.0d / minutes; // older - lower priority - minor impact
	}

	public Task setPriority(double priority)
	{
		this.priority = priority;
		return this;
	}

	public Task setParentTask(Task task)
	{
		parentTaskId = task.id;
		debug = task.debug;
		if (parentTaskId != null && parentTaskId < 0 && task.parentTaskId != null)
			parentTaskId = task.parentTaskId;
		priority = task.priority;

		if (task.getPreferredServer() != null && task.getPreferredServer().startsWith("re:"))
			setPreferredServer(task.getPreferredServer());

		user = task.user;
		return this;
	}

	public long getReferenceSize(){
		return (data == null?0:data.length) 
				+ (result == null?0:result.length) 
				+ (configuration ==null? 0: configuration.length);
	}

	public long getDataSize()
	{
		return getSize(data);
	}

	private long getResultSize(){
		return getSize(result);
	}

	public long getConfigurationSize() 
	{
		long size = getSize(configuration);
		return size > cfgSize ? size : cfgSize;
	}

	private long getSize(byte data[]){
		if(data == null)return 0; // this is localTask, we will not try to determine its size at this time
		Serializable reference = ClassCompressor.byteToObject(data);
		if(referencer == null) referencer = DataReferenceFactory.createReferencer();
		if(reference instanceof DataReference)
			return referencer.getDataSize((DataReference)reference);
		return data.length;
	}

	public long getByteSize()
	{
		return getDataSize() + getResultSize() + getConfigurationSize();
	}

	public Serializable getConfiguration() throws IOException
	{
		if(localConfig != null) return localConfig;
		return resolveReference(configuration);
	}

	public Serializable getResult() throws IOException
	{
		if(localResult != null) return localResult;
		return resolveReference(result);
	}

	public Serializable getData() throws IOException
	{
		if(localData != null) return localData;
		return resolveReference(data);
	}

	private Serializable resolveReference(byte []dat){
		if(dat == null) return null;
		Serializable reference = ClassCompressor.byteToObject(dat);
		if(referencer == null) referencer = DataReferenceFactory.createReferencer();
		if(reference instanceof DataReference)
			return referencer.getReference((DataReference)reference);
		return reference;
	}

	public void keepReferenceInTask(DataReference ref, REFERENCE_TYPE type){
		StringBuffer buf = new StringBuffer(REFERENCE_LENGTH*4);
		buf.append(referenceId == null?"":referenceId);

		int n = 0;
		switch(type){
		case CFG:
			n = 0;
			break;
		case DATA:
			n = 1;
			break;
		case RESULT:
			n = 2;
			break;
		}

		if(buf.length()<n*REFERENCE_LENGTH)
			for(int i = buf.length();i<n*REFERENCE_LENGTH;i++)buf.append(' ');
		if(buf.length() == n*REFERENCE_LENGTH)
			buf.append(ref.getReference());
		else
			buf.replace(n*REFERENCE_LENGTH,(n+1)*REFERENCE_LENGTH, ref.getReference());

		referenceId = buf.toString();
	}

	private byte[] createReference(Serializable obj, REFERENCE_TYPE type) throws IOException{

		byte data[] = ClassCompressor.objectToByte(obj);
		if(obj instanceof DataReference) return data; // work is already done

		tmpSize = data.length;

		if(referencer == null) referencer = DataReferenceFactory.createReferencer();
		if (data != null && data.length > 250){
			DataReference ref = referencer.saveReference(obj, TASK_DATABASE);
			keepReferenceInTask(ref,type);
			return ClassCompressor.objectToByte(ref);
		}
		return data;
	}

	public void setConfiguration(Serializable configuration) throws IOException
	{
		this.configuration = null;
		localConfig = null;

		if(configuration == null) return;

		if(isLocalTask && !(configuration instanceof DataReference))
			localConfig = configuration;
		else
			this.configuration = createReference(configuration, REFERENCE_TYPE.CFG);

		cfgSize = tmpSize; // saving local size

	}

	public void setResult(Serializable result) throws IOException
	{
		this.result = null;
		localResult = null;

		if(result == null)return;

		if(isLocalTask && !(result instanceof DataReference))
			localResult = result;
		else
			this.result = createReference(result, REFERENCE_TYPE.RESULT);
	}

	public void setData(Serializable data) throws IOException
	{
		this.data = null;
		localData = null;

		if(data == null) return;

		setTaskProperties(data); // it is assumed that setData are normally called after setConfiguration()

		if(isLocalTask && !(data instanceof DataReference))
			localData = data;
		else
			this.data = createReference(data, REFERENCE_TYPE.DATA);
	}

	public boolean isReady()
	{
		return isError() || !Task.isAliveStatus(status);
	}

	public boolean isError()
	{
		return status != null && (status.startsWith(KILL)  || status.startsWith(ERROR));
	}

	public void check() throws CalculationException
	{

		if (ERROR.equals(status))
			// Identify whether this is a message (UserFriendlyException thrown by server)
			//  or a severe Exception. That is simplified but working way
			if (detailedStatus.contains("Exception"))
				throw new CalculationException("Calculation server error in "+this+" by "+calcServerId+": \n"+detailedStatus);
			else
				throw new UserFriendlyException(detailedStatus);

		else if (status != null && status.startsWith(KILL))
			throw new CalculationException("Task " +id+ " has been killed at "+calcServerId);
	}

	public String toString()
	{
		return taskType + "/" + id + (parentTaskId != null ? "." + parentTaskId : "");
	}

	public void mergeConfiguration(Mergable mergeConfiguration) throws Exception
	{
		Mergable currentConfiguration = (Mergable) getConfiguration();
		if (mergeConfiguration != null)
			if (currentConfiguration == null)
				this.setConfiguration((Serializable) mergeConfiguration);
			else
			{
				currentConfiguration.mergeWith(mergeConfiguration);
				this.setConfiguration((Serializable) currentConfiguration);
			}
	}

	public long getCalculationTime()
	{
		if (timeAssigned == null)
			return 0;
		if (timeCompleted != null)
			return (timeCompleted.getTime() - timeAssigned.getTime()) / 1000;
		else if (Task.isAssigned(status))
			return (Calendar.getInstance().getTimeInMillis() - timeAssigned.getTime()) / 1000;
		else 
			return 0;
	}

	public boolean canBeCalculatedBy(String serverId) // Is there any preference on who should calculate this task?
	{
		String preferredServer = getPreferredServer();
		if (preferredServer == null)
			return true;
		if (!preferredServer.startsWith("re:"))
			return preferredServer.equals(serverId);
		preferredServer = preferredServer.replaceAll("re:", "");
		Pattern p = Pattern.compile(preferredServer,Pattern.DOTALL);
		return p.matcher(serverId).matches();
	}

	public Task setPreferredServer(String serverId)
	{
		this.preferredServer = serverId;
		return this;
	}

	public String getPreferredServer()
	{
		return preferredServer;
	}

	public boolean hasHigherPriorityThan(Task task)
	{
		// if priorities differ - simply compare them
		if (priority != task.priority)
			return priority > task.priority;

			// if priorities are the same, prefer the task the with earliest-posted parent
			if (task.parentTaskId != null && parentTaskId != null && task.parentTaskId > 0 && parentTaskId > 0)
				return parentTaskId < task.parentTaskId;

			return false;
	}

	public void createReferences() throws IOException
	{
		boolean resolved = !isLocalTask && minRequiredMemory != null;
		logger.info("Creating and resolving references "+resolved);
		if(resolved) return; // task is set and minRequiredMemory is determined
		isLocalTask = false; // it is time to recalculate all sizes before sending the task to metaserver
		setConfiguration(getConfiguration());
		setData(getData());
		setResult(getResult());
	}

	public Task setMinRequiredMemory(int requiredMemory)
	{
		minRequiredMemory = requiredMemory;
		return this;
	}


	private void setTaskProperties(Serializable data){

		if(data == null) return; // to be determined later

		if (data instanceof Dimensionable)
		{
			this.datarows = ((Dimensionable) data).getRows();
			this.cols = ((Dimensionable) data).getCols();
			this.nonzero = 0; //((Dimensionable) data).getMaxNonZero();
		}

		if (data instanceof MemoryEstimate)
			minRequiredMemory = ((MemoryEstimate) data).getMinRequiredMemory(this);
	}

	public int getMinRequiredMemory(){

		if(minRequiredMemory == null)
			try{  // datasizes are not yet set
				setTaskProperties(getData());
			}catch(IOException e){ // not yet ready, will try later
				logger.info("There was an exception in getMinRequiredMemory "+e.getMessage());
				return 0;
			}

		return minRequiredMemory;
	}

	public String getUser() {
		return user;
	}

	public Task setUser(String user) {
		if (user != null && user.length() > 50)
			this.user = user.substring(0, 50);
		else
			this.user = user;

		return this;
	}

	public void setError(String error)
	{
		status = "error";
		detailedStatus = error;
	}

	public boolean isCachable()
	{
		if (disableTaskLevelCache)
			return false;
		return (!Workflow.equals(taskType) && !"MolOptimiser".equalsIgnoreCase(taskType) && (user == null || !user.startsWith(TEST_USER_PREFIX)));
	}

	public static class TaskPriority
	{
		public static int EXTRA_HIGH = 10; // only for testing and really urgent tasks of super users
		public static int HIGH = 2; // Application tasks, descriptor calculation tasks
		public static int NORMAL = 1; // Majority of tasks submitted by users
		public static int LOW = 0; // Batch calculation tasks submitted by CM

		public static int lowerPriority(int priority)
		{
			if (priority > EXTRA_HIGH)
				return EXTRA_HIGH;
			if (priority > HIGH)
				return HIGH;
			if (priority > NORMAL)
				return NORMAL;
			return LOW;
		}
	}

	public static class DebugLevel {
		public static final int NONE = 0;
		public static final int ERRORS = 1;
		public static final int ALL = 2;
	}

	public boolean hasResult() {
		return localResult != null || result != null;
	}

	public String getMD5() {
		return OCHEMUtils.getMD5(data, configuration);
	}

	public boolean isLocalTask() {
		return isLocalTask;
	}

	public void setLocalTask() {
		isLocalTask = true;
	}

	public void clearData() {
		localData = null;
		data = null;
	}

	public void clearConfig() {
		localConfig =  null;
		configuration = null;		
	}

	public void clearResult() {
		localResult =  null;
		result = null;		
	}

	public void clearAll() {
		localData = localConfig = localResult = null;
		data = configuration = result = null;
	}

	public void updateResult(Task task) {
		result = task.result;
		referenceId = task.referenceId;
	}

	public static boolean isAliveStatus(String status) {
		return isAssigned(status) || INIT.equals(status);
	}

	public static boolean isAssigned(String status) {
		return ASSIGNED.equals(status) || STOP.equals(status);
	}

	public String setDetailedStatus(String status)
	{
		if(status != null) {
			status = status.replaceAll(" +", " ");
			status = status.replaceAll("[^:; \\n0-9a-zA-Z\\/\\.\\+\\-\\%]*", "");
		}
		return detailedStatus=status;
	}

	public static void main(String[] args) throws Exception{
		Task t=new Task();
		System.out.println(t.setDetailedStatus(" a§aa:;   \t§BBB"));
	}

}
