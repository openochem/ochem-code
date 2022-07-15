<%@ page import="qspr.metaserver.protocol.Task.TaskPriority,qspr.metaserver.*,qspr.metaserver.protocol.*,java.util.*,java.text.DecimalFormat"%><%
	MetaServer ms = (MetaServer) request.getSession().getAttribute("metaserver");
	DecimalFormat fmt = new DecimalFormat("#.#");
%>


<html>
	<head>
		<style type="text/css">
			@import url("css/main.css");
			BODY {font-family: Courier; padding: 30px 30px 30px 30px;}
			H1 {font-family: Verdana; font-size: 16px; background-color: #0037a8; color: #FF0; padding: 3px;}
			H1 I {color: #88FF88; font-weight: normal;}
			H2 {font-family: Verdana; font-size: 12px; margin-top: 20px; background: #FFC; padding: 2px;}
			B {color: #008;}
			B.failure {color: #800;}
			I {font-weight: bold; font-style: normal; color: #0A0;}
			TABLE.online-servers TD {text-align: left; padding: 5px 10px 5px 10px; background-color: #FFA;}
			TABLE.online-servers TR.disabled TD {color: #666 !important; background-color: #CCC !important;}
			TABLE.online-servers TR.busy TD {background-color: #FFEEAA !important;}
			.disabled B, .disabled I, .disabled {color: #666;}
			SMALL {font-size: 8px;}
			A IMG {border: none; vertical-align: middle; margin-right: 5px;}
			.red {color: red !important;}
			.invisible {visibility: hidden; display: none;}
			.upper-bar {position: fixed; right: 10px; top: 5px; }
			.search {border: 1px solid gray;}
			.offline, .offline I, .offline B, .offline TD {color: gray !important;}
			.configuration-xml {position: absolute; padding: 5px 30px 30px 30px; background-color: #DDDDDD; border: 1px solid black;}
			A.close-details {float: right; text-decoration: none; border: 1px solid #AAAAAA; position: relative; top: -4px; left: 28px; font-size: 30px; color: #500;}
			
			.grey TD, .grey TD B, .grey TD * {
			    background-color: #EEEEEE !important;
			    color: #999999 !important;
			}
		</style>
		<script language="javascript" src="js/jquery-1.4.1.min.js"></script>
		<script language="javascript">
			$(document).ready(function(){
				var search = $("#search");
				search.keyup(function(){
					$("[key]").addClass("invisible");
					$("[key*="+$(this).val()+"]").removeClass("invisible");
				});

				if (search.val())
					search.keyup();


				$(".close-details").click(function(){
					$(this).parents(".configuration-xml").addClass("invisible");
					return false;
				});	

				$(".open-details").click(function(){
					$(".configuration-xml").addClass("invisible");
					$(this).parents("td").find(".configuration-xml").removeClass("invisible");
					return false;
				});
			});
		</script>
	</head>
	<body>
		<div class="upper-bar">
			<a href="javascript:location.reload(true);"><img src="img/refresh.gif"/></a>
			<input type="text" class="search" id="search"></input>
		</div>
		<h1>Metaserver status <i>(<%=ms.connections.get() %> connections, uptime <%=ms.niceTime(ms.getUptime()) %>)</i></h1>
		<a href="?action=charts" class="fb-button">Charts</a>
		<a href="?action=serversAvailability" class="fb-button">Servers availability report</a>
		<a href="?action=debugfiles" class="fb-button">Debug files</a>
		<a href="?action=resettasks" class="fb-button">Reset tasks</a>
		<%
			if (ms.stopTaskAssignment)
			{
				%>
				<div class="warning">
					The MongoDB storage has exceeded its capacity. Assigning new tasks has been suspended for a while.<br/>
					No action is required: the task assignment will resume automatically when MongoDB is cleaned up. 
				</div>
				<%		
			}
		
			for (String taskType : ms.forcedTaskBindings.keySet())
			{
				%><br/><b><%=taskType%></b> tasks are bounded to the server <b><%=ms.forcedTaskBindings.get(taskType)%></b> <a href="?action=removeTaskBinding&taskType=<%=taskType%>">delete binding</a> <% 
			}
		%>
		<h2>Queued tasks in cache (<%=ms.getQueuedTasksCount() %> total, <%=ms.getPrimaryQueuedTasksCount() %> root tasks)</h2>
		<%
			synchronized (ms.taskQueue)
			{
				List<String> keys = new ArrayList<String>();
				keys.addAll(ms.taskQueue.tasksQueue.keySet());
				Collections.sort(keys);
				for (String key : keys)
				{ 
					List<Task> tasks = ms.taskQueue.tasksQueue.get(key);
					if (tasks == null)
						continue;
					for (Task task : tasks)
					{
					     %><b title="<%=task.datarows%> rows, Priority <%=task.priority%>"><%=task%></b> in queue for <%=ms.niceTime((Calendar.getInstance().getTimeInMillis() - task.time.getTime())/1000)%>
					    <% if (task.minRequiredMemory != null){ %>, requires <%=fmt.format(task.minRequiredMemory / 1024D)%>GB of memory <%} %>
					    <% if (task.getPreferredServer() != null){ %>, exclusively for <%=task.getPreferredServer()%><%} %>
					    <br/><%
					 }
				}
			}
		%>
		<h2>Assigned tasks pending (<%=ms.assignedTasks.size() %>)</h2>
		<table border="0" cellspacing="2" cellpadding="1">
		<%
			synchronized(ms.assignedTasks)
			{
				List<String> keys = new ArrayList<String>();
				keys.addAll(ms.assignedTasks.keySet());
				Collections.sort(keys);
				
				for (String key : keys)
				{
				    Task task = ms.assignedTasks.get(key);
				    OnlinePeer server = ms.onlinePeers.get(key);
				    boolean offline = server == null || server.getPing() > MetaServer.METASERVER_TASK_NOT_RESPONING;
				    %>
				    <tr key="<%=task%>.<%=key%>" <% if (offline) {%>class="offline"<%} %>>
				    <td>
					    <% if (task.parentTaskId == null) { %>
					    	<a href="?killtask=<%=task.id%>" title="Kill task"><img src="img/icon_axe.gif"/></a>
					    <% } else { %>
				    	&nbsp;&nbsp;&nbsp;&nbsp;<a href="?restarttask=<%=task.id%>" title="Restart task"><img src="img/refresh.gif"/></a>
				    	<a href="?stoptask=<%=task.id%>" title="Request task to stop"><img src="img/check-16.gif"/></a>
				    	<% }
					    %>
				    </td>
				    <td><i><%=key%></i></td>
				    <td> calculates</td>
				    <td> 
				     <% if (task.parentTaskId == null) { %>
				     	<a href="?action=taskProfile&id=<%=task.id %>" target="_blank"><b title="<%=task.datarows %> rows, priority <%=task.priority%>"><%=task%></b></a>
				     <% } else { %>
				     	<b title="<%=task.datarows %> rows, priority <%=task.priority%>"><%=task%></b>
				     <% } %>
				     <% if (task.taskName != null) {%><br/><small><%=task.taskName %></small><% } %>
				    	
				    </td>
				    <td> for <%=ms.niceTime((Calendar.getInstance().getTimeInMillis() - task.timeAssigned.getTime())/1000)%></td>
				    <td>
				    	<% if (task.parentTaskId == null && task.client != null) { %>
					    	submitted by <i><%=task.client %></i> 
					    <% }  %>
				    </td>
				    </tr>
				    <%
				}
			}
		%>
		</table>
		<%
		// Calculate the servers (total, free, by priority)
		int totalServers = 0;
		int freeServers = 0;
		int[] freeByPriority = new int[]{0, 0, 0, 0};
		int notResponding = 0;
		int busy = 0;
		synchronized(ms.onlinePeers)
		{
			for (String key : ms.onlinePeers.keySet())
			{
				OnlinePeer peer = ms.onlinePeers.get(key);
				if (peer.isClient || peer.serverInfo == null)
					continue;
				
				if (peer.notResponding() && peer.currentTask == null)
					continue;
				
				totalServers++;
				if (peer.currentTask == null)
				{
					freeServers++;
					if (peer.serverInfo.minimumPriority == null || peer.serverInfo.minimumPriority <= TaskPriority.LOW)
						freeByPriority[0]++;
					else if (peer.serverInfo.minimumPriority <= TaskPriority.NORMAL)
						freeByPriority[1]++;
					else if (peer.serverInfo.minimumPriority <= TaskPriority.HIGH)
						freeByPriority[2]++;
					else
						freeByPriority[3]++;
				}
				else if (peer.notResponding())
					notResponding++;
				else
					busy++;
			}
		}
		%>
		<h2>Online servers 
			(<%=totalServers %> total, 
			<% if (busy > 0) { %><%=busy %> busy, <%} %>
			<% if (notResponding > 0) { %><%=notResponding %> not responding, <%} %>
			<%=freeServers %> free = <%=freeByPriority[0] %> low, <%=freeByPriority[1] %> normal + <%=freeByPriority[2] %> high + <%=freeByPriority[3] %> extra-high priority servers)
		</h2>
		<table class="online-servers">
		<%
			long currentTime = Calendar.getInstance().getTimeInMillis();
			synchronized(ms.onlinePeers)
			{
				List<String> keys = new ArrayList<String>();
				keys.addAll(ms.onlinePeers.keySet());
				Collections.sort(keys);
				for (String key : keys)
				{ 
				    OnlinePeer peer = ms.onlinePeers.get(key);
				    if (peer.notResponding() && peer.currentTask == null)
				    	continue;
				    if (!peer.isClient && peer.serverInfo != null) 
				    {
				    %><tr<% if (peer.disabled) {%> class="disabled"<%} else if (peer.currentTask != null) {%> class="busy"<%} else if (peer.notResponding()) {%> class="grey"<%}  %> key="<%=key%>.<%=peer.currentTask%>.<%=peer.serverInfo.owner %>">
				    	<td>
				    		<i title="<%=peer.serverInfo.owner %>"><%=key%>
					    	<% 
					    	if (peer.serverInfo.minimumPriority != null) 
					    	{%>
					    		<% if (peer.serverInfo.minimumPriority >= 10) {%>
					    			<img src="img/priority2.png" title="Accepts only extra high priority tasks"/>
					    		<% } else if (peer.serverInfo.minimumPriority >= 1) {%>
					    			<img src="img/priority1.png" title="Accepts only high priority tasks"/>
					    		<%}
					    	}%>
				    		<br/>
				    		<small>
				    			<%=peer.serverInfo.ipAddress %><br/><%=peer.serverInfo.platform %>, <%=fmt.format(peer.serverInfo.availableMemory / 1024D) %>GB<br/>
				    			<% if (peer.serverInfo.version == null || ms.currentVersion == null) { %>
				    				<span class="red">unknown version</span>
				    			<% } else if (!peer.serverInfo.version.equals(ms.currentVersion)) { %>
				    				<span class="red"><%=peer.serverInfo.version %></span>
				    			<% } else { %>
				    				<%=peer.serverInfo.version %>
				    			<% } %>
				    			
				    			<% if (peer.serverInfo.configurationXml != null) { %>
									<a href="#" class="open-details" title="Show version.xml">&gt;&gt;</a>
								<% } %>
				    		</small>
				    		</i>
				    		<div class="invisible configuration-xml">
				    			<a href="#" class="close-details">x</a>
				    			<pre><% if (peer.serverInfo.configurationXml != null) { %><%=ms.htmlAngleBrackets(peer.serverInfo.configurationXml) %><% } %>
				    			</pre>
				    		</div>
				    	</td>
				    	<td><%=(currentTime - peer.lastActive)/1000%> sec</td>
				    	<td>
				    		<% if (currentTime  - peer.conflictTime < 2000) { %>
				    			<img src="img/warning.gif"/>
				    		<% } %>
				    		<b title="Supported task types"><%=peer.serverInfo.supportedTaskTypes %></b>
				    		<% if (peer.serverInfo.failures != null) { %><br/><b class="failure" title="Disabled task types"><%=peer.serverInfo.failures %></b><% } %>
				    		<% if (peer.disabledTasks.size() != 0) { %><br/><b class="failure" title="Tasks, disabled manually for this server">Disabled manually: <%=peer.disabledTasks.toString() %></b><% } %>
				    	</td>
				    	<td>
				    		<nobr>CPU load: <%=String.format("%.1f", 1.0 * peer.serverInfo.cpuUsage*100) %>%</nobr><br/>
				    		<nobr>Disc space: <%=String.format("%.1f", 1.0 * peer.serverInfo.diskSpaceLeft / 1024) %>GB</nobr><br/>
				    		<nobr>Used RAM: <%=String.format("%.1f", 1.0 * peer.serverInfo.usedMemory) %>MB</nobr><br/>
				    		<nobr>Peak RAM: <%=String.format("%.1f", 1.0 * peer.serverInfo.peakMemory) %>MB</nobr><br/>
				    	</td>
				    	<td>
				    		<% if (peer.currentTask != null) { %><b><%=peer.currentTask%></b><br/><% } %>
				    		<%=peer.status %>
				    		<% if (peer.currentTask == null) {%><a href="?server=<%=key%>&restart=1" title="Restart and re-run tests"><img src="img/restart.png"/></a><%} %>
				    	</td>
				    	</tr><%
				    }
				}
			}
		%>
		</table>
		<h2>Online clients</h2>
		<table class="online-servers">
		<%
			synchronized(ms.onlinePeers)
			{
				List<String> keys = new ArrayList<String>();
				keys.addAll(ms.onlinePeers.keySet());
				Collections.sort(keys);
				for (String key : keys)
				{ 
				    OnlinePeer peer = ms.onlinePeers.get(key);
				    if (peer.isClient) 
				    {
				    %><tr key="<%=peer.currentTask %>">
				    	<td><i><%=key%><br/><small><%=peer.ipAddress %></small></i></td>
				    	<td><%=(Calendar.getInstance().getTimeInMillis() - peer.lastActive)/1000%> sec</td>
				    	<td><% if (peer.currentTask != null) { %><b><%=peer.currentTask%></b>:<% } %><%=peer.status %></td>
				    	</tr><%
				    }
				}
			}
		%>
		</table>
	</body>
</html>