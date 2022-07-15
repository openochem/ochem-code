<%@page import="qspr.metaserver.frontend.ServersAvailbilityData.TaskTypePriorityInfo"%>
<%@page import="qspr.metaserver.protocol.Task.TaskPriority"%>
<%@ page import="qspr.metaserver.*,qspr.metaserver.frontend.*,qspr.metaserver.protocol.*,java.util.*" %>
<%
	ServersAvailbilityData data = (ServersAvailbilityData) request.getAttribute("data");
	List<Integer> priorities = new ArrayList<Integer>();
	priorities.addAll(data.priorities);
	Collections.sort(priorities);
	
	List<String> taskTypes = new ArrayList<String>();
	taskTypes.addAll(data.taskTypes);
	Collections.sort(taskTypes);
%>
<html>
	<head>
		 <link rel="stylesheet" type="text/css" href="css/tooltip.css" /> 
		<style type="text/css">
			@import url("css/main.css");
			TABLE.priorities TD, TABLE.priorities TH {
				background-color: #EEF;
				padding: 3px 10px;
				font-family: Helvetica;
			}
			
			TD.free-servers {text-align: center;}
			SPAN.total-servers {color: #555;}
			SPAN.queued-tasks {color: #B33; font-size: 80%;}
		</style>
	</head>
	<body>
		<h1>Availability of calculation servers</h1>
		<a href="/" class="fb-button">Go back to Metaserver status</a><br/><br/>
		A summary of the free servers by task types and priorities:
		<table class="priorities">
			<tr>
				<th>Task type / Priority</th>
				<% for (Integer priority : priorities) {
				%><th><%=priority %></th><% 	
				}
				%>
			</tr>
			<% for (String taskType : taskTypes) {
				%><tr><td><%=taskType %></td><% 
					for (Integer priority : priorities) {
						%><td class='free-servers'>
						<% 
							TaskTypePriorityInfo availability = data.getAvailability(taskType, priority); 
							if (availability.totalServers > 0 || availability.queue > 0) { 
								if (availability.totalServers  != availability.freeServers) { %>
								<%=availability.freeServers %> <span class='total-servers'>(of <%=availability.totalServers %>)</span>
								<% } else { %>
									<%=availability.freeServers %>
								<% }
									if (availability.queue > 0) {
										%><br/><span class='queued-tasks'><%=availability.queue %> tasks</span><% 
									}
								%>
						</td><%	
						}
					}
				%></tr><% 
			} %>
		</table>
	</body>
</html>
