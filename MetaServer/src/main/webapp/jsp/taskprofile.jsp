<%@page import="java.io.StringWriter"%>
<%@page import="qspr.metaserver.frontend.TaskTreeNode"%>
<%@page import="qspr.metaserver.protocol.Task.TaskPriority"%>
<%@ page import="qspr.metaserver.*,qspr.metaserver.protocol.*,java.util.*" %>
<%
	TaskTreeNode root = (TaskTreeNode) request.getAttribute("rootNode");
%>
<html>
	<head>
		<script language="javascript" src="js/jquery-1.4.1.min.js"></script>
		<script language="javascript" src="js/jquery.flot.js"></script>
		<script language="javascript" src="js/plotting.js"></script>
		<script language="javascript" src="js/tooltip.js"></script>
		 <link rel="stylesheet" type="text/css" href="css/tooltip.css" /> 
		<style type="text/css">
			@import url("css/main.css");
			TABLE.nodes TD {border: 3px solid white;}
			.node {margin-left: 10px; clear: both; padding-left: 10px; font-family: Arial;}
			.node.ready TD:first-child {background: url(img/check-16.gif); background-repeat: no-repeat;}
			.node.assigned TD:first-child {background: url(img/inprogress-16.gif); background-repeat: no-repeat;}
			.node.ready {color: #060; background-color: #F0FFF0;}
			.node.level-0 > TD:first-child {padding-left: 20px;}
			.node.level-1 > TD:first-child {padding-left: 50px;}
			.node.level-2 > TD:first-child {padding-left: 60px;}
			
			.level-0 {}
			.level-1 {font-size: 90%;}
			.level-2 {font-size: 80%;}
		</style>
	</head>
	<body>
		<h1>Task profile for <%=root.task.id %></h1>
		<table class="nodes">
		<%=printNode(root, 0) %>
		</table>
	</body>
</html>

<%!

public String printNode(TaskTreeNode node, int level)
{
	StringWriter sw = new StringWriter();
	String ownTime = MetaServer.niceTime(node.task.getCalculationTime());
	String totalTime = MetaServer.niceTime(node.totalTime);
	sw.append("<tr class='node "+node.task.status+" level-"+level+"'><td>" + node.task.id + "</td><td>"+node.task.datarows+" rows</td><td>" + (node.task.calcServerId != null ? node.task.calcServerId : "") + "</td><td>" + node.task.taskType + "</td><td align='center'>" 
	+ node.task.status + "</td><td>time " + ownTime + (!ownTime.equals(totalTime) ? "<br/>(" + totalTime + " including sub-tasks)" : "") + "</td><td>" + (node.task.getDetailedStatus() != null ? node.task.getDetailedStatus() : "")  + "</td></tr>");
	for (TaskTreeNode child : node.children)
		sw.append(printNode(child, level + 1));
	
	return sw.toString();
}

%>