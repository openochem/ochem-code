<%@page import="java.io.StringWriter"%>
<%@page import="qspr.metaserver.frontend.TaskTreeNode"%>
<%@page import="qspr.metaserver.protocol.Task.TaskPriority"%>
<%@page import="qspr.metaserver.*,qspr.metaserver.protocol.*,java.util.*" %>
<%
	Map<Integer, List<Map<String, String>>> filest  = (Map<Integer, List<Map<String, String>>>)request.getAttribute("files_tasks");
	Map<String, List<Map<String, String>>> filess  = (Map<String, List<Map<String, String>>>)request.getAttribute("files_servers");
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
			.delete {color:#A00;}
			.download {color:#777;}
			table th {
				background-color:#DDD;
				align: center;
				padding: 5px;
			}
			table td {
				padding: 5px;
			}
		</style>
	</head>
	<body>
		<h1>Debug files</h1>
		<a class="fb-button" href="/" >Go back to Metaserver status</a>
		<a class='fb-button' href='?action=debugfiles&saction=clear'>Delete all debug files</a><br/><br/>
		<table>
		<tr><th></th><th>Filename</th><th>Size</th><th>Task ID</th><th>Parent Task ID</th><th>Server</th><th>Upload date</th></tr>
		<%
		List<String> servers = new ArrayList<String>();
		servers.addAll(filess.keySet());
		Collections.sort(servers);
		for (String server_name : servers)
			out.print(printServerFileList(server_name, filess.get(server_name)));

		List<Integer> tasks = new ArrayList<Integer>();
		tasks.addAll(filest.keySet());
		Collections.sort(tasks);
		for (Integer task_id : tasks)
			out.print(printTaskFileList(task_id, filest.get(task_id)));
		%>
		</table>
	</body>
</html>

<%!
public String printTaskFileList(Integer task_id, List<Map<String, String>> taskList)
{
	StringWriter sw = new StringWriter();
	sw.append("<tr><td colspan='7'><h2>Task "+task_id+"</h2></td></tr>");
	Collections.sort(taskList, new FileComparator());
	for (Map<String, String> file : taskList)
		sw.append(printSingleFile(file));
	return sw.toString();
}

public String printServerFileList(String server_name, List<Map<String, String>> taskList)
{
	StringWriter sw = new StringWriter();
	sw.append("<tr><td colspan='7'><h2>Server "+server_name+"</h2></td></tr>");
	Collections.sort(taskList, new FileComparator());
	for (Map<String, String> file : taskList)
		sw.append(printSingleFile(file));
	return sw.toString();
}

public String printSingleFile(Map<String, String> file)
{
	return  String.format("<tr><td><a class='delete' href='/?action=debugfiles&saction=delete&id=%s'>[x]</a></td><td><a class='download' href='/?action=debugfiles&saction=download&id=%s'>%s</a></td><td>%s kb</td><td> %s</td><td>%s</td><td>%s (%s)</td><td>%s</td></tr>", 
			file.get("_id"), file.get("_id"), file.get("filename"), (Integer.valueOf(file.get("size")) / 1024)+"", 
			file.get("task_id") == null ? "none" : file.get("task_id"), 
			file.get("parent_task_id") == null ? "none" : file.get("parent_task_id"), 
			file.get("server_name") == null ? "none" : file.get("server_name"),
			file.get("server_address") == null ? "none" : file.get("server_address"),
			file.get("date"));
}

class FileComparator implements Comparator<Map<String, String>>
{
	public int compare(Map<String, String> a, Map<String, String> b)
	{
		String name1 = a.get("filename");
		String name2 = b.get("filename");
		return name1.compareTo(name2);
	}
}
%>