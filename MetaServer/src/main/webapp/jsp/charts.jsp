<%@page import="qspr.metaserver.protocol.Task.TaskPriority"%>
<%@ page import="qspr.metaserver.*,qspr.metaserver.protocol.*,java.util.*" %>
<%
	MetaServer ms = (MetaServer) request.getSession().getAttribute("metaserver");
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
			#tasks {width: 800px; height: 600px;}
		</style>
	</head>
	<body>
		<h1>Metaserver Charts</h1>
		<a href="/" class="fb-button">Go back to Metaserver status</a><br/><br/>
		Show last <select name="">
			<option value="60">60 minutes</option>
			<option value="120">2 hours</option>
			<option value="240">4 hours</option>
			<option value="1440">1 day</option>
			<option value="2880">2 days</option>
		</select>
		<div id="tasks">
			
		</div>
		<script language="javascript">
			function load(interval)
			{
				$.ajax({
					url: "?action=getStatistics",
					dataType: "json",
					data: "interval=" + interval,
					success: function(obj)
					{
						var plot = new Plot();
						plot.options.legend.position = "ne";
						plot.options.legend.show = true;
						var dt;
						
						var add = function(name, label, color)
						{
							dt = new Array();
							for (var i = 0; i < obj.logs.length; i++)
								dt.push([i, obj.logs[obj.logs.length - i - 1][name]]);
							plot.addData(dt, color).setLabel(label);
							plot.currentData.lines = {show: true};
							plot.currentData.points.show = false;
							plot.currentData.hoverable = true;
						}
						
						add("completedTasks", "Completed tasks", "green");
						plot.currentData.bars = {show: true};
						plot.currentData.lines.show = false;
						add("newTasks", "New tasks", "blue");
						add("errors", "Errors", "red");
						
						plot.render("#tasks");
						$("#tasks").bind("plothover", function (event, pos, item) {
					        if (item)
					          	tooltip.show("" + item.datapoint[1] + " " + item.series.label, 200);
					        else
					        	tooltip.hide();
					    });
						
						$("#tasks").bind("mouseout", function (event, pos, item) {
					       tooltip.hide();
					    });
					},
					error: function(obj)
					{
						window.alert("Error");
					}
				});
			}
			
			$("select").change(function(){
				load($(this).val());
			});
			
			$("select").change();
		</script>
		<div id="tt">
			<div id="tttop"> </div>
			<div id="ttcont"> </div>
			<div id="ttbot"> </div>
		</div>
	</body>
</html>