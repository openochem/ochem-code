<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	<xsl:template name="content">
		<script language="javascript" src="js/lib/jquery.flot.js"></script>
		<script language="javascript" src="js/commons/plotting.js"></script>
		<script language="javascript" src="js/lib/tooltip.js"></script>
		 <link rel="stylesheet" type="text/css" href="css/tooltip.css" />
		<style type="text/css">
			#memory, #load{width: 800px; height: 350px;}
			.monitors H1 {font-family: Arial;}
			.monitors > DIV {float: left; width: 800px;}
		</style>
		<table width="100%">
	  		<tr>
	  			<td class="itunes-up">
	  				<img src="img/icons/monitor32.png"/><h1>System monitoring</h1>
	  				This page shows different system monitors (currently, only memory usage)
	  			</td>
	  		</tr>
	  		<tr>
	  			<td class="itunes-right">
	  			Show last <select name="length">
					<option value="600">10 minutes</option>
					<option value="3600" selected="1">1 hour</option>
					<option value="7200">2 hours</option>
					<option value="14400">4 hours</option>
					<option value="172800">2 days</option>
				</select>
				<a href="javascript:load('gc')" class="fb-button">Run garbage collection</a>
				<div class="monitors">
					<div>
					<h1>OCHEM memory usage</h1>
		  			<div id="memory" class="plot1">
				
					</div>
					</div>
					
					<div>
					<h1>System load (one minute average)</h1>
					<div id="load" class="plot2">
				
					</div>
					</div>
				</div>
				</td>
			</tr>
		</table>
		
		<script language="javascript">
					var loading = false;
	  				function load(arg)
					{
						if (loading)
							return;
						var params = new Array();
						if (arg == "gc")
							params.push("gc=1");
						if (arg == "fin")
							params.push("fin=1");							
						params.push("length=" + $("[name=length]").val());
						loading = true;
						$.ajax({
							url: "/systemstatus/getMonitorData.do?out=json&amp;nodb=1",
							dataType: "json",
							data: params.join("&amp;"),
							success: function(obj)
							{
								loading = false;
								var plot = new Array(new Plot(), new Plot());
								var dt;
								
								plot.minY = 100000000;
								
								var add = function(plotNum, label, color)
								{
									dt = new Array();
									for (var i = 0; i &lt; obj[plotNum].length; i++)
										dt.push([i * 2, obj[plotNum][i]]);
									plot[plotNum].addData(dt, color).setLabel(label);
									plot[plotNum].currentData.lines = {show: true};
									plot[plotNum].currentData.points.show = false;
									plot[plotNum].currentData.hoverable = true;
								}
								
								add(0, "Memory usage", "green");
								add(1, "System load", "red");
								
								for (var i = 0; i &lt; plot.length; i++)
								{
									plot[i].currentData.bars = {show: false};
									plot[i].currentData.lines.show = true;
									plot[i].render(".plot" + (i + 1));
								}
								
								
								$("#memory").bind("plothover", function (event, pos, item) {
							        if (item)
							          	tooltip.show("" + item.datapoint[1] + "Mb used ", 200);
							        else
							        	tooltip.hide();
							    });
								
								$("#memory").bind("mouseout", function (event, pos, item) {
							       tooltip.hide();
							    });
							},
							error: function(obj)
							{
								loading = false;
								//window.alert("Error");
							},
							after: function() {
								loading = false;
							}
						});
					}
					
					$(load);
					
					$("[name=length]").change(load);
					
					setInterval(load, 3000);
				</script>
	</xsl:template>
	
</xsl:stylesheet>
