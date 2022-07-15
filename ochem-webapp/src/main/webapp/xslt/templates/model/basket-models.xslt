<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">

	<xsl:include href="../../helper.xslt" />
	
	<xsl:template name="content">
		<title>Basket models overview</title>
		<style type="text/css">
			.modelblock {background-color: #EEF; padding: 10px; margin: 5px;}
			table.confusion-matrix {margin: 10px 10px 10px 10px;}
	  		table.confusion-matrix TD, table.tiny TH {padding: 5px; border: 1px solid black;}
	  		table.confusion-matrix TD {text-align: center;}
	  		.classificationSummary TD {padding: 15px; text-align: center;}
	  		.matrix{background-color: #C1E1FB}
			.real{background-color: #BBEEBB}
			.predicted{background-color: #EEAAAA}
	  	</style>
	  	<script language="javascript" type="text/javascript" src="js/lib/excanvas.pack.js"></script>
	  	<script language="javascript" type="text/javascript" src="js/lib/jquery.flot-0.6.min.js"></script>
	  	<script language="javascript" type="text/javascript" src="js/commons/plotting.js"></script>
	  	<script language="javascript" type="text/javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
	  	<script language="javascript" type="text/javascript" src="js/commons/browser.js"></script>
	  	<script language="javascript">
	  		var modelCount = 0;
	  	</script>
	  	<table width="100%">
	  		<tr>
	  			<td class="itunes-up">
	  				<h1>Basket's models overview</h1>
	  				Review all the models. based on selected basket
	  			</td>
	  		</tr>
	  		<tr><td class="itunes-right">
	  			Filter by name: <input type="text" id="name"/><br/>
					<xsl:for-each select="//others/modelMapping">
						<div style="float: left; width: 450px; overflow:auto;" id="model{position()}" class="modelblock" mapping_id="{@id}" model_name="{model/@name}">
							<span><a tab="Model profile" href="model/profile.do?mapping_id={@id}"><xsl:value-of select="model/@name"/></a></span><br/>
							<xsl:value-of select="property/@name"/> by <xsl:value-of select="model/template/@name"/><br/> 
							<span class="info">&#160;</span>
							<div id="placeholder{position()}" style="width:430px; height:200px;"></div>
							<script language="javascript">
								modelCount++;
								$(document).bind("reloadCharts", null, function(){
									loadModel(<xsl:value-of select="position()"/>, <xsl:value-of select="@id"/>, <xsl:value-of select="model/@id"/>);
									adOverview.addModel(<xsl:value-of select="@id"/>, "<xsl:value-of select="model/@name"/>");
								});
							</script>
						</div>
					</xsl:for-each>
				</td>
			</tr>
		</table>
		
		<table><tr>
			<td>
			Cumulative accuracy chart:
			<div id="ad-overview" style="width: 700px; height: 500px;"></div>
			</td>
			<td>
			<span id="percentage">20</span>% (<a href="javascript:adOverview.changePercentage(10);">+</a> / <a href="javascript:adOverview.changePercentage(-10);">-</a>) best vs. overall accuracy 
			<div id="ad-overview-cross" style="width: 500px; height: 500px;"></div>
			</td>
			<td id="legend"></td>
			</tr>
		</table>
		
		<script language="javascript">
		
		var plotSet = new PlotSet();
		include.plugins('view');
		var view = 0;
		var filterPoints = false;
		var filterCriteria = "";
		var filters = new Filters();
		filters.useUrlParameters = false;
		
		if (getParams["articles"])
			filters.setValue("articles", true);
		
		function loadModel(num, mappingId, modelId)
		{
			var ajax = new QSPR.Ajax();
			ajax.url = 'model/loadstatistics.do';
			var placeholder = $("#placeholder" + num);
			placeholder.html("&lt;img src='img/roller.gif'/&gt;");
			
			ajax.send({
				data: 'mm_id=' + mappingId + "&amp;" + filters.getQueryString(),
				success: function(json)
				{
					errors = 0;
					var sets = array(json.statistics.set);
					if (view == 0)
							view = new View({url: "js/templates/models/model-stats.ejs"});
						$("#model" + num+" span.info").html(view.render(json.statistics));
					$("#model" + num+" a").attr("title", json.others.modelMapping.model.description.replace(/\n/g, "<br/>"));
					
					if (sets[0].classificationSummary)
					{
						// function drawConfusionMatrix(classificationSummary, placeholder) {}
						var mainTable = $('<table class="classificationSummary"></table>');
						var maintr;
						mainTable.append(maintr = $("<tr></tr>"));
						placeholder.html("");
						for (var set = 0; set &lt; sets.length; set++)
						{
							// This is a classification model. Draw a confision matrix
							var maintd;
							var table ;
							var tr;
							var td;
							
							var nodes = array(sets[set].classificationSummary.nodes);
							var nodes = array(sets[set].classificationSummary.nodes);
							if(nodes.length &lt; 1)
								continue;
							
							var options = array(json.others.option);
							var maxOption = 0;
							for (var i = 0; i &lt; nodes.length; i++)
								if (nodes[i].real &gt; maxOption)
									maxOption = parseInt(nodes[i].real);
							
							//one additional for table header		
							maxOption = maxOption + 1;
							maintr.append(maintd = $("<td></td>"));
							maintd.append(table = $('<table class="confusion-matrix"></table>'));
							for (var i = 0; i &lt;= maxOption; i++)
							{
								table.append(tr = $("<tr></tr>"));
								for (var k = 0; k &lt;= maxOption; k++)
										tr.append(td = $('<td>0</td>'));
							}
							placeholder.append(mainTable);
							//set table row header
							table.find("tr").eq(0).find("td").eq(0).html("Real&#8595;/Predicted&#8594;");
							for (var i = 1; i &lt;= maxOption; i++)
							{
								table.find("tr").eq(0).find("td").eq(i).html(options[(i*1)-1].name).addClass('Predicted');
								table.find("tr").eq(i).find("td").eq(0).html(options[(i*1)-1].name).addClass('real');
								table.find("tr").eq(i).find("td").eq(i).addClass('matrix');
							}
							
							for (var i = 0; i &lt; nodes.length; i++)
								table.find("tr").eq(parseInt(nodes[i].real)+1).find("td").eq(parseInt(nodes[i].predicted)+1).html(nodes[i].count);
						
						}
					}
					else
					{
						// This is a regression model. Draw a real-predicted chart
						 
						var plot = new Plot();
						
						for (var set = 0; set &lt; sets.length; set++)
						{
							//var dt = new Array();
							if(sets[set].point)
							for (var i = 0; i &lt; sets[set].point.length; i++)
								if (!sets[set].point[i].error &amp;&amp; (!filterPoints || sets[set].point[i].selected == "true"))
								{
									var series = sets[set].point[i].articleId ? sets[set].point[i].articleId : set;
									plot.addPoint(series, [1.0*sets[set].point[i].real, 1.0*sets[set].point[i].predicted, 
											{absolutePointNum: i, recordId: 1.0*sets[set].point[i].id, selected: sets[set].point[i].selected, set: set}]);
								}
								else
									errors++;
								
							//plot.addData(dt);
						}
						
						plot.addBaseLine().render("#placeholder"+num);
						plot.pointClicked = function(series, point, data)
			       		{
			       			var absPointNum = data[2].absolutePointNum;
			       			var setNum = data[2].set;
			       			openTab("Model's compound", webRoot+"modeldot/show.do?render-mode=popup&amp;id="+modelId+"&amp;mm_id="+mappingId+"&amp;statnum=1&amp;ep_id=" + data[2].recordId);
			       		}
			       		
			       		plotSet.plots.push(plot);
			       		
			       		if (plotSet.plots.length == modelCount)
			       		{
			       			$(".modelblock a[title]").tooltip({showURL: false});
			       			//plotSet.alignRanges();
			       		}
		       		}
				},
				error: function(msg)
				{
					$("#placeholder" + num).html("Error");
				}
			});
		}
		
		
		
		// AD overview charts
		var ADOverview = function(placeholder, crossChartPlaceholder, legendPlaceholder)
		{
			var self = this;
			this.percentage = 20;
			this.plot = new Plot();
			this.plot.options.legend.position = "se";
			this.crossPlot = new Plot();
			this.crossPlot.options.legend.show = false;
			this.cnt = 0;
			this.placeholder = placeholder; 
			this.crossChartPlaceholder = crossChartPlaceholder;
			this.legendPlaceholder = legendPlaceholder;
			this.models = [];
			
			
			this.highlightModel = function(series)
			{
				for (var s in this.plot.dataArray)
					if (this.plot.dataArray[s].lines)
					if (this.plot.dataArray[s].num != series.num)
					{
						this.plot.dataArray[s].lines.lineWidth = 1;
						this.plot.dataArray[s].highlight = false;
					}
					else
					{	
						this.plot.dataArray[s].lines.lineWidth = 5;
						this.plot.dataArray[s].highlight = true;
					}
				series.lines.lineWidth = 5;
				series.highlight = true;
				self.plot.render(self.placeholder);
			}
			
			this.plot.onLegendClick = function(series)
			{
			}
			
			this.addModel = function(modelId, modelName)
			{
				var model;
				this.models = model = {id: modelId, name: modelName, accuracy: 0, accuracy: 0}; 
				this.cnt++;	
				var ajax = new QSPR.Ajax();
				ajax.url = 'model/loadAd.do';
				ajax.send({
					data: "mapping_id=" + modelId,
					
					success: function(obj)
					{
						var a;
						var dt = [];
						var adConf = obj["ad-configuration"];
						for (var i = 0; i &lt; adConf.percents.length; i++)
						{
							dt.push(a = [1.0*adConf.percents[i], 1.0*adConf.error[i]]);
							if (a[0] &lt; 20) 
								model.accuracy20 = a[1];
							model.accuracy = a[1];
						}
						self.plot.addData(dt).setLabel(modelName);
						self.plot.setPoints({line: true});
			       		self.plot.currentData.lines = {lineWidth: 1};
			       		self.plot.currentData.modelId = modelId;
			       		
			       		self.crossPlot.addData([[model.accuracy, model.accuracy20]]);
			       		self.crossPlot.currentData.original = self.plot.currentData;
			       		self.crossPlot.currentData.label = self.plot.currentData.label;
			       		self.crossPlot.currentData.hoverable = true;
					},
					
					
					error: function()
					{
						// Ignore silently
					},
					
					after: function()
					{
						self.cnt--;
						
						// The last model has arrived. Draw the chart.
						if (self.cnt == 0)
						{
							self.crossPlot.tooltips = true;
							self.plot.dataArray.sort(function(a, b){return a.label > b.label ? 1: -1});
							self.crossPlot.dataArray.sort(function(a, b){return a.original.label > b.original.label ? 1: -1});
							self.plot.render(self.placeholder);
							self.crossPlot.render(self.crossChartPlaceholder);
							
							$(self.crossChartPlaceholder).bind("plotclick", function (event, pos, item) {
						        if (item)
						          self.highlightModel(item.series.original);	
   							 });
   							 
   							 var legend = $(self.legendPlaceholder);
   							 for (var i = 0; i &lt; self.plot.dataArray.length; i++)
   							 {
   							 	var input = $('<input type="checkbox" checked="true"/>');
   							 	input.attr("series", self.plot.dataArray[i].num);
   							 	var div = $("<div>" + self.plot.dataArray[i].label + "</div>");
   							 	div.prepend(input);
   							 	legend.append(div);
   							 }
   							 
   							 legend.find("input").click(function(){
   							 	var series = self.plot.getSeriesByNum($(this).attr("series"));
   							 	series.lines.show = $(this).is(":checked");
   							 	series.label = $(this).is(":checked") ? series.title : undefined;
   							 	$("div[mapping_id=" + series.modelId + "]").setClass("invisible", !$(this).is(":checked"));
   							 	self.plot.render();
   							 });
						}
					}
				});	
			}
			
			this.refreshCrossChart = function(percentage)
			{
				var dt;
				this.crossPlot.dataArray = [];
				for (var i in this.plot.dataArray)
				{
					var series = this.plot.dataArray[i];
					if (!series.data)
						continue;
					var overallAccuracy = series.data[series.data.length - 1][1];
					var accuracy = 0;
					for (var k in series.data)
					{
						if (series.data[k][0] &lt; percentage)
							accuracy = series.data[k][1];	
					}
					
					this.crossPlot.addData([[overallAccuracy, accuracy]]);
					this.crossPlot.currentData.color = series.color;
					this.crossPlot.currentData.original = series;
					this.crossPlot.currentData.label = series.label;
				}	

								
				self.crossPlot.render(self.crossChartPlaceholder);
				$(self.crossChartPlaceholder).bind("plotclick", function (event, pos, item) {
		        	if (item)
		          		self.highlightModel(item.series.original);	
					});
			}
			
			this.changePercentage = function(inc)
			{
				this.percentage += inc;
				if (this.percentage &lt; 10)
					this.percentage = 10;
				if (this.percentage &gt; 90)
					this.percentage = 90;
					
				$("span#percentage").html("" + this.percentage);
				this.refreshCrossChart(this.percentage);
			}
		}
		
		var adOverview = new ADOverview("#ad-overview", "#ad-overview-cross", "#legend");
		
		
		var ArticleBrowser = function()
		{
			this.controller = "article";
			Browser.call(this);
			this.filters.setValue("basket", getParams["basket"]);
			this.scope = this.container = "article-browser";
			this.itemTemplate = "js/templates/article-inline.ejs";
			this.itemElement = "article";
			this.filters.setValue("pagesize", 100);
			
			this.doChoose = function()
			{
				$(document).trigger("reloadCharts");
			}
		}
		
		//var articleBrowser = new ArticleBrowser();
		
		$(document).ready(function(){
			$(document).trigger("reloadCharts");
			
			$("input#name").keydown(function(){
				$("div[model_name]").addClass("invisible");
				$("div[model_name*=" + $(this).val() + "]").removeClass("invisible");
			});
		});
	</script>
	</xsl:template>
</xsl:stylesheet>