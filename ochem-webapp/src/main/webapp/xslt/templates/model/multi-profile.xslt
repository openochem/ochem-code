<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE xsl:stylesheet [
  <!ENTITY return "&#13;">
]>

<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt"/>
	<xsl:include href="profile-common.xslt"/>
	<xsl:template name="content">
		
		<style type="text/css">
	  		table.tiny {margin: 10px 10px 10px 10px;}
	  		table.tiny TD, table.tiny TH {padding: 5px; border: 1px solid black;}
	  		table.tiny TD {text-align: center;}
	  		.classificationSummary TD {padding: 15px; text-align: center;}
	  		.ad B {display: block; margin-top: 10px;}
			.modelblock {background-color: #EEF; padding: 10px; margin: 5px;}
			.confusion-matrix{background-color: #C1E1FB}
			.real{background-color: #BBEEBB}
			.predicted{background-color: #EEAAAA}
			.model-table td
			{
				border: 1px solid #AAA;
				background-color: transparent; 
				padding: 2px;
				vertical-align: middle;
				font-size: 80%; 
				padding: 1px 5px 1px 5px;
				text-align: right;
			}
			
			.model-table
			{
				border: 1px solid #AAA; 
				background-color: transparent; 
				padding: 2px;
				align: left;
				overflow: hidden; 
				font-size: 110%; 
			}
	  	</style>
	  	<script language="javascript" type="text/javascript" src="js/lib/excanvas.pack.js"></script>
	  	<script language="javascript" type="text/javascript" src="js/lib/jquery.flot-0.6.min.js"></script>
	  	<script language="javascript" src="js/commons/plotting.js"></script>
	  	<title>Model profile</title>
		
			<table width="100%">
				<xsl:if test="not(model/taskId)">
					<tr><td class="itunes-up"><h1>Overview of a multi-model</h1>
					This is a multi-model that predicts several properties simultaneously
					</td></tr>
				</xsl:if>
				<tr><td class="itunes-right">
							<!--MODEL NAME-->
							Model name: <xsl:value-of select="model/@name"/> <a href="modelapplier/apply.do?model={model/@id}" tab="Apply a model">[apply to new compounds]</a>
							<br/>
							
							<!--Training method-->
							Training method: <xsl:value-of select="model/template/@name"/><br/>
							
							<br/>
							
							<!--Loop over properties-->
							<!--<script language="javascript">	var modelCount = 0;	</script>-->
							<xsl:for-each select="model/modelMappings">
								<div  style="float: left; width: 520px; height: 350px; overflow:auto;" id="model{position()}" class="modelblock">
								<!--<span><a tab="Model profile" href="model/profile.do?id={@id}"><xsl:value-of select="@name"/></a></span><br/>-->
								Property: <xsl:value-of select="property/@name"/> <!-- docs_info --> <!-- <a tab="Property description" href="wikipage/action.do?entities=property&amp;id={property/@id}"><xsl:value-of select="property/@name"/></a> -->
								
								<xsl:if test="unit/@id != '9'">
									measured in <xsl:value-of select="unit/@name"/>
								</xsl:if>
								(<a tab="{property/@name}" href="model/profile.do?mapping_id={@id}">Details..</a>)
								<br/>
									<span class="info">&#160;</span>
									<div id="placeholder{position()}" style="width:400px;height:190px;"></div>
									<script language="javascript">
										$(document).ready(function(){
											loadModel(<xsl:value-of select="position()"/>, <xsl:value-of select="@id"/>);
										});
									</script>
									</div>
							</xsl:for-each>
							
					<br style="clear: both;"/>		
					<div class="warning invisible">
					</div>
						
		<br style="clear: both;"/>		
		<xsl:if test="not(model/taskId)"><xsl:call-template name="recalculate"/></xsl:if>
		
		</td></tr>
		</table>
		
				
		<script language="javascript">
		
		var plotSet = new PlotSet();
		include.plugins('view');
		var view = 0;
		
		function warn(msg)
		{
			$(".warning").html(error_msg);
			$(".warning").removeClass("invisible");
		}
		
		function loadModel(num, mm_Id)
		{
			var modelId = <xsl:value-of select="model/@id"/>;
			var ajax = new QSPR.Ajax();
			ajax.url = 'model/loadstatistics.do';
			var placeholder = $("#placeholder" + num);
			placeholder.html("&lt;img src='img/roller.gif'/&gt;");
			ajax.send({
				data: 'mm_id='+mm_Id,
				success: function(json)
				{
					
					errors = 0;
					var sets = array(json.statistics.set);
					var _modelMapping = json.others.modelMapping;
					
					//model-table
					if (view == 0)
							view = new View({url: "js/templates/models/model-overview.ejs"});
						$("#model" + num+" span.info").html(view.render(json.statistics));
					$("#model" + num+" a").attr("title", json.others.modelMapping.model.description.replace(/\n/g, "<br/>"));
					if (sets[0].classificationSummary)
					{
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
							
							if(sets[set].setId == "validation")
							{
								var string = '<a tab="Test set" href="epbrowser/show.do?basket-select={model/validation-set/@id}&amp;property='+_modelMapping.property.id+'">
									<xsl:value-of select="model/validation-set/@name"/></a>   ('+_modelMapping.validationSize+')	';
								$("#model"+num).find("#validation").html(string);
							}else if(sets[set].setId == "excluded")
							{
								var string = "Excluded from training set";
								$("#model"+num).find("#excluded").html(string);
							}else
							{
								var string = '<a tab="Training set" href="epbrowser/show.do?basket-select={model/training-set/@id}&amp;property='+_modelMapping.property.id+'">';
								string = string + '<xsl:value-of select="model/training-set/@name"/></a>   ('+_modelMapping.trainingSize+')	';
								$("#model"+num).find("#training").html(string);
							}
							
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
							maintd.append(table = $('<table class="tiny"></table>'));
							for (var i = 0; i &lt;= maxOption; i++)
							{
								table.append(tr = $("<tr></tr>"));
								for (var k = 0; k &lt;= maxOption; k++)
									tr.append('<td>0</td>');
								
							}
							placeholder.append(mainTable);
							//set table column header name
							table.find("tr").eq(0).find("td").eq(0).html("Real&#8595;/Predicted&#8594;");
							for (var i = 1; i &lt;= maxOption; i++)
							{
								table.find("tr").eq(0).find("td").eq(i).html(options[(i*1)-1].name).addClass('Predicted');
								table.find("tr").eq(i).find("td").eq(0).html(options[(i*1)-1].name).addClass('real');
								table.find("tr").eq(i).find("td").eq(i).addClass('confusion-matrix');
							}
							
							for (var i = 0; i &lt; nodes.length; i++)
								table.find("tr").eq(parseInt(nodes[i].real)+1).find("td").eq(parseInt(nodes[i].predicted)+1).html(nodes[i].count);
						
						}
					}
					else
					{
						var plot = new Plot();
						plot.options.grid.autoHighlight = false;
						for (var set = 0; set &lt; sets.length; set++)
						{
							if(sets[set].setId == "validation"){
								var string = '<a tab="Test set" href="epbrowser/show.do?basket-select={model/validation-set/@id}&amp;property='+_modelMapping.property.id+'">
									<xsl:value-of select="model/validation-set/@name"/></a>('+_modelMapping.validationSize+')';
								$("#model"+num).find("#validation").html(string);
							}else if(sets[set].setId == "excluded"){
								var string = "Excluded from training set";
								$("#model"+num).find("#excluded").html(string);
							}else{
								var string = '<a tab="Training set" href="epbrowser/show.do?basket-select={model/training-set/@id}&amp;property='+_modelMapping.property.id+'">';
								string = string + '<xsl:value-of select="model/training-set/@name"/></a>('+_modelMapping.trainingSize+')	';
								$("#model"+num).find("#training").html(string);
							}
							
							var dt = new Array();
							if(sets[set].point)
								for (var i = 0; i &lt; sets[set].point.length; i++)
									if (!sets[set].point[i].error)
										dt.push([1.0*sets[set].point[i].real, 1.0*sets[set].point[i].predicted, i, 1.0*sets[set].point[i].id, 1.0*sets[set].point[i].moleculeId]);
									else
										errors++;
								
							plot.addData(dt, setColors[sets[set].setId]);
						}
						
						plot.addBaseLine().render("#placeholder"+num);
						plot.pointClicked = function(setNum, pointNum, data)
			       		{
			       			var moleculeId = data[4];
			       			for (var i in plotSet.plots)
			       			{
			       				var otherPlot = plotSet.plots[i];
			       				if (otherPlot.selector != undefined)
			       				{
			       					otherPlot.flot.unhighlight();
			       					for (var set = 0; set &lt; otherPlot.dataArray.length; set++)
			       						for (var p = 0; p &lt; otherPlot.dataArray[set].data.length; p++)
			       						{
			       							if (otherPlot.dataArray[set].data[p][4] == moleculeId)
			       								otherPlot.flot.highlight(set, p);
			       						}
			       				}
			       			}
			       			
			       			var epId = this.dataArray[setNum].data[pointNum][3];
			       			openTab("Model's compound", webRoot+"modeldot/show.do?render-mode=popup&amp;showall=true&amp;id="+modelId+"&amp;mm_id="+mm_Id+"&amp;statnum=1&amp;ep_id=" + epId);
			       		}
			       		
			       		plotSet.plots.push(plot);
		       		}
		       		//var error_msg = 'Number of compounds ignored because of <a href="javascript:errorsShow()">errors</a> in model = '+errors+'.';
		       		if(errors > 0)
		       			warn(error_msg);
		       			
		       		$(document).trigger("DOM_updated", $(document));
				},
				error: function(msg)
				{
					$("#placeholder" + num).html("Error");
				}
			});
		}
	</script>
	</xsl:template>
</xsl:stylesheet>