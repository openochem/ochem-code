<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />

	<xsl:template name="content">
		<style type="text/css">
			TABLE.summary TD {padding: 0px;}
			TABLE.summary TH {background-color: #DDF; font-weight: bold; padding: 0px;}
			TABLE.summary TD, TABLE.summary TH {border-bottom: 1px solid white; border-right: 1px solid white; text-align: center;}
			TABLE.summary TD DIV, TABLE.summary TH DIV {position: relative; padding: 5px 15px 5px 15px;}
			TD.status_published {background-color: #DDF;}
			TD.status_saved {background-color: #AFA;}
			TD.status_ready {background-color: #EFE;}
			TD.status_running, TD.status_init {background-color: #EEE;}
			TD.status_error {background-color: #FEE;}
			TD.status_invalid {background-color: #FDD;}
			TD.status_ A {color: #888;}
			TD.status_in {color: #666;}
			TABLE.summary TD IMG, TABLE.summary TH IMG {position: absolute; bottom: 0px; right: 0px;}
			.command-panel A {margin-right: 10px;}
			.command-panel IMG {margin-right: 3px;}
		</style>
	  	<script language="javascript" type="text/javascript" src="js/lib/excanvas.min.js"></script>
	  	<script language="javascript" type="text/javascript" src="js/lib/jquery.flot-0.6.min.js"></script>
	  	<script language="javascript" src="js/commons/plotting.js"></script>
	  	<script language="javascript" src="js/blocks/model-profile.js"></script>
		<title>Multiple Models Builder</title>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<table width="100%">
			<tr><td class="itunes-up">
				<h1>Multiple models overview</h1>
			</td></tr>
			<tr><td class="itunes-right">
				Predicted property: <a href="properties/edit.do?id={/model/others/multiple-models[1]/trainingSet/properties[1]/@id}" tab="Predicted property"><xsl:value-of select="/model/others/multiple-models[1]/trainingSet/properties[1]/@name"/></a><br/>
				Training set: <a href="basket/edit.do?id={/model/others/multiple-models[1]/trainingSet/@id}" tab="Training set"><xsl:value-of select="/model/others/multiple-models[1]/trainingSet/@name"/></a>
				<xsl:if test="/model/others/multiple-models[1]/trainingSetVersionCount &gt; 1">
					(<xsl:value-of select="/model/others/multiple-models[1]/trainingSetVersionCount"/> different versions detected)
					<a class="infolink" title="The models overviewed here are based on several different versions of this training set. That is caused by the changes in the training set, e.g. manual exclusion of some records."></a>
				</xsl:if>
				<br/><br/>
							
				Metrics <select name="metrics">
					<xsl:choose>
						<xsl:when test="/model/others/multiple-models[1]/modelType = 'regression'">
							<option value="rmse" class="regression">RMSE - Root Mean Square Error</option>
							<option value="mae" class="regression">MAE - Mean Absolute Error</option>
							<option value="r2" class="regression">R2</option>
							<option value="q2" class="regression">Q2</option>	
						</xsl:when>
						<xsl:when test="/model/others/multiple-models[1]/modelType = 'classification'">
							<option value="auc" class="classification">AUC</option>
							<option value="balancedAccuracy" class="classification">Balanced accuracy</option>
							<option value="accuracy" class="classification">Accuracy</option>
							<option value="r2" class="regression">MCC</option>
						</xsl:when>
						<xsl:when test="/model/others/multiple-models[1]/modelType = 'mixed'">
							<option value="rmse" class="regression">RMSE - Root Mean Square Error</option>
							<option value="mae" class="regression">MAE - Mean Absolute Error</option>
							<option value="r2" class="regression">R2 or MCC</option>
							<option value="q2" class="regression">Q2</option>	
							<option value="auc" class="classification">AUC</option>
							<option value="balancedAccuracy" class="classification">Balanced accuracy</option>
							<option value="accuracy" class="classification">Accuracy</option>
						</xsl:when>									
					</xsl:choose>
					<option value="size">Model size</option>
					<option value="records">Records</option>
					<option value="errors">Errors</option>
				</select>
				for
				<select name="set">
					<option value="training">Training set</option>
					<option value="validation">Validation set</option>
					<option value="excluded">Excluded set</option>
					<option value="whole">Training + Excluded</option>
				</select>
				Validation: <select name="validation">
					<xsl:for-each select="/model/others/multiple-models">
						<xsl:if test="hasAnyModels = 'true'">
							<option value="{validation}"><xsl:value-of select="validationTitle"/> (<xsl:value-of select="models-count"/> models)</option>
						</xsl:if>
					</xsl:for-each>
					<option value="all">All validation protocols</option>
				</select><br/><br/>
				<xsl:for-each select="/model/others/multiple-models">
				<xsl:variable name="table" select="."/>
				<table class="summary {validation} invisible all" cellspacing="2" cellpadding="2">
					<tr>
						<th></th>
						<xsl:for-each select="methodCodes">
							<th col="{position()}">
								<div>
									<xsl:value-of select="."/>
									<a action="groupMenu" id="method-{position()}-{$table/validation}"><img src="img/icons/corner.gif"/></a>
								</div>
							</th>
						</xsl:for-each> 
					</tr>
					
					<xsl:for-each select="descriptorCodes">
						<xsl:variable name="nRow" select="position()"/>
						<tr>
							<th>
								<div>
									<xsl:value-of select="."/>
									<a action="groupMenu" id="descriptors-{position()}-{$table/validation}"><img src="img/icons/corner.gif"/></a>
								</div>
							</th>
							<xsl:for-each select="$table/methodCodes">
							<xsl:variable name="nCol" select="position()"/>
							<xsl:variable name="mData" select="$table/modelData[$nCol]/item[$nRow]"/>
							
							<td class="status_{$mData/status}" col="{$nCol}">
								<div>
									<xsl:if test="$mData/modelId">
										<xsl:attribute name="model-id"><xsl:value-of select="$mData/modelId"/></xsl:attribute>
									</xsl:if>
									<span title="{$mData/modelName} - {$mData/detailedStatus}">
									<xsl:choose>
										<xsl:when test="$mData/stats/r2">
											<a tab="Model profile">
												<xsl:attribute name="href">
													<xsl:choose>
														<xsl:when test="$mData/pendingTaskId">pendingtasks/fetchnew.do?id=<xsl:value-of select="$mData/pendingTaskId"/></xsl:when>
														<xsl:otherwise>model/profile.do?id=<xsl:value-of select="$mData/modelId"/></xsl:otherwise>
													</xsl:choose>
												</xsl:attribute>
												
												<xsl:for-each select="$mData/stats">
													<span class="m rmse {@id}"><xsl:value-of select="rmse"/></span>
													<span class="m mae invisible {@id}"><xsl:value-of select="mae"/></span>
													<span class="m average invisible {@id}"><xsl:value-of select="average"/></span>
													<span class="m r2 invisible {@id}"><xsl:value-of select="r2"/></span>
													<span class="m mcc invisible {@id}"><xsl:value-of select="r2"/></span>
													<span class="m q2 invisible {@id}"><xsl:value-of select="q2"/></span>
													<span class="m accuracy invisible {@id}"><xsl:value-of select="accuracy"/></span>
													<span class="m balancedAccuracy invisible {@id}"><xsl:value-of select="balancedAccuracy"/></span>
													<span class="m auc invisible {@id}"><xsl:value-of select="auc"/></span>
													<span class="m size invisible {@id}"><xsl:value-of select="size"/></span>
													<span class="m records invisible {@id}"><xsl:value-of select="records"/></span>
													<span class="m errors invisible {@id}"><xsl:value-of select="errors"/></span>
												</xsl:for-each>
											</a>
										</xsl:when>
										<xsl:otherwise>
											<xsl:choose>
												<xsl:when test="$mData/status = 'ready'">
													<a href="pendingtasks/fetchnew.do?id={$mData/pendingTaskId}" tab="Model profile"><xsl:value-of select="$mData/status"/></a>
												</xsl:when>
												<xsl:otherwise>
													<xsl:value-of select="$mData/status"/>
												</xsl:otherwise>
											</xsl:choose>
											
										</xsl:otherwise>
									</xsl:choose>
									</span>
									<xsl:if test="$mData/moreModels">
										<span class="more-models">(<a action="showmodels">+<xsl:value-of select="count($mData/moreModels)"/></a> models)</span>
									</xsl:if>
									
									<xsl:if test="$mData/modelId">
										<a action="menu" id="model-{$mData/modelId}"><img src="img/icons/corner.gif"/></a>
									</xsl:if>
									<div class="invisible description">
										<xsl:value-of select="$mData/modelName"/>
										<hr/>
										<xsl:value-of select="$mData/description"/>	
										<xsl:if test="$mData/status = 'error'">
											<hr/><b>Error:</b><xsl:value-of select="$mData/detailedStatus"/>
										</xsl:if>
										<xsl:if test="$mData/status = 'running'">
											<hr/><b>The model is running: </b> <xsl:value-of select="$mData/detailedStatus"/>
										</xsl:if>
									</div>
								</div>
								<xsl:for-each select="$mData/moreModels">
									<xsl:apply-templates select="."/>
								</xsl:for-each>
								<xsl:if test="not($mData/modelId)">
									<a action="cross_over" title="Create a model with this configuration">+</a>
								</xsl:if>
							</td>
						</xsl:for-each>
						</tr>
					</xsl:for-each>
				</table>
				</xsl:for-each><br/>
				<div class="command-panel">

					<xsl:choose>
						<xsl:when test="/model/param = 'consider-predicates'">
							<a class="long" href="multiplemodels/show.do?set={/model/others/multiple-models[1]/trainingSet/@id}"><img src="img/icons/refresh.png"/>Predicates ("&lt;", "&gt;") and/or optimal classification thresholds are used for statistics (click to change)</a>         
		    			</xsl:when>
			    	 	<xsl:otherwise>
							<a class="long" href="multiplemodels/show.do?set={/model/others/multiple-models[1]/trainingSet/@id}&amp;consider-predicates=1"><img src="img/icons/refresh.png"/>Predicates ("&lt;", "&gt;") and/or optimal classification thresholds are NOT used for statistics (click to change)</a>         
						</xsl:otherwise>
	      			</xsl:choose>
      			<br/><br/>
					<a class="long" href="javascript:window.location.reload()"><img src="img/icons/refresh.png"/>Refresh</a>
					<xsl:if test="//item[status='invalid']">
						<a class="long" href="multiplemodels/show.do?set={/model/others/multiple-models[1]/trainingSet/@id}&amp;deleteInvalidModels=1">Delete <xsl:value-of select="count(//item[status='invalid'])"/> invalid model(s)</a>
					</xsl:if><br/>
					<br/><a href="multiplemodels/exportReportAsXls.do"><img src="img/icons/xls.gif"/>Export as Excel file</a><br/>
					<a href="multiplemodels/exportReportAsR.do"><img src="img/icons/R.gif"/>Export as R script</a>
				</div>
			</td></tr>
		</table>
		
		<div id="waitingDialog"> 
		    <div class="hd">Please wait</div> 
		    <div class="bd" style="text-align: center;"> 
		        <span class="message">Please wait until action is completed.<br/>
		        It may take a while.
		        </span>
		        <br/>
		        <img src="img/roller_small.gif"/> 
		    </div> 
		</div>
		
		<script language="javascript">
			var slMetrics = $("select[name=metrics]");
			var slSet = $("select[name=set]");
			var summary = new Actionable();
			summary.ajax.url = "multiplemodels/action.do";
			
			$("select[name=metrics]").change(function(){
				$(".summary span.m").addClass("invisible");
				$(".summary span." + slMetrics.val() + "." + slSet.val()).removeClass("invisible");
			});
			
			$("select[name=set]").change(function(){
				$(".summary span.m").addClass("invisible");
				$(".summary span." + slMetrics.val() + "." + slSet.val()).removeClass("invisible");
			});
			
			$("select[name=validation]").change(function(){
				$(".summary").addClass("invisible");
				$(".summary." + $(this).val()).removeClass("invisible");
			});
			
			
			$("select[name=validation]").change();
			slMetrics.change();
			
			$(document).ready(function() {
				summary.cellMenu = new YAHOO.widget.Menu("cellMenu");
				summary.cellMenu.render();
				
				summary.groupMenu = new YAHOO.widget.Menu("groupMenu");
				summary.groupMenu.render();
				
				summary.xmlDialog = new YAHOO.widget.Dialog("xmlDialog", { 
			    	width:"700px", 
					fixedcenter:true, 
					modal:true, 
					visible:false 
			    });
			    
			   summary.iterationsDialog = new YAHOO.widget.Dialog("iterationsDialog", { 
			    	width:"500px", 
					fixedcenter:true, 
					modal:true, 
					visible:false 
			    });
			    
			    summary.xmlDialog.setContent = function(content)
			    {
			    	$("#configuration-xml").html(content);
			    }
			    
			    summary.xmlDialog.render();
			    
			    summary.iterationsDialog.setContent = function(content)
			    {
			    	$("#rocarea").html(content);
			    }
			    
			    summary.iterationsDialog.render();
				
				$(".summary DIV[model-id] SPAN").each(function()
				{
					var desc = $(this).parent("div").find("div.description").html();
					if (desc)
						$(this).attr("title", (desc).replace(/\n/g, "<br/>"));
					$(this).tooltip();
				});
				
				$(".command-panel A.long").click(function(){
					summary.waitingDialog.setContent("It may take a while to generate a report for a large number of models<br/>(from seconds to several minutes).<br/>Please, wait patiently");
					summary.waitingDialog.show();
					return true;
				});
			});
			
			
			summary.doMenu = function(link)
			{
				summary.cell = $(link).parents("div[model-id]");
				summary.currentModelId = summary.cell.attr("model-id");
				var status = summary.cell.find('span').html();
				$("a[action=terminate]").setClass("invisible", !(status == "init" || status == "running"));
				$("a[action=save]").setClass("invisible", (status == "ready" || status == "saved"));
				$("a[action=export_model]").setClass("invisible", (status == "ready" || status == "saved"));
				$("a[action=pendingtasks]").setClass("invisible", !(status == "init" || status == "running"));
				summary.cellMenu.cfg.setProperty("context", [$(link).attr('id'), "tl", "tr"]);
				summary.cellMenu.show();
			}
			
			summary.doExport_configuration = function(link)
			{
				openTab("Export XML configuration", "model/exportModelXml.do?id=" + summary.currentModelId);
			}

			summary.doExport_model = function(link)
			{
				openTab("Export model", "model/exportModel.do?id=" + summary.currentModelId);
			}
						
			summary.doGroupmenu = function(link)
			{
				var parent = $(link).parent().parent();
				
				var tds = parent.attr("col") ? 
					$(link).parents("table.summary").find("td[col="+parent.attr("col")+"]") :
					parent.parent().find("td");
				summary.matchingTDs = tds;
				var matchingModelsCount = tds.find("div[model-id]").length;
				var missingModelsCount = tds.filter(".status_").length;
				
				var matchingModelIds = [];
				$(tds).find("div[model-id]").each(function(){
					matchingModelIds.push($(this).attr("model-id"));	
				});
				
				summary.currentModelId = matchingModelIds.join(",");
				
				$("a[action=create_missing_models]").setClass("invisible", missingModelsCount == 0);
				
				$("#missing-models-count").html(missingModelsCount);
				$(".matching-models-count").html(matchingModelsCount);
				
				summary.groupMenu.cfg.setProperty("context", [$(link).attr('id'), "tl", "tr"]);
				summary.groupMenu.show();
			}
			
			summary.doDelete_matching_models = function(link)
		    {
		    	if (!window.confirm("Are your sure you want to delete these models?\nThis cannot be undone"))
		    		return;
		    		
		    	summary.callAction("delete", link, {success: function(){
		    		summary.matchingTDs.find("div[model-id]").find("span").html("deleted");
		    	}});
		    }

		    summary.doCreate_missing_models_ready = function()
		    {
		   		summary.doCreate_missing_models_do("true")
	    	}
		    
		    summary.doCreate_missing_models = function(){
		    	summary.doCreate_missing_models_do(false)
		    }
		    
		    summary.doCreate_missing_models_do = function(saved)
		    {
		    	var descriptorModels = [];
		    	var methodModels = [];
		    	
		    	var filter = function(div)
				{
					return $(this).attr('model-id') > 0;
				}
		    	
		    	summary.matchingTDs.filter(".status_").each(function(){
		    		descriptorModels.push($(this).parents("tr").eq(0).find("div[model-id]").filter(filter).attr("model-id"));
		    		methodModels.push($(this).parents("table.summary").find("td[col="+$(this).attr("col")+"]").find("div[model-id]").filter(filter).eq(0).attr("model-id"));
		    	});
		    	
		    	new LongOperation({
		    		url: "multiplemodels/createMultipleCrossOverModels.do",
		    		data: "method-from=" + methodModels.join(",") + "&amp;descriptors-from=" + descriptorModels.join(",")+"&amp;saved="+saved,
		    		finished: function() {
		    			window.alert("Models have been successfully started");
		    		}
		    	}).start();
		    }
			
			summary.doCross_over = function(link)
			{
				var filter = function(div)
				{
					return $(this).attr('model-id') > 0;
				}
				
				var modelDescriptors = $(link).parents("tr").eq(0).find("div[model-id]").filter(filter).attr("model-id");
				var td = $(link).parents("td").eq(0);
				var modelMethod = $(link).parents("table.summary").find("td[col="+td.attr("col")+"]").find("div[model-id]").filter(filter).eq(0).attr("model-id");

				openTab("Create a model", "multiplemodels/createCrossOverModel.do?method-from=" + modelMethod +  "&amp;descriptors-from=" + modelDescriptors);				
			}
			
			summary.doRecalculate = function(link)
			{
				window.alert("Not implemented yet");
			}
			
			summary.doConfiguration = function(link)
			{
				summary.xmlDialog.setContent('<img src="img/roller_transparent.gif"/>');
				summary.xmlDialog.show();
				$.ajax({
					url: "multiplemodels/getModelXml.do",
					data: "out=json&amp;model=" + summary.currentModelId,
					dataType: "json",
					success: function(response){
						var xml = response.message.message;
						xml = xml.replace(/&gt;/g, "&amp;gt;");
						xml = xml.replace(/&lt;/g, "&amp;lt;");
						summary.xmlDialog.setContent(xml);	
					},
					error: function(){
						summary.xmlDialog.setContent("Failed to load configuration XML.");
					}
				});
			}
			
			summary.doIterations = function(link)
			{
				summary.iterationsDialog.setContent('<img src="img/roller_transparent.gif"/>');
				summary.iterationsDialog.show();
      			var ajax = new QSPR.Ajax("model/iterationsPlot.do?id=" + summary.currentModelId);
				ajax.send({
					data: "",
					success: function(reply)
					{
						$("#rocarea").html('<div id="roc" style="width: 450px; height: 450px;"></div>');
						var rocCurve = new Plot();
       					rocCurve.selector = "#roc";
       					rocCurve.options.legend.position = "ne";
       					rocCurve.drawLine(0, 0, 0.0001, 0.0001);
				       	rocCurve.currentData.lines = {lineWidth: 1};
				       	rocCurve.currentData.shadowSize = 0;
				       	var rocData = array(reply.others.rocCurve);
       					for (var set = 0; set &lt; rocData.length; set++)
       					{
	       					var dt = [];
	       					var rocPoints = array(rocData[set].points);
		       				for (var i = 0; i &lt; rocPoints.length; i++)
				       			dt.push([parseFloat(rocPoints[i].x), parseFloat(rocPoints[i].y), parseInt(rocPoints[i].id)]);
		       				rocCurve.addData(dt, setColors[rocData[set].setId]).setLabel(rocData[set].setId);
		       				rocCurve.currentData.setId = set;
		       				rocCurve.setPoints({line: true});
				       		rocCurve.currentData.lines = {lineWidth: 1};
			       		}
	       				rocCurve.render("#roc");
					}
				});

			}
			
			
			summary.doShowmodels = function(link)
			{
				$(link).parents("td").eq(0).find("div[model-id]").removeClass("invisible");
				$(link).parents("div[model-id]").find("SPAN.more-models").addClass("invisible");
			}
			
			summary.onTerminateSuccess = function()
			{
				summary.cell.find("span").html("terminate requested");
			}
			
			summary.onSaveSuccess = function()
			{
				summary.cell.find("span").html("saved");
			}

			summary.getActionQuery = function(action)
			{
				return "modelId=" + summary.currentModelId;
			}
			
			summary.onDeleteSuccess = function(action)
			{
				summary.cell.find("span").html("deleted");
			}
			
			summary.doCreate = function()
			{
				openTab("Create a model by template", "modelconfigurator/createFromTemplate.do?model=" + summary.currentModelId);
			}
			
			summary.doPendingtasks = function()
			{
				openTab("Pending tasks for a model", "pendingtasks/tasks.do?model=" + summary.currentModelId)
			}
			
			summary.waitingDialog = new YAHOO.widget.Dialog("waitingDialog", { 
		    	width:"325px", 
				fixedcenter:true, 
				modal:true, 
				visible:false, 
				close: false
		    });
		    
		    summary.waitingDialog.setContent = function(content)
			{
			   	$("#waitingDialog SPAN.message").html(content);
			}
		    
		    summary.waitingDialog.render();
		    
		    summary.ajax.beforeRequest = function(){
		    	summary.waitingDialog.show();
		    }
		    
		    summary.ajax.afterRequest = function(){
		    	summary.waitingDialog.hide();
		    }
		    
		    
		</script>
		
		<div id="cellMenu" class="yuimenu">
	    	<div class="bd">
		        <ul class="first-of-type">
		        	<li class="yuimenuitem no-trash">
		                <a action="pendingtasks" class="yuimenuitemlabel">
		                    Open in the pending tasks list
		                </a>
		            </li>
		            <li class="yuimenuitem no-trash">
		                <a action="create" class="yuimenuitemlabel">
		                    Create another model using this configuration
		                </a>
		            </li>
		            <li class="yuimenuitem no-trash">
		                <a action="save" class="yuimenuitemlabel">
		                    Save this model
		                </a>
		            </li>
		            <li class="yuimenuitem no-trash">
		                <a action="export_model" class="yuimenuitemlabel">
		                    Export this model
		                </a>
		            </li>
		            <li class="yuimenuitem no-trash">
		                <a action="configuration" class="yuimenuitemlabel">
		                    Show XML configuration
		                </a>
		            </li>
		            <li class="yuimenuitem no-trash">
		                <a action="iterations" class="yuimenuitemlabel">
		                    Training/Validation performance
		                </a>
		            </li>
		            <li class="yuimenuitem no-trash">
		                <a action="export_configuration" class="yuimenuitemlabel">
		                    Export XML configuration
		                </a>
		            </li>
		            <li class="yuimenuitem no-trash">
		                <a action="terminate" class="yuimenuitemlabel">
		                    Terminate the task
		                </a>
		            </li>
		            <li class="yuimenuitem no-trash">
		                <a action="delete" class="yuimenuitemlabel">
		                    Delete the model
		                </a>
		            </li>
		        </ul>
	        </div>            
    	</div>
    	
    	<div id="groupMenu" class="yuimenu">
	    	<div class="bd">
		        <ul class="first-of-type">
		        	<li class="yuimenuitem no-trash">
		                <a action="delete_matching_models" class="yuimenuitemlabel">
		                    Delete <span class="matching-models-count"></span> matching models
		                </a>
		            </li>
		            <li class="yuimenuitem no-trash">
		                <a action="create_missing_models" class="yuimenuitemlabel">
		                    Create <span id="missing-models-count"></span> missing models
		                </a>
		            </li>
		            <li class="yuimenuitem no-trash">
		                <a action="create_missing_models_ready" class="yuimenuitemlabel">
		                    Create models for those saved in the first column
		                </a>
		            </li>
		        </ul>
	        </div>            
    	</div>
    	
    	<div id="xmlDialog">
			<div class="hd">XML Configuration of the model</div> 
   		 	<div class="bd">
   		 		<div style="overflow: auto; height: 300px;">
					<pre id="configuration-xml">
					</pre>
				</div>
			</div>
		</div>
		
		<div id="iterationsDialog">
			<div class="hd">Training of the model</div> 
   		 	<div class="bd">
   		 		<div style="overflow: auto; height: 500px;" id="rocarea">
				</div>
			</div>
		</div>
		
	</xsl:template>
	
	<xsl:template match="moreModels">
		<div model-id="{modelId}" class="invisible">
			<span>
			<xsl:choose>
				<xsl:when test="stats/r2">
					<a tab="Model profile">
						<xsl:attribute name="href">
							<xsl:choose>
								<xsl:when test="pendingTaskId">pendingtasks/fetchnew.do?id=<xsl:value-of select="pendingTaskId"/></xsl:when>
								<xsl:otherwise>model/profile.do?id=<xsl:value-of select="modelId"/></xsl:otherwise>
							</xsl:choose>
						</xsl:attribute>
						
						<xsl:for-each select="stats">
						
						<span class="m rmse {@id}"><xsl:value-of select="rmse"/></span>
							<span class="m mae invisible {@id}"><xsl:value-of select="mae"/></span>
							<span class="m average invisible {@id}"><xsl:value-of select="average"/></span>
							<span class="m r2 invisible {@id}"><xsl:value-of select="r2"/></span>
							<span class="m q2 invisible {@id}"><xsl:value-of select="q2"/></span>
							<span class="m accuracy invisible {@id}"><xsl:value-of select="accuracy"/> %</span>
							<span class="m balancedAccuracy invisible {@id}"><xsl:value-of select="balancedAccuracy"/> %</span>
							<span class="m auc invisible {@id}"><xsl:value-of select="auc"/> %</span>
							<span class="m size invisible {@id}"><xsl:value-of select="size"/></span>
							<span class="m records invisible {@id}"><xsl:value-of select="records"/></span>
							<span class="m errors invisible {@id}"><xsl:value-of select="errors"/></span>
						</xsl:for-each>
					</a>
				</xsl:when>
				<xsl:otherwise>
					<xsl:choose>
						<xsl:when test="status = 'ready'">
							<a href="pendingtasks/fetchnew.do?id={pendingTaskId}" tab="Model profile"><xsl:value-of select="status"/></a>
						</xsl:when>
						<xsl:otherwise>
							<xsl:value-of select="status"/>
						</xsl:otherwise>
					</xsl:choose>
				</xsl:otherwise>
			</xsl:choose>
			<xsl:if test="moreModels">
				<span class="more-models">(+<xsl:value-of select="count(moreModels)"/> models)</span>
			</xsl:if>
			</span>
			<xsl:if test="modelId">
				<a action="menu" id="model-{modelId}"><img src="img/icons/corner.gif"/></a>
			</xsl:if>
			<div class="invisible description">
				<xsl:value-of select="modelName"/>
				<hr/>
				<xsl:value-of select="description"/>
			</div>
		</div>
	</xsl:template>
</xsl:stylesheet>