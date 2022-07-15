<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />

	<xsl:template name="content">
		<style type="text/css">
			TABLE.torefactor TD, TABLE.torefactor TH {background-color: #FFD; padding: 10px; text-align: left; border-bottom: 1px solid white;}
			TABLE.torefactor TH {font-style: bold; border-bottom: 1px solid black;}
			TR.ready TD {background-color: #dafcb3;}
			SPAN.filter {padding-right: 10px;}
		</style>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<title>Published tasks</title>
		<table height="100%" width="100%">
			<tr><td class="itunes-up">
				<h1 id="page-title">Published tasks</h1>
				The overview of all published tasks
				</td></tr>
			<tr>
				<td class="itunes-right">
				<span class="filter">
				<select filter="1" name="task-type">
					<option value="">All tasks types</option>
					<option value="MODEL_TRAINING">Model training tasks</option>
					<option value="MODEL_APPLICATION">Prediction tasks</option>
					<option value="DESCRIPTOR_CALCULATION">Descriptor calculation tasks</option>
					<option value="TOXALERT_SCREENING">ToxAlert screening tasks</option>
					<option value="SET_COMPARISON">SetCompare tasks</option>
					<option value="EXPERIMENTAL_DESIGN">Experimental design tasks</option>
				</select>
				</span>
				<span class="filter">
				Partial task name: <input name="task-name" type="text" filter="1"/>
				</span>
				<span class="filter">
					Publication ID: <input name="article-id" type="text" filter="1"/> 
				</span>
				<span class="filter">
					<input type="checkbox" name="own-tasks-only" filter="1"/> My tasks only 
				</span>
				<input name="published" type="hidden" filter="1" value="1"/>
				<a href="javascript:sampleBrowser.request(true)">[Refresh]</a><br/><br/>
					<div class="pager-strip">
						<span><b class="showed">none</b> of <b class="total">none</b></span>
						<div class="pgr">
						</div>
					</div>
					
					<table class="torefactor">
						<tr>
							<th>User</th>
							<th>Task type / <br/>Time started</th>
							<th>Task name / <br/>Publication</th>
							<th>Model</th>
							<th>Property /<br/>Set</th>
							<th>Method</th>
							<th></th>
						</tr>
						<tbody id="Browser">
						</tbody>
					</table>
					
					<div class="pager-strip">
						<span><b class="showed">none</b> of <b class="total">none</b></span>
						<div class="pgr">
						</div>
					</div>
				</td>
			</tr>
		</table>
		
		<div id="xmlDialog">
			<div class="hd">XML Configuration of the model</div> 
   		 	<div class="bd">
   		 		<div style="overflow: auto; height: 300px;">
					<pre>
						<xsl:value-of select="configurationXml"/>
					</pre>
				</div>
			</div>
		</div>
		
		<script language="javascript">
			
			var xmlDialog = new YAHOO.widget.Dialog("xmlDialog", { 
			    	width:"700px", 
					fixedcenter:true, 
					modal:true, 
					visible:false 
		    });
		    
		    xmlDialog.render();
		
			var PendingTasksBrowser = function()
			{
				var self = this;
				this.controller = "pendingtasks";
				Browser.call(this);
				this.itemTemplate = "js/templates/published-pending-task.ejs"
				this.itemElement = "pending-task";
				this.useTableRows = true;
				this.pager.selectors.pager = ".pgr";
				
				this.onItemDrawn = function()
				{
					// Replace HTML tags
					if (this.currentEntity.model.configurationXml)
					{
						this.currentEntity.model.configurationXml = this.currentEntity.model.configurationXml.replace(/&gt;/g, "&amp;gt;");
						this.currentEntity.model.configurationXml = this.currentEntity.model.configurationXml.replace(/&lt;/g, "&amp;lt;");
					}
					if (this.currentEntity.status == "ready")
						this.currentBlock.addClass("ready");
						
					var sel = this.currentBlock.find('select');
					var curId =this.currentEntity.id;	
					sel.val(this.currentEntity.priority);
					sel.change(function()
					{
						self.setPositionById(curId);
						self.callAction("set_priority");	
					});
				}
				
				this.listenEvent("items_loaded", function(){
					self.filters.setValue("published", 1);
				});
				
				
				this.doShowxml = function()
				{
					$("#xmlDialog pre").html(this.currentEntity.model.configurationXml);
					xmlDialog.show();
				}
				
				this.listenEvent("items_loaded", function(){
					updateVisibility();
				});
			}
		
			include.plugins('view');
			var sampleBrowser = new PendingTasksBrowser();
			$(document).ready(function() {sampleBrowser.initialize();});
			
			var ajax = new QSPR.Ajax();
			
			
			
			function updateVisibility()
			{
			}
			
			var taskTypeMap = new Object();
			taskTypeMap["MODEL_TRAINING"] = "Model training";
			taskTypeMap["MODEL_APPLICATION"] = "Prediction";
			taskTypeMap["DESCRIPTOR_CALCULATION"] = "Calculation of descriptors";
			taskTypeMap["DM_CALCULATION"] = "Calculation of DM";
			taskTypeMap["TOXALERT_SCREENING"] = "ToxAlert screening";
			taskTypeMap["SET_COMPARISON"] = "SetCompare task";
			taskTypeMap["EXPERIMENTAL_DESIGN"] = "Experimental Design";
		</script>
	</xsl:template>
</xsl:stylesheet>