<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />

	<xsl:template name="content">
		<style type="text/css">
			TABLE.torefactor TD, TABLE.torefactor TH {background-color: #FFD; padding: 10px; text-align: left; border-bottom: 1px solid white;}
			TABLE.torefactor TH {font-style: bold; border-bottom: 1px solid black;}
			TR.status-ready TD {background-color: #c0e793;}
			TR.status-error TD {color: #400;}
		</style>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<title>Pending tasks</title>
		<table height="100%" width="100%">
			<tr><td class="itunes-up">
				<h1 id="page-title">Pending tasks<a class="infolink" target="_blank" href="https://docs.ochem.eu/display/MAN/Pending+tasks+browser"></a></h1>
				The overview of all running tasks and all completed tasks awaiting your action
				</td></tr>
			<tr>
				<td class="itunes-right">
				<select filter="1" name="task-type">
					<option value="">All tasks types</option>
					<option value="MODEL_TRAINING">Model training tasks</option>
					<option value="MODEL_APPLICATION">Prediction tasks</option>
					<option value="DESCRIPTOR_CALCULATION">Descriptor calculation tasks</option>
					<option value="TOXALERT_SCREENING">ToxAlert screening tasks</option>
					<option value="SET_COMPARISON">SetCompare tasks</option>
					<option value="EXPERIMENTAL_DESIGN">Experimental design tasks</option>
				</select>
				<select filter="1" name="status">
					<option value="">All tasks statuses</option>
					<option value="init">Waiting in queue</option>
					<option value="assigned,stop">Running tasks</option>
					<option value="error,kill,killed">
						<xsl:if test="//pendingTaskFilter/status = 'error,kill,killed'">
							<xsl:attribute name="selected">true</xsl:attribute>
						</xsl:if>
						Failed tasks
					</option>
					<option value="ready">Ready tasks</option>
				</select>
				<input name="task-name" type="text" filter="1"/>
				<a href="javascript:sampleBrowser.request(true)">[Refresh]</a>
				<a action="delete_all">[Delete all matching tasks]</a>
				<a action="recalculate_errors" id="recalculate-errors" class="invisible">[Recalculate all failed tasks]</a>
				<a action="fetch_ready" id="fetch-ready" class="invisible">[Pre-fetch all ready tasks]</a>
				<xsl:if test="//user/superuser = 'true'">
					<input type="checkbox" name="other-users" filter="1">
						<xsl:if test="//pendingTaskFilter/unpublishedTasksOfOtherUsers = 'true'">
							<xsl:attribute name="checked">checked</xsl:attribute>
						</xsl:if>
					</input> See tasks of other users
				</xsl:if>
				<xsl:if test="/model/session/user/group">
					<input type="checkbox" name="group" filter="1"/><xsl:value-of select="/model/session/user/group/name"/> members' models
				</xsl:if>
					
					<table class="torefactor">
						<tr>
							<th class="user invisible">User</th>
							<th>Task type / <br/>Time started</th>
							<th>Model / Task name</th>
							<th>Property /<br/>Set</th>
							<th>Method</th>
							<th>Status</th>
							<th>Priority</th>
							<th>Details</th>
							<th></th>
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
				this.itemTemplate = "js/templates/pending-task.ejs";
				if (getParams["published"])
					this.itemTemplate = "js/templates/pending-task-published.ejs"
				this.itemElement = "pending-task";
				this.useTableRows = true;
				this.pager.selectors.pager = ".pgr";
				
				this.parentGetActionQuery = this.getActionQuery;
				this.filters.setFromUrl(); 
				
				this.getActionQuery = function()
				{
					if (getParams["debug"] != undefined)
						return this.parentGetActionQuery()+"&amp;debug=true";
					else
						return this.parentGetActionQuery();
				}
				
				this.onItemDrawn = function()
				{
					// Replace HTML tags
					if (this.currentEntity.model.configurationXml)
					{
						this.currentEntity.model.configurationXml = this.currentEntity.model.configurationXml.replace(/&gt;/g, "&amp;gt;");
						this.currentEntity.model.configurationXml = this.currentEntity.model.configurationXml.replace(/&lt;/g, "&amp;lt;");
					}
					this.currentBlock.addClass("status-" + this.currentEntity.status);
						
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
					if (maximumTaskPriority &lt; 10)
						self.mainDiv.find('option[value=10]').remove();
						
					if (isSuperUser || getParams["published"])
						$(".user").removeClass("invisible");
						
					if (getParams["published"])
					self.filters.setValue("published", 1);
					
				});
				
				this.onSet_prioritySuccess = function()
				{
					this.request(false);
				}
				
				this.beforeDelete = function()
				{
					return window.confirm("Do you really want to delete this model?");
				}
				
				this.beforeRecalculate = function()
				{
					return window.confirm("Do you really want to recalculate this task? This will start the calculation process from the very beginning");
				}
				
				this.beforeKill = function()
				{
					return window.confirm("Do you really want to terminate this task?");
				}
				
				this.onDeleteSuccess = function()
				{
					this.deleteRecord();
				}
				
				this.onRecalculateSuccess = this.onKillSuccess = this.onPublishSuccess = this.onUnpublishSuccess = function()
				{
					this.request();
				}
				
				this.doShowxml = function()
				{
					$("#xmlDialog pre").html(this.currentEntity.model.configurationXml);
					xmlDialog.show();
				}
				
				this.beforeRecalculate_errors = function()
				{
					return window.confirm("Are you sure you want to recalculate the matching models? This can be a lot of calculations.");
				}
				
				this.doRecalculate_errors = function()
				{
					new LongOperation({
			    		url: "pendingtasks/recalculateAllFailed.do",
			    		data: "out=json&amp;" + self.getActionQuery(),
			    		finished: function() {
			    			window.alert("Models have been resubmitted for recalculation.");	
			    		}
			    	}).start();
				}
				
				this.beforeDelete_all = function()
				{
					return window.confirm("Are you sure you want to delete all the matching pending tasks? This cannot be undone!");
				}
				
				this.doDelete_all = function()
				{
					new LongOperation({
			    		url: "pendingtasks/deleteAll.do",
			    		data: "out=json&amp;" + self.getActionQuery(),
			    		finished: function() {
			    			window.alert("Tasks have been deleted.");
			    			window.location.reload();	
			    		}
			    	}).start();
				}
				
				this.doFetch_ready = function()
				{
					new LongOperation({
			    		url: "pendingtasks/fetchReady.do",
			    		data: "out=json&amp;" + self.getActionQuery(),
			    		finished: function() {
			    			window.alert("Ready models have been fetched.");	
			    		}
			    	}).start();
				}
				
				this.listenEvent("items_loaded", function(){
					updateVisibility();
				});
			}
		
			include.plugins('view');
			var sampleBrowser = new PendingTasksBrowser();
			$(document).ready(function() {sampleBrowser.initialize();});
			
			
			
			var ajax = new QSPR.Ajax();
			
			function updateModel(link, type)
			{
				if (!window.confirm('Do you really want to ' + type + ' this model?'))
					return false;
				var id = $(link).attr('id');
				ajax.send(
				{ 
					url: "pendingtasks/"+type+".do?id=" + id, 
					success: function(){
						window.alert("Your request was successfull. \nRefresh this page to see actual current status");
						$(link).addClass("invisible");
					}
				});
			}
			
			function updateVisibility()
			{
				$("#recalculate-errors").setClass("invisible", $("[name=status]").val().indexOf("error") == -1 || sampleBrowser.pager.totalNum == 0);
				$("#fetch-ready").setClass("invisible", $("[name=status]").val().indexOf("ready") == -1 || sampleBrowser.pager.totalNum == 0);
			}
			
			function detailedStatus(link)
			{
				var link = $(link).parent("td").find("span");
				var more = $(link).parent("td").find("a");
				
				var expanded = more.get(0).expanded;
				if (!expanded)
					expanded = false;
				more.get(0).expanded = expanded = !expanded;
				
				if (expanded)
					more.html("[less&lt;&lt;]");
				else
					more.html("[more&gt;&gt;]");
				
				var html = link.html();
				var title = link.attr('title');
				link.html(title);
				link.attr('title', html);
			}
			
			function confirmDelete()
			{
				return window.confirm("Do you really want to delete this model?");
			}
			
			function autoRefresh()
			{
				if ($("input[name=autorefresh]:checked").length > 0)
					sampleBrowser.request(false);
			}
			
			setInterval('autoRefresh()', 1000*60); // Refresh the page every minute
			
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