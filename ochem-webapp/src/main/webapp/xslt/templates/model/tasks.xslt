<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"  version="1.0">
	<xsl:include href="../../helper.xslt" />

	<xsl:template name="content">
		<style type="text/css">
			TABLE.torefactor TD, TABLE.torefactor TH {background-color: #FFD; padding: 10px; text-align: left;}
			TABLE.torefactor TH {font-style: bold; border-bottom: 1px solid black;}
			TR.ready TD {background-color: #c0e793;}
		</style>
		<script language="javascript">
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
			
			setInterval('window.location.reload()', 1000*60); // Refresh the page every minute
		</script>
		<title>Pending tasks</title>
		
		<table width="100%">
			<tr><td class="itunes-up">
						<img src="img/icons/pendingtask.png"/>
				<h1>Tasks</h1>
				Here goes the list of pending model creation tasks. The most fresh 10 tasks are shown.
			</td></tr>
			<tr><td class="itunes-right">
			<a href="javascript:location.reload(true);">[Refresh]</a>
				<table class="torefactor" cellspacing="10">
				<tr>
					<xsl:if test="session/user/rank = 10">
						<th>User</th>
					</xsl:if>
					<th>Started</th>
					<th>Model name</th>
					<th>Property</th>
					<th>Dataset</th>
					<th>Method</th>
					<th>Status</th>
					<th>Details</th>
					<th>...</th>
					<th>...</th>	
				</tr>
				<xsl:for-each select="others/model">
					<tr model="{@id}">
						<xsl:if test="status = 'ready'">
							<xsl:attribute name="class">ready</xsl:attribute>
						</xsl:if>
						<xsl:if test="/model/session/user/rank = 10">
							<td><xsl:value-of select="session/user/@login"/></td>
						</xsl:if>
						<td><i><xsl:value-of select="@datetime"/></i></td>
						<td>
							<xsl:choose>
								<xsl:when test="recalculation = 'true'">
									<a tab="Model profile" href="model/profile.do?id={@id}"><xsl:value-of select="@name"/></a>
								</xsl:when>
								<xsl:otherwise>
									<xsl:value-of select="@name"/>
								</xsl:otherwise>
							</xsl:choose>
							
						</td>
						<td><xsl:value-of select="modelMappings/property/@name"/></td>
						<td><a tab="Training set" href="basket/edit.do?id={training-set/@id}"><xsl:value-of select="training-set/@name"/></a></td>
						<td><xsl:value-of select="template/@name"/></td>
						<td>
							<xsl:value-of select="status"/>
							<xsl:if test="status = 'init' or status='assigned'">
								<img src="img/roller_small.gif"/>
							</xsl:if>
						</td>
						<td>
							<span title="{detailedStatus}" onclick="detailedStatus(this);"><xsl:value-of select="substring(detailedStatus, 1, 30)"/></span>
							<xsl:if test="string-length(detailedStatus) &gt; 30">
								<a href="#" onclick="detailedStatus(this); return false;">[more&gt;&gt;]</a>
							</xsl:if>
						</td>
						
						<td>
							<xsl:if test="status = 'init' or status='assigned'">
								<a href="javascript:void(0)" onclick="updateModel(this, 'kill'); return false;" id="{@id}">terminate</a>
							</xsl:if>
							<xsl:if test="status = 'error' or status='ready' or status='kill' or status='killed'">
								<a href="javascript:void(0)" onclick="updateModel(this, 'recalculate'); return false;" id="{@id}">recalculate</a>
							</xsl:if>
						</td>
						<td>
							<a href="#" onclick="showXml(this); return false;"><img src="img/icons/xml.jpg"/></a>
							<a href="pendingtasks/tasks.do?delete={@id}" onclick="return confirmDelete();"><img src="img/icons/delete.gif"/></a>
							<xsl:if test="status = 'ready'">
								<a tab="Review a ready model" href="pendingtasks/fetch.do?id={@id}"><img src="img/icons/save.gif"/></a>
							</xsl:if>
							<div class="invisible">
								<xsl:value-of select="model/description"/>
							</div>
						</td>
					</tr>
				</xsl:for-each>
				<xsl:if test="not(others/model)">
					<tr>
						<td colspan="5">
							<i>No pending tasks</i>
						</td>
					</tr>
				</xsl:if>	
				</table>
			</td></tr>
		</table>
		
		<div id="xmlDialog">
			<div class="hd">XML Configuration of the model</div> 
   		 	<div class="bd">
   		 		<xsl:for-each select="others/model">
	   		 		<div style="overflow: auto; height: 300px;" model="{@id}">
						<pre>
							<xsl:value-of select="configurationXml"/>
						</pre>
					</div>
				</xsl:for-each>
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
		    
		    showXml = function(link)
			{
				var id = $(link).parents("tr[model]").attr("model");
				$("#xmlDialog [model]").addClass("invisible");
				$("#xmlDialog [model=" + id + "]").removeClass("invisible");
				xmlDialog.show();
			}
		</script>
		
		
		
	</xsl:template>
</xsl:stylesheet>