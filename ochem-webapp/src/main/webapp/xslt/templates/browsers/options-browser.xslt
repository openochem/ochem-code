<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />

	<xsl:template name="content">
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<link rel="stylesheet" type="text/css" href="css/smoothness/jquery-ui-1.9.1.custom.css" />
		<style type="text/css">
			.ui-widget {font-size: 90% !important;}
		</style>
		<script language="javascript" src="js/lib/jquery-ui-1.9.1.custom.min.js"></script>
		<script language="javascript">
			include.plugins('view');
			var sampleBrowser = new Browser("options");
			sampleBrowser.url = "propertyoptions/list.do";
			sampleBrowser.actionURL = "propertyoptions/*.do";
			sampleBrowser.itemElement = "option";
			sampleBrowser.itemTemplate = 'js/templates/poption.ejs';
			sampleBrowser.pager.selectors.pager = ".pgr";
			
			sampleBrowser.doAdd = function()
			{
				var str = window.prompt("Please, enter the name of the property option:");
				sampleBrowser.addOptions(str);
			}
			
			sampleBrowser.doAddmultiple = function()
			{
				$("#add-options").dialog("open");	
			}
			
			sampleBrowser.addOptions = function(str)
			{
				if (str &amp;&amp; str != "")
				{
					sampleBrowser.ajax.send({
						url: "propertyoptions/add.do",
						data: "property=<xsl:value-of select="//property/@id"/>&amp;name=" + URLEncode(str),
						success: function()
						{
							sampleBrowser.request(true);
						}
					});
				}
			}
			
			sampleBrowser.doPrompt_edit = function(link)
			{
				var str = window.prompt("Please, enter the new name of the property option:", sampleBrowser.currentEntity.name);
				if (str &amp;&amp; str != ""  &amp;&amp; str != sampleBrowser.currentEntity.name){
					sampleBrowser.ajax.send({
						url: "propertyoptions/edit.do",
						data: "id=" +  sampleBrowser.currentEntity.id + "&amp;name=" + URLEncode(str),
						success: function()
						{
							sampleBrowser.request(true);
						}
					});
				}
				
			}
			
			sampleBrowser.doEdit = undefined; // Disable the default edit handler from the upstream class
			
			sampleBrowser.onEditSuccess = function()
			{
				sampleBrowser.request(false);	
			}
			
			sampleBrowser.onDeleteSuccess = function()
			{
				sampleBrowser.deleteRecord();
			}
			
			$(document).ready(function() {
				sampleBrowser.initialize();
				$("#add-options").dialog({
					model: true,
					autoOpen: false,
					width: 450,
					buttons: {
						"Add options": function(){
							sampleBrowser.addOptions($("#multiple-options").val());
							$(this).dialog("close");
							sampleBrowser.request(true);
						},
						Cancel: function() {
							$( this ).dialog("close");
						}
					}
				});
			});
		</script>
		<title>
			Options for property <xsl:value-of select="//property/@name"/>
		</title>
		<table height="100%" width="100%">
			<tr><td class="itunes-up silver">
					<h1>Options for <a href="/properties/edit.do?id={//property/@id}" tab="Property profile"><xsl:value-of select="//property/@name"/></a></h1>
					This page allows you to view, add, modify and delete the options for a qualitative property "<b><xsl:value-of select="//property/@name"/></b>"
				</td></tr>
			<tr>
				<td class="itunes-right">
					<!--  <a class="rounded-button" action="add">Add a new option</a> -->
					<a class="rounded-button" action="addmultiple">Add multiple options</a>
					Filter options by name: <input type="text" name="name" filter="1"/><br/><br/>
					<input type="hidden" name="property" value="{//property/@id}" filter="1"/>
					<div class="pager-strip">
						<span><b class="showed">none</b> of <b class="total">none</b></span>
						<div class="pgr">
						</div>
					</div>
					
					<div id="Browser">
					</div>
					
					<div class="pager-strip">
						<span><b class="showed">none</b> of <b class="total">none</b></span>
						<div class="pgr">
						</div>
					</div>
				</td>
			</tr>
		</table>
		
		<div id="add-options" title="Add multiple options, one per line"> 
	       <textarea style="width: 400px; height: 300px;" id="multiple-options">
	       </textarea>
		</div>
	</xsl:template>
</xsl:stylesheet>