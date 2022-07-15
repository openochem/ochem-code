<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />

	<xsl:template name="content">
		<style type="text/css">
			TABLE.torefactor TD, TABLE.torefactor TH {background-color: #FFD; padding: 10px; text-align: left; border-bottom: 1px solid white;}
			TABLE.torefactor TH {font-style: bold; border-bottom: 1px solid black;}
			TR.ready TD {background-color: #c0e793;}
		</style>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<title>Model templates</title>
		
		<table height="100%" width="100%">
			<tr><td class="itunes-up">
				<h1>Model templates</h1>
				
				</td></tr>
			<tr>
				<td class="itunes-right">
				In the list below, you can find the <b>predefined model configuration templates</b>. You can use these templates and skip the detailed model configuration.<br/><br/>
					<div class="pager-strip">
						<span><b class="showed">none</b> of <b class="total">none</b></span>
						<div class="pgr">
						</div>
					</div>
					<div id="Browser"></div>
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
					<pre id="configuration-xml">
					</pre>
				</div>
			</div>
		</div>

		
		<script language="javascript">
			
			include.plugins('view');
			var TemplateBrowser = function()
			{
				var self = this;
				this.controller = "modeltemplate";
				Browser.call(this);
				this.itemTemplate = "js/templates/models/model-template.ejs";
				this.itemElement = "model-template";
				this.pager.selectors.pager = ".pgr";
			}
			
			var templateBrowser = new TemplateBrowser();
			
			templateBrowser.doConfiguration = function(link)
				{
					templateBrowser.xmlDialog.setContent('<img src="img/roller_transparent.gif"/>');
					templateBrowser.xmlDialog.show();
					$.ajax({
						url: "modeltemplate/getXml.do",
						data: "out=json&amp;id=" + templateBrowser.currentEntity.id,
						dataType: "json",
						success: function(response){
							var xml = response.message.message;
							xml = xml.replace(/&gt;/g, "&amp;gt;");
							xml = xml.replace(/&lt;/g, "&amp;lt;");
							templateBrowser.xmlDialog.setContent(xml);	
						},
						error: function(){
							templateBrowser.xmlDialog.setContent("Failed to load configuration XML.");
						}
					});
				}
			
			$(document).ready(function(){
				templateBrowser.initialize();
				
				templateBrowser.xmlDialog = new YAHOO.widget.Dialog("xmlDialog", { 
			    	width:"700px", 
					fixedcenter:true, 
					modal:true, 
					visible:false 
			    });
			    templateBrowser.xmlDialog.render();
			    
			    templateBrowser.xmlDialog.setContent = function(content)
			    {
			    	$("#configuration-xml").html(content);
			    }
			});
			
		</script>
	</xsl:template>
</xsl:stylesheet>