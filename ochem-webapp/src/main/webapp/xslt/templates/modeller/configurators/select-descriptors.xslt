<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Model builder - Descriptors</title>
		<style type="text/css">
			.warning {background-color: #FCC;}
			.outlined {border: 1px solid gray; padding: 15px; position: relative; margin: 2px; float: left; height: 3em;}
			.outlined SPAN {position: absolute; left: 10px; top: -0.5em; background-color: white;}
		</style>
		<h1>Select pre-filtered descriptors for the model</h1>
			<div id="BrowserScope">
					Filter by name: 
					<input type="text" filter="1" name="name" value=""/>
					<input type="checkbox" name="selected" filter="1">
						<xsl:if test="//param[@key='show-selected'] = 'true'">
							<xsl:attribute name="checked">checked</xsl:attribute>
						</xsl:if>
					</input> Show only selected descriptors
				<div class="pager-strip">
					<b class="showed">none</b> of <b class="total">none</b>
				</div>
				
				<div id="pager">
				</div>
				
				<div id="Browser">
				</div>
				
				<div class="pager-strip">
					<b class="showed">none</b> of <b class="total">none</b>
				</div>
				<br/><br/>
			</div>
			
			<script language="javascript">
			include.plugins('view');
			var sampleBrowser = new UnifiedBrowser();
			sampleBrowser.scope = "#BrowserScope";
			sampleBrowser.url = "descriptors/list.do";
			sampleBrowser.actionURL = "descriptors/action.do";
			sampleBrowser.itemTemplate = "js/templates/descriptor.ejs";
			sampleBrowser.itemElement = "descriptor";
			sampleBrowser.onToggleSuccess = function(xml)
			{
				var newValue = (this.currentEntity.selected == "true") ? "false" : "true";
				this.currentEntity.selected = newValue;
				this.currentBlock.find('img[name="checked"]').attr('src', 
					(newValue == "true") ? "img/icons/checked.gif" : "img/icons/unchecked.gif");
			}
			
			$(document).ready(function() {sampleBrowser.initialize();});
		</script>
	</xsl:template>
</xsl:stylesheet>