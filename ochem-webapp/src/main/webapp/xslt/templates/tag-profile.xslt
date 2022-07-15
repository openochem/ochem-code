<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />

	<xsl:template name="content">
		<title>Tag profile: <xsl:value-of select="tag/@name"/></title>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/browsers/property-browser.js"></script>
		<script language="javascript">
			include.plugins('view');
			var sampleBrowser = new PropertyBrowser("properties", "property");
			sampleBrowser.filters.useUrlParameters = false;
			$(document).ready(function() {sampleBrowser.initialize();});
		</script>
		<table height="100%" width="100%">
			<tr>
				<td class="itunes-up" colspan="2">
					<xsl:call-template name="area-of-interest"/>
					<h1>
						Tag profile: <xsl:value-of select="tag/@name"/>
					</h1>
					List of properties related to this tag
				</td>
			</tr>
			<tr>
				<td class="itunes-right">
					<!-- wiki_info -->
<!-- 					<a tab="Wiki" href="wikipage/action.do?name={tag/@name}" border="0"><img src="img/icons/wiki.gif"/></a> -->
<!-- 					You can visit or edit our <a tab="Wiki" href="wikipage/action.do?name={tag/@name}" border="0">wiki page</a> for tag <i>"<xsl:value-of select="tag/@name"/>"</i> -->
					<br/>Manage tags of each individual property by visiting its profile.
					<br/><br/>
					<input type="hidden" name="tag" value="{tag/@id}" filter="1"/>
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
				</td>
			</tr>
		</table>
	</xsl:template>
</xsl:stylesheet>