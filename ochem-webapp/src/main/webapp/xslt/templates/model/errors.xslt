<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />

	<xsl:template name="content">
		<title>Error records browser</title>
		<style type="text/css">
			.conditions {color: green; font-size: 80%; float: right; clear: right; width: 500px; text-align: right;}
			.article-data {font-size: 8pt; margin-top: 7px; margin-left: 0px; margin-bottom: 7px;}
			img {vertical-align: middle;}
		</style>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/browsers/modelrecord-browser.js"></script>
		<script language="javascript">		
			include.plugins('view');
			var errorBrowser = new ModelRecordBrowser();
			errorBrowser.type = "error";
			$(document).ready(function() 
			{
				errorBrowser.initialize();
				errorBrowser.filters.setFromUrl();
			});	
			var recalculated = false;
			<xsl:if test="/model/param = 'true'">
				recalculated = true;
			</xsl:if>
			
			function openInBrowser()
			{
				var url = "epbrowser/show.do?modelerrors=<xsl:value-of select="model/@id"/>&amp;modelError=" + errorBrowser.filters.getValue("errorTitle");
				if (recalculated)
					url += "&amp;recalculated=1";
				openTab("Model errors", url);
			}
		</script>
		<table width="99%">
			<tr><td class="itunes-up">
				<h1><doc term="Records+with+errors">Records with errors</doc></h1>
				A list of the records, not processed by the modeling workflow because of errors
			</td></tr>
			<tr>
				<td class="errorbrowser itunes-right">
					<input type="hidden" filter="1" name="id"/>
					<input type="hidden" filter="1" name="mm_id"/>
					<input type="hidden" filter="1" name="recalculated"/>
					<a href="#" onclick="openInBrowser(); return false;">
						[Open these records in the browser of experimental measurements]
					</a>			
					<br/>
					Filter by the error message:<input type="text" filter="1" name="errorTitle"/>
		    		<div class="pager-strip">
					<b class="showed">none</b> of <b class="total">none</b>
					</div>
					<div id="pager">
					</div>		
					<div id="ErrorBrowser">
					</div>	
				</td>
			</tr>		
			<tr>
				<td align="right">
					<a href="javascript:window.closeTab()">close</a>
				</td>
			</tr>
		</table>
	</xsl:template>
	
</xsl:stylesheet>