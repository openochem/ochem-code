<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="2.0">
<xsl:include href="../../helper.xslt" />
<xsl:template name="content">
	<title>Batch upload (reloaded)</title>
	<link rel="stylesheet" type="text/css" href="css/batch.css" />
	<link rel="stylesheet" type="text/css" href="css/main.css" />
	<link rel="stylesheet" type="text/css" href="css/smoothness/jquery-ui-1.10.3.custom.min.css" />
	<script type="text/javascript" src="js/lib/jquery-ui-1.10.3.custom.min.js" />
	<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
	<script language="javascript" src="js/commons/browser.js"></script>
	<script language="javascript" src="js/blocks/batch-upload-30.js"></script>
	<script language="javascript">
		include.plugins('view');
		browser = new BatchUploadBrowser();
		browser.options.highlight_row = "no";
		$(document).ready(function(){
			browser.dialogs();
			browser.initialize();
		});
	</script>
	<table width="100%">
	<tr><td class="itunes-up silver">
		<img src="img/icons/batchupload.png"/>
		<h1><doc term="Batch+data+upload">Batch upload 3.0 - records preview</doc></h1>
		Preview the records you are about to upload, select the desired actions
	</td></tr>
	<tr><td class="itunes-right big-padding ui-widget">
		<h1>Batch upload preview browser</h1>
		<table class="bu-item compact-item p500">
		<tr><td colspan="2"><b>Summary:</b></td></tr>
		<tr><td>All rows in the sheet</td><td>Count: <b><xsl:value-of select="//previewSummary/@total"/></b></td></tr>
		<tr><td></td><td></td></tr>
		<xsl:for-each select="//previewSummary/summary">
		<tr><td>Status: <i><xsl:value-of select="@status"/>, <xsl:value-of select="@uploadStatus"/></i></td><td>Count: <b><xsl:value-of select="@count"/></b></td></tr>
		</xsl:for-each>
		</table>
		<table width="100%"><tr><td class="nowrap">
		Filter by row number:
		<input type="text" name="rownum" filter="1"/>
		and row type:
		
		<select name="type" filter="1">
			<option value="all">all</option>
			<option value="valid">valid</option>
			<option value="invalid">all invalid</option>
			<option value="error">errors</option>
			<option value="fatal_error">fatal errors</option>
			<option value="warning">warnings</option>
			<option value="duplicate_internal">internal duplicates</option>
			<option value="duplicate_external">external duplicates</option>
		</select>
		</td>
		<td class="narrow nowrap bubutton">
		<a action="batchmenu">Batch operations</a>	
		</td></tr></table>
		<div class="pager-strip">
			<span><b class="showed">none</b> of <b class="total">none</b></span>
			<div id="pager" class="pgr">
			</div>
		</div>
		<div id="Browser">
		</div>
		<div class="pager-strip">
			<span><b class="showed">none</b> of <b class="total">none</b></span>
			<div id="pager" class="pgr">
			</div>
		</div>
		<form action="batchupload30/browser_submit.do" method="post">
			<input type="submit" name="submit" value="Proceed with upload"/>
		</form>
	</td></tr>
	<tr><td class="bubutton itunes-right big-padding ui-widget right">
		<a href="batchupload30/cancel.do">Cancel Batch Upload</a>
		<a href="batchupload30/report.do">Download Excel file</a>
	</td></tr>
	</table>
	<div id="batchoperations" title="Batch errors" style="display: none;">
		<table>
		<xsl:for-each select="//previewSummary/summary">
		<tr class="narrow nowrap bubutton">
			<td><xsl:value-of select="@status"/>:</td> 
			<xsl:choose>
			 <xsl:when test="@status = 'warning' or @status = 'valid'">
			 	<td><a action="batch" operation="save"><xsl:attribute name="type"><xsl:value-of select="@status"/></xsl:attribute>save</a></td>
			 	<td><a action="batch" operation="skip"><xsl:attribute name="type"><xsl:value-of select="@status"/></xsl:attribute>skip</a></td>
			 	<td></td>
			 </xsl:when>
			 <xsl:when test="@status = 'duplicate_internal' or @status = 'duplicate_external'">
			 	<td><a action="batch" operation="save"><xsl:attribute name="type"><xsl:value-of select="@status"/></xsl:attribute>save</a></td>
			 	<td><a action="batch" operation="skip"><xsl:attribute name="type"><xsl:value-of select="@status"/></xsl:attribute>skip</a></td>
			 	<td><a action="batch" operation="merge"><xsl:attribute name="type"><xsl:value-of select="@status"/></xsl:attribute>merge</a></td>
			 	<td><a action="batch" operation="put_original_to_basket"><xsl:attribute name="type"><xsl:value-of select="@status"/></xsl:attribute>use existing records</a></td>
			 </xsl:when>
			 <xsl:when test="@status = 'error'">
			 	<td><a action="batch" operation="save"><xsl:attribute name="type"><xsl:value-of select="@status"/></xsl:attribute>save as error</a></td>
			 	<td><a action="batch" operation="skip"><xsl:attribute name="type"><xsl:value-of select="@status"/></xsl:attribute>skip</a></td>
			 	<td></td>
			 </xsl:when>
			 <xsl:when test="@status = 'fatal_error' or @status = 'undefined'">
			 	<td></td>
			 	<td><a action="batch" operation="skip"><xsl:attribute name="type"><xsl:value-of select="@status"/></xsl:attribute>skip</a></td>
			 	<td></td>
			 </xsl:when>
			</xsl:choose>
		</tr>
		</xsl:for-each>
		</table>			
	</div>
	<div id="selenium-batch-upload-page-4" class="invisible"/>
</xsl:template>	
</xsl:stylesheet>