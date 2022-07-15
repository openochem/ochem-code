<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="2.0">
<xsl:include href="../../helper.xslt" />
<xsl:template name="content">
	<title>Batch upload (reloaded)</title>
	<link rel="stylesheet" type="text/css" href="css/batch.css" />
	<link rel="stylesheet" type="text/css" href="css/smoothness/jquery-ui-1.10.3.custom.min.css" />
	<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
	<script language="javascript" src="js/commons/browser.js"></script>
	<script language="javascript" src="js/blocks/batch-upload.js"></script>
	<table width="100%">
	<tr><td class="itunes-up silver">
		<img src="img/icons/batchupload.png"/>
		<h1><doc term="Batch+data+upload">Batch upload 3.0 - finished</doc></h1>
		Your upload has been finished
	</td></tr>
	<tr><td class="itunes-right big-padding ui-widget">
		<h1>Batch upload results</h1>
		Batch upload is finished. You can download the <a href="batchupload30/report.do">detailed upload report</a>.
		<table class="bu-item compact-item p500">
		<tr><td colspan="2"><b>Summary:</b></td></tr>
		<tr><td>All rows in the sheet</td><td>Count: <b><xsl:value-of select="//previewSummary/@total"/></b></td></tr>
		<tr><td></td><td></td></tr>
		<xsl:for-each select="//previewSummary/summary">
		<tr><td>Status: <i><xsl:value-of select="@status"/>, <xsl:value-of select="@uploadStatus"/></i></td><td name="{@status},{@uploadStatus}">Count: <b><xsl:value-of select="@count"/></b></td></tr>
		</xsl:for-each>
		</table>
		N.B.! If basket already existed, new records were added to it.  In the case of duplicates, the final number of records can be smaller.
	</td></tr>
	<tr><td class="bubutton itunes-right big-padding ui-widget right">
		<a href="batchupload30/cancel.do">New Batch Upload</a>
		<a href="batchupload30/report.do">Download Excel file</a>
	</td></tr>
	</table>
	<div id="selenium-batch-upload-page-5" class="invisible"/>
</xsl:template>	
</xsl:stylesheet>