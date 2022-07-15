<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="2.0">
<xsl:include href="../../helper.xslt" />
<xsl:template name="content">
	<title>Batch upload (reloaded)</title>
	<link rel="stylesheet" type="text/css" href="css/batch.css" />
	<link rel="stylesheet" type="text/css" href="css/main.css" />
	<link rel="stylesheet" type="text/css" href="css/smoothness/jquery-ui-1.10.3.custom.min.css" />
	<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
	<script language="javascript" src="js/commons/browser.js"></script>
	<script language="javascript" src="js/blocks/batch-upload.js"></script>
	<script language="javascript">
		include.plugins('view');
		browser = new BatchUploadBuBrowser();
		$(document).ready(function(){
			browser.initialize();
		});
	</script>
	<table width="100%">
	<tr><td class="itunes-up">
		<img src="img/icons/batchupload.png"/>
		<h1><doc term="Batch+data+upload">Your batch upload history</doc></h1>
		Upload of data via CSV, SDF or Excel files.
	</td></tr>
	<tr><td class="itunes-right big-padding ui-widget">
		<h1>Batch upload history</h1>
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
	</td></tr>
	<tr><td class="bubutton itunes-right big-padding ui-widget">
		<a href="batchupload30/cancel.do?render-mode=popup">Batch Upload</a>
	</td></tr>
	</table>
</xsl:template>	
</xsl:stylesheet>