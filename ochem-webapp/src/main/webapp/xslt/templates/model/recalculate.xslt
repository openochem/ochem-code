<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />

	<xsl:template name="content">
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<style type="text/css">
			.itunes-right {text-align: center; vertical-align: center;}
		</style>
		
		<title>Recalculate model data</title>
		
		<table width="100%">
			<tr><td class="itunes-up">
				<h1>Provide molecule</h1>
			</td></tr>
			<tr><td class="itunes-right">
				Status: <span id="status">Initializing...</span>		
			</td></tr>
		</table>
	</xsl:template>	
</xsl:stylesheet>