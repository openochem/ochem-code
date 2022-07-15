<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />

	<xsl:template name="content">
		<script language="javascript" src="js/commons/midnighter.js"></script>
		<table width="100%">
			<iframe width="1px" height="1px" src="basket/doindices.do?id={basket/@id}"/>
			<tr><td class="itunes-up">
				<h1>Calculation of indices</h1>
			</td></tr>
			<tr><td class="itunes-right">
				Please, wait a bit. Your indices are being calculated and will be available in excel format.
			</td></tr>
		</table>
	</xsl:template>
</xsl:stylesheet>