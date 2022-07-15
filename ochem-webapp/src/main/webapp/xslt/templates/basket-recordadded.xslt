<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />

	<xsl:template name="content">
		<script language="javascript" src="js/commons/midnighter.js"></script>
		<table width="100%">
			<tr><td class="itunes-up">
				<h1>Basket successfully updated</h1>
			</td></tr>
			<tr>
				<td class="itunes-right">
				The basket <b><xsl:value-of select="basket/@name"/></b> has been successfully updated 
				and now contains <b><xsl:value-of select="basket/@size"/></b> records (<xsl:value-of select="basket/excludedCount"/> thereof excluded)
				<br/><br/>
				<div class="popup-footer">
				<a href="javascript:window.closeTab()" class="button-link">close</a>
				</div>
				</td>
			</tr>
		</table>
	</xsl:template>
</xsl:stylesheet>