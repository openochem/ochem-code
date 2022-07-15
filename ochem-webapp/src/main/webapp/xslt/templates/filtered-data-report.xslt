<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xs="http://www.w3.org/2001/XMLSchema" version="2.0">
	<xsl:include href="../helper.xslt" />
	<xsl:include href="inc/data-report.xslt" />
	<xsl:template name="content">
		<title>OCHEM data report</title>
		
		<table height="100%" width="100%">
			<tr>
				<td class="itunes-up" colspan="2">
					<img src="img/icons/analytics-48.png"/>
					<h1>OCHEM data report</h1>
					A detailed report on a particular subset of OCHEM data records.
				</td>
			</tr>
			<tr>
				<td class="itunes-right">
					<xsl:call-template name="data-report"/>
					<script language="javascript">
						if (getParams["property"])
							report.filters.setValue("property", getParams["property"]);
					</script>
				</td>
			</tr>
		</table>
	</xsl:template>
</xsl:stylesheet>