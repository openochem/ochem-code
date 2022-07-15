<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xs="http://www.w3.org/2001/XMLSchema" version="2.0">
	<xsl:include href="../helper.xslt" />
	<xsl:include href="inc/data-report.xslt" />
	<xsl:template name="content">
		<title>User contributions summary</title>
		
		<table height="100%" width="100%">
			<tr>
				<td class="itunes-up" colspan="2">
					<img src="img/icons/editor.png"/>
					<h1>User contributions summary</h1>
					This page provides a detailed report on the data and models contributed by a particular user
				</td>
			</tr>
			<tr>
				<td class="itunes-right">
					<xsl:call-template name="data-report"/>
					<script language="javascript">
						if (!getParams["introducer"])
							report.filters.setValue("introducer", "<xsl:value-of select="/model/session/user/@login"/>");
						$("input[name=introducer]").closest("div").remove();
					</script>
				</td>
			</tr>
		</table>
	</xsl:template>
</xsl:stylesheet>