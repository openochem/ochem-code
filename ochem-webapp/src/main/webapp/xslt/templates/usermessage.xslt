<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	<xsl:template name="content">
	<div><br/><p>
		<xsl:apply-templates select="log" mode="error"/>
		<xsl:apply-templates select="log" mode="normal">
		</xsl:apply-templates>
		</p>
	</div>
	</xsl:template>
</xsl:stylesheet>
