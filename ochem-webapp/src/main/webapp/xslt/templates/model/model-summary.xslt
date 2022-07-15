<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />

	<xsl:template match="model-attachment">
		<b>Descriptors:</b><br/>
		<xsl:for-each select="configuration/descriptors/types">
			<xsl:value-of select="type"/><br/>
		</xsl:for-each><br/>
		
		<b>Descriptors filtering:</b><br/>
			Filter descriptors with pairwise correlation more than <xsl:value-of select="selection/correlationThreshold"/>
			Filter descriptors with less than <xsl:value-of select="selection/numDifferentValues"/> unique values
	</xsl:template>
</xsl:stylesheet>