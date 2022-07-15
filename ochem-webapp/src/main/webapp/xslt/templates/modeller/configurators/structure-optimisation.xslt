<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	<xsl:include href="../../model/inc/structure-optimisation-commons.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Model builder - structure optimisation options</title>
		<xsl:call-template name="structure-optimisation"/>
	</xsl:template>
	
</xsl:stylesheet>