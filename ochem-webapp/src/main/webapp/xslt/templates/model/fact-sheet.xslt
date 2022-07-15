<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />

	<xsl:template name="content">
		<title>Model fact sheet</title>
		<style type="text/css">
			
		</style>
		
		<h1>Model fact sheet: <xsl:value-of select="/model/model/modelMappings[1]/property/@name"/></h1>
		
		<p>
		Full model profile is available at <a href="https://ochem.eu/model/{/model/model/publicId}">https://ochem.eu/model/<xsl:value-of select="/model/model/publicId"/></a>
		</p>
		
		<h2>Predicted endpoint</h2>
		<p>
		<b><xsl:value-of select="/model/model/modelMappings[1]/property/@name"/></b>. <xsl:value-of select="/model/model/modelMappings[1]/property/description"/>
		</p>
		
		<h2>Dataset description</h2>
		<p>
		  <xsl:call-template name="break">
		    <xsl:with-param name="text" select="/model/model/training-set/description" />
		  </xsl:call-template>
		</p>
		
		<h2>Data preprocessing</h2>
		
		<h2>Molecular descriptors</h2>
		
		<h2>Machine learning method</h2>
		
		<h2>Validation</h2>
		
	</xsl:template>
	
	<xsl:template name="break">
	  <xsl:param name="text" select="."/>
	  <xsl:choose>
	    <xsl:when test="contains($text, '&#xa;')">
	      <xsl:value-of select="substring-before($text, '&#xa;')"/>
	      <br/>
	      <xsl:call-template name="break">
	        <xsl:with-param 
	          name="text" 
	          select="substring-after($text, '&#xa;')"
	        />
	      </xsl:call-template>
	    </xsl:when>
	    <xsl:otherwise>
	      <xsl:value-of select="$text"/>
	    </xsl:otherwise>
	  </xsl:choose>
</xsl:template>
	
</xsl:stylesheet>