<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:template match="model">
		<model web-root="{@web-root}" render-mode="{@render-mode}">
			<xsl:apply-templates/>
			<content>
				<xsl:call-template name="content"/>
			</content>
		</model>
	</xsl:template>
	
	<xsl:template match="@*|node()" priority="-1">
		<xsl:copy>
			<xsl:apply-templates select="@*|node()" />
		</xsl:copy>
	</xsl:template>
	
	
	<xsl:template match="log" mode="error">
		<div class="err notification">
			<xsl:for-each select="log-line[@type='error']">
				<xsl:apply-templates /><br />
			</xsl:for-each>
		</div>
	</xsl:template>
	
	<xsl:template match="log" mode="normal">
		<div class="notification">
			<xsl:for-each select="log-line[@type='normalhigh']">
				<b><xsl:apply-templates /></b><br/>
			</xsl:for-each>		
			<xsl:for-each select="log-line[@type='normal']">
				<xsl:apply-templates /><br/>
			</xsl:for-each>
		</div>
	</xsl:template>
	
	<xsl:template match="log" mode="stats">
		<div class="stats notification">
			<xsl:for-each select="log-line[@type='stats']">
				<xsl:apply-templates /><br/>
			</xsl:for-each>
		</div>
	</xsl:template>
	
	<xsl:template name="area-of-interest">
		<xsl:choose>
			<xsl:when test="//others/tag">
				<div class="right">
					<doc term="Area+of+interest">Area of your interest:</doc>
					<xsl:for-each select="//others/tag">
<!-- 						<a href="wikipage/action.do?name={@name}" class="greeny" tab="Wiki"><xsl:value-of select="@name"/></a>,  -->
						<a class="greeny"><xsl:value-of select="@name"/></a>,
					</xsl:for-each><br/>
					<a tab="Select area of interest" href="globalfilters/show.do">[change]</a>
				</div>	
			</xsl:when>
			<xsl:otherwise>
				<div class="right">
					<doc term="Area+of+interest">Area of your interest:</doc>
						no tags selected
					<a tab="Select area of interest" href="globalfilters/show.do">[change]</a>
				</div>		
			</xsl:otherwise>
		</xsl:choose>
	</xsl:template>
	
</xsl:stylesheet>