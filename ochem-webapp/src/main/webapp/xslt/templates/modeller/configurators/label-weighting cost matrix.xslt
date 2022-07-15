<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Model builder - Property weighting</title>
		<h1>Property weighting options</h1>
		<style type="text/css">
			.weights TD, .weights TH {padding: 6px; border-right: 3px solid white;}
			.weights TD {background-color: #FFE;}
			.weights TH {background-color: #FFE; font-weight: bold;}
			
			.cost-matrix INPUT {width: 40px; border: 1px solid gray; text-align: right; margin: 3px; font-size: 14pt;}
			
		</style>
		(N.B.! This option is under development!) This dialog allows to give more weight to a particular property or class while learning the model.
		This could be useful for training sets that contain multiple properties (e.g, "LogS" and "LogP") and/or multiple classes (e.g., "active" or "inactive").<br/>
		The default option is to give equal weight to each property/class.<br/>
		The property weigh for classes is a multiplication of the property and class weights.
		<br/>
		
		<table class="weights">
		<tr style="border-bottom: 1px solid black;">
			<th>Property/class name</th>
			<th>Weight</th>
		</tr>
		<xsl:for-each select="//others/labelWeighting/property-weight">
			<tr style="padding-top: 5px;">
				<td><xsl:value-of select="@name"/></td><td><input type="text" value="{@weight}" name="prop-{@name}"/></td>
			</tr>
			<xsl:for-each select="property-class">
				<tr>	
					<td style="padding-left: 20px;">•   <xsl:value-of select="@name"/></td><td><input type="text" value="{@weight}" name="prop-{../@name}$$$$class-{@name}"/></td>	
				</tr>
			</xsl:for-each>
		</xsl:for-each>
		</table>
	<!--	
		<xsl:if test="//others/labelWeighting/property-weight/property-class">
			<xsl:for-each select="//others/labelWeighting/property-weight">
				<br/><input type="checkbox" name="use-cm-{@name}" property="{position()}"/>Use cost matrix for <xsl:value-of select="@name"/><br/>
				<div id="cm-{position()}" class="block">
				<xsl:variable name="propertyName" select="@name"/>
				<table class="cost-matrix">
					<tr>
						<td></td>
						<xsl:for-each select="//others/labelWeighting/property-weight/property-class">
							<td><xsl:value-of select="@name"/></td>
						</xsl:for-each>	
					</tr>
					<xsl:for-each select="//others/labelWeighting/property-weight/property-class">
						<xsl:variable name="outerName" select="@name"/>
						<tr>
							<td><xsl:value-of select="@name"/></td>
							<xsl:for-each select="//others/labelWeighting/property-weight/property-class">
								<td><input type="text" name="cm-{$propertyName}--{$outerName}--{@name}" value="0.0"/></td>
							</xsl:for-each>	
						</tr>
					</xsl:for-each>
				</table>
				</div>
			</xsl:for-each>
		</xsl:if>
	-->	
		<script language="javascript">
			$("[property]").change(function(){
				$("#cm-" + $(this).attr("property")).setClass("invisible", !$(this).is(":checked"));
			});
			$(document).ready(function(){
				setTimeout("$('[property]').change()", 100);
			});
			
		</script>
	</xsl:template>
	
</xsl:stylesheet>