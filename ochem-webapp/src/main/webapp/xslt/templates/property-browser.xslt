<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />

	<xsl:template name="content">
		<style type="text/css">
			#filters { margin: 8px 0px;}
			
		</style>
		<xsl:variable name="entity">
			<xsl:choose>
				<xsl:when test="//param[@key='condition']">condition</xsl:when>
				<xsl:otherwise>property</xsl:otherwise>
			</xsl:choose>
		</xsl:variable>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/browsers/property-browser.js"></script>
		<script language="javascript">
			include.plugins('view');
			var sampleBrowser = new PropertyBrowser("properties", "property");
			sampleBrowser.listenEvent("items_loaded", function(){
				if (!$("input[name=show-counts]").is(":checked"))
					$(".count-records-link").html("Show records");
			});
			$(document).ready(function() {sampleBrowser.initialize();});
		</script>
		<title>
			<xsl:choose>
				<xsl:when test="//param[@key='condition']">Condition browser</xsl:when>
				<xsl:otherwise>Property browser</xsl:otherwise>
			</xsl:choose>
		</title>
		<table height="100%" width="100%">
			<tr><td class="itunes-up">
				<xsl:call-template name="area-of-interest"/>
				<xsl:choose>
					<xsl:when test="//param[@key='condition']">
						<img src="img/icons/condition.png"/>
						<h1><doc term="Property+browser#Propertybrowser-conditionBrowserBrowserofconditions">Conditions browser</doc></h1>
						Please search condition database before creating new ones
					</xsl:when>
					<xsl:otherwise>
						<img src="img/icons/property.png"/>
						<h1><doc term="Property+browser#Propertybrowser-propertyBrowserBrowserofproperties">Properties browser</doc></h1>
						Please search property database before creating new ones
					</xsl:otherwise>
				</xsl:choose>
				</td></tr>
			<tr>
				<td class="itunes-right" a="{$entity}">
				<a action="edit" title="Add new {$entity}" class="rounded-button">
					<xsl:if test="//param[@key='condition']">
						<xsl:attribute name="query">condition=true</xsl:attribute>
					</xsl:if>				
					Create new <xsl:value-of select="$entity"/>
				</a>
				<br class="push"/>
				<div id="filters">
				Type part of name to filter:
					<input id="namefilter" type="text" name="query" filter="1"/>
					
					<a href="javascript:sampleBrowser.request(true)" class="rounded-button">
						search
					</a>
					
					<img src="img/icons/awaiting-approval.png"/>
						Properties from other users:&#160;&#160;
						<select name="approval-status" filter="1">
							<option value="only-approved">Only approved properties</option>
							<option value="all">All properties</option>
							<option value="only-awaiting-approval">Only awaiting approval</option>
						</select>
					
					
					
					
					&#160;&#160;&#160;<input type="checkbox" filter="1" name="show-counts"/>&#160;Calculate counts
					<div style="float: right">
						
					</div>
				</div>
					<div class="pager-strip">
						<span><b class="showed">none</b> of <b class="total">none</b></span>
						<div class="pgr">
						</div>
					</div>
					
					<div id="Browser">
					</div>
					
					<div class="pager-strip">
						<span><b class="showed">none</b> of <b class="total">none</b></span>
						<div class="pgr">
						</div>
					</div>
				</td>
			</tr>
		</table>
	</xsl:template>
</xsl:stylesheet>