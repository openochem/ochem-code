<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />

	<xsl:template name="content">
		<xsl:variable name="entity">unit</xsl:variable>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/browsers/unit-browser.js"></script>
		<script language="javascript">
			include.plugins('view');
			var sampleBrowser = new UnitBrowser("unit", "unit");
			$(document).ready(function() {sampleBrowser.initialize();});
		</script>
		<title>Unit browser</title>
		<table height="100%" width="100%">
			<tr><td class="itunes-up">
				<xsl:call-template name="area-of-interest"/>
					<img src="img/icons/property.png"/>
					<h1>Units browser</h1>
					Please search unit database before creating new ones
				</td></tr>
			<tr>
				<td class="itunes-right" a="{$entity}">
					Filter by name:
					<input type="hidden" name="lightweight" value="1" filter="1"/>
					<input type="text" name="query" filter="1"/>
					and System of units:
					<select name="category" filter="1" onlynew="1">
						<xsl:if test="count(others/unitcategory) &gt; 1">
							<option value="0">Select All</option>
						</xsl:if>	
						<xsl:for-each select="others/unitcategory">
						<option value="{@id}">
							<xsl:if test="@selected = 'true'">
								<xsl:attribute name="selected">selected</xsl:attribute>
							</xsl:if>
							<xsl:value-of select="@name"/>
						</option>
						</xsl:for-each>
					</select>
					<a href="javascript:sampleBrowser.request(true)">
						[search]
					</a>
					<a action="edit" title="Add new {$entity}">
						[Create new <img src="img/icons/new.gif"/>]
					</a>
					
					<div class="pager-strip">
						<span><b class="showed">none</b> of <b class="total">none</b></span>
						<div id="pager" class="pgr">
						</div>
					</div>
					
					<div id="Browser">
					</div>
					
					<div class="pager-strip">
						<span><b class="showed">none</b> of <b class="total">none</b></span>
						<div id="pager" class="pgr">
						</div>
					</div>
					
				</td>
			</tr>
		</table>
	</xsl:template>
</xsl:stylesheet>