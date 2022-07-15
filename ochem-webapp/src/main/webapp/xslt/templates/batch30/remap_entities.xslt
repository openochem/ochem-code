<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="2.0">
<xsl:include href="../../helper.xslt" />
<xsl:template name="content">
	<title>Batch upload (reloaded)</title>
	<link rel="stylesheet" type="text/css" href="css/batch.css" />
	<link rel="stylesheet" type="text/css" href="css/smoothness/jquery-ui-1.10.3.custom.min.css" />
	<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
	<script language="javascript" src="js/blocks/batch-upload-30.js"></script>
	<script language="javascript">
		remapping = new PropertyRemappingActionable();
		$(document).ready(function(){
			
		});
	</script>
	<table width="100%">
	<tr><td class="itunes-up silver">
		<img src="img/icons/batchupload.png"/>
		<h1><doc term="Batch+data+upload">Batch Upload 3.0 - Entity remapping</doc></h1>
		Review and remap the properties, conditions, units, articles and baskets involved in the data upload
	</td></tr>
	<tr><td class="itunes-right big-padding ui-widget">
		<h1>Database entities remapping</h1>
		<form action="batchupload30/remap_entities_submit.do?render-mode=popup" method="post">
		<div class="clear">
		<xsl:for-each select="//remapping/properties">
		 	<xsl:variable name="pprefix">property<xsl:value-of select="position()"/></xsl:variable>
		 	
			<div class="container p500">
				Property: <xsl:choose><xsl:when test="@name = 'Dummy'">molecules without property</xsl:when><xsl:otherwise><a action="property"><xsl:attribute name="name"><xsl:value-of select="$pprefix"/></xsl:attribute><xsl:value-of select="@name"/></a></xsl:otherwise></xsl:choose><br/>
				<input type="hidden" value="{@name}"><xsl:attribute name="name"><xsl:value-of select="$pprefix"/></xsl:attribute></input>
				<xsl:call-template name="message-block"/>
				<xsl:if test="options">
					<div class="container">
						Options<br/>
						<ul>
						<xsl:for-each select="options">
							<xsl:variable name="oprefix"><xsl:value-of select="$pprefix"/>_option<xsl:value-of select="position()"/></xsl:variable>
							<li><xsl:attribute name="name"><xsl:value-of select="$oprefix"/></xsl:attribute><xsl:value-of select="@name"/></li>
							<input type="hidden" value="{@name}"><xsl:attribute name="name"><xsl:value-of select="$oprefix"/></xsl:attribute></input>
							<xsl:call-template name="message-block"/>
						</xsl:for-each>
						</ul>
					</div>
				</xsl:if>	
				<xsl:if test="@name != 'Dummy'">
				<xsl:if test="units">
					<div class="container">
						Values<br/>
						<ul>
						<xsl:for-each select="units">
							<xsl:variable name="uprefix"><xsl:value-of select="$pprefix"/>_unit<xsl:value-of select="position()"/></xsl:variable>
							<li>Unit: <a action="unit"><xsl:attribute name="name"><xsl:value-of select="$uprefix"/></xsl:attribute><xsl:value-of select="@name"/></a>, min value: <xsl:value-of select="minValue"/>, max value:  <xsl:value-of select="maxValue"/></li>
							<input type="hidden" value="{@name}"><xsl:attribute name="name"><xsl:value-of select="$uprefix"/></xsl:attribute></input>
							<xsl:call-template name="message-block"/>
						</xsl:for-each>
						</ul>
					</div>
				</xsl:if>
				</xsl:if>
				<xsl:for-each select="conditions">
					<xsl:variable name="cprefix"><xsl:value-of select="$pprefix"/>_condition<xsl:value-of select="position()"/></xsl:variable>
					<div class="container">
						Condition: <a action="condition"><xsl:attribute name="name"><xsl:value-of select="$cprefix"/></xsl:attribute><xsl:value-of select="@name"/></a><br/>
						<input type="hidden" value="{@name}"><xsl:attribute name="name"><xsl:value-of select="$cprefix"/></xsl:attribute></input>
						<xsl:call-template name="message-block"/>
						<xsl:if test="options">
							<div class="container">
								Options<br/>
								<ul>
								<xsl:for-each select="options">
									<xsl:variable name="coprefix"><xsl:value-of select="$cprefix"/>_option<xsl:value-of select="position()"/></xsl:variable>
									<li><xsl:attribute name="name"><xsl:value-of select="$coprefix"/></xsl:attribute><xsl:value-of select="@name"/></li>
									<input type="hidden" value="{@name}"><xsl:attribute name="name"><xsl:value-of select="$coprefix"/></xsl:attribute></input>
									<xsl:call-template name="message-block"/>
								</xsl:for-each>
								</ul>
							</div>
						</xsl:if>
						<xsl:if test="units">
							<div class="container">
								Values<br/>		
								<ul>						
								<xsl:for-each select="units">
									<xsl:variable name="cuprefix"><xsl:value-of select="$cprefix"/>_unit<xsl:value-of select="position()"/></xsl:variable>
									<li>Unit: <a action="unit"><xsl:attribute name="name"><xsl:value-of select="$cuprefix"/></xsl:attribute><xsl:value-of select="@name"/></a>, min value: <xsl:value-of select="minValue"/>, max value:  <xsl:value-of select="maxValue"/></li>
									<input type="hidden" value="{@name}"><xsl:attribute name="name"><xsl:value-of select="$cuprefix"/></xsl:attribute></input>
									<xsl:call-template name="message-block"/>
								</xsl:for-each>
								</ul>
							</div>
						</xsl:if>
					</div>
				</xsl:for-each>			
			</div>
		</xsl:for-each>
		<xsl:for-each select="//remapping/articles">
			<xsl:variable name="aprefix">article<xsl:value-of select="position()"/></xsl:variable>
			<div class="container p500">
				Article: 
				<a action="article"><xsl:attribute name="name"><xsl:value-of select="$aprefix"/></xsl:attribute>
				<xsl:choose>
					<xsl:when test="@name = 'A1551'">
						Unpublished
					</xsl:when>
					<xsl:otherwise>
						<xsl:value-of select="@name"/>
					</xsl:otherwise>
				</xsl:choose>				
				</a><br/>
				<input type="hidden" value="{@name}"><xsl:attribute name="name"><xsl:value-of select="$aprefix"/></xsl:attribute></input>
				<xsl:call-template name="message-block"/>
			</div>
		</xsl:for-each>
		<xsl:for-each select="//remapping/baskets">
			<xsl:variable name="bprefix">basket<xsl:value-of select="position()"/></xsl:variable>
			<div class="container p500">
				Molecule set: <a action="basket"><xsl:attribute name="name"><xsl:value-of select="$bprefix"/></xsl:attribute><xsl:value-of select="@name"/></a><br/>
				<input type="hidden" value="{@name}"><xsl:attribute name="name"><xsl:value-of select="$bprefix"/></xsl:attribute></input>
				<xsl:call-template name="message-block"/>
			</div>
		</xsl:for-each>		
		</div>
		<xsl:if test="count(//sheetSchema/messages) > 0">
			<div class="container p500">
			<ul>
			<xsl:for-each select="//sheetSchema/messages">
			<li><xsl:value-of select="."/></li>
			</xsl:for-each>
			</ul>
			</div>
		</xsl:if>
		<xsl:if test="count(//*[type = 'error']) > 0">
			<div class="container p500"><b class="error">There are errors about the uploded sheet. Please resolve them before you continue.</b></div>
		</xsl:if>		
		<xsl:if test="count(//*[type = 'warning']) > 0">
			<div class="container p500">
				<input type="checkbox" name="ignorewarnings" id="ignorewarnings"/><label for="ignorewarnings">Ignore warnings</label><br/>
				<b class="warn">There are warnings in the uploded sheet. Please check the above checkbox to confirm that you would like to ignore them.</b>
			</div>
		</xsl:if>
		<input type="submit" name="submit" value="submit"/>
		</form>
	</td></tr>
	<tr><td class="bubutton itunes-right big-padding ui-widget right">
		<a href="batchupload30/cancel.do">Cancel Batch Upload</a>
		<a href="batchupload30/report.do">Download Excel file</a>
	</td></tr>
	
	</table>
	<div id="selenium-batch-upload-page-3" class="invisible"/>
</xsl:template>	
<xsl:template name="message-block">
	<xsl:for-each select="messages">
		<xsl:choose>
			<xsl:when test="type = 'error'"><div><b class="error">Error: <xsl:value-of select="message"/></b></div></xsl:when>
			<xsl:when test="type = 'warning'"><div><b class="warn">Warning: <xsl:value-of select="message"/></b></div></xsl:when>
			<xsl:when test="type = 'notice'"><div><b class="notice">Notice: <xsl:value-of select="message"/></b></div></xsl:when>
			<xsl:otherwise><div><b class="notice">Notice: <xsl:value-of select="message"/></b></div></xsl:otherwise>
		</xsl:choose>
	</xsl:for-each>
</xsl:template>
</xsl:stylesheet>