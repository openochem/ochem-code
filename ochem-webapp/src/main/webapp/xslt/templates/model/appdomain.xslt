<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">

	<xsl:include href="../../helper.xslt" />
	
	<xsl:template name="content">
		<style type="text/css">
	  	</style>
	  	<script language="javascript" type="text/javascript" src="js/lib/excanvas.pack.js"></script>
	  	<script language="javascript" type="text/javascript" src="js/lib/jquery.flot-0.6.min.js"></script>
	  	<script language="javascript" type="text/javascript" src="js/commons/plotting.js"></script>
	  	<table width="100%">
	  		<tr>
	  			<td class="itunes-up">
	  				<h1>Applicability domain</h1>
	  				Estimates applicability of selected models to provided test set.
	  			</td>
	  		</tr>
	  		<tr><td class="itunes-right">
				<div id="demo" class="yui-navset"> 
				    <ul class="yui-nav"> 
				    	<xsl:for-each select="others/ad">
					    	<li class="selected">
					       		<a href="#tab{position()}"><em><xsl:value-of select="model/@name"/></em></a>
					    	</li>
				    	</xsl:for-each>
				    </ul> 
					<div class="yui-content">
						<xsl:for-each select="others/ad">
							<div>
								Applicability domain for model <i><xsl:value-of select="model/@name"/></i> 
								using dataset <i><xsl:value-of select="model/training-set/@name"/></i><br/>
								Distance to model used: <i><xsl:value-of select="@name"/></i>
								<div id="placeholder{position()}" style="width:500px;height:250px;"></div>
							</div>
						</xsl:for-each>
					</div>
				</div>
				</td>
			</tr>
		</table>
				
		<script language="javascript">
			var dt;
			var plot;
			var plots = new Array();
			<xsl:for-each select="others/ad">
				plot = new Plot();
				dt = ([
				<xsl:for-each select="training">
					[<xsl:value-of select="@dm"/>, <xsl:value-of select="@error"/>]
					<xsl:if test="position() != last()">,</xsl:if>
				</xsl:for-each>
				]);
				plot.addData(dt, "#C00");
				dt = ([
				<xsl:for-each select="test">
					[<xsl:value-of select="@dm"/>, <xsl:value-of select="@error"/>]
					<xsl:if test="position() != last()">,</xsl:if>
				</xsl:for-each>
				]);
				plot.addData(dt, "#0C0");
				var AD = undefined;
	       		<xsl:apply-templates select="ad-configuration" mode="adData"/>
	       		if (AD)
	       		if (AD.intervals)
	       		{
	       			plot.drawLine(plot.minX, AD.errors[0], AD.intervals[0], AD.errors[0]);
	       			for (var i = 0; i &lt; AD.intervals.length - 1; i++)
	       			{
	       				plot.drawLine(AD.intervals[i], AD.errors[i+1], AD.intervals[i+1], AD.errors[i+1]);
	       				plot.currentData.color = "#000";	
	       			}
	       		}
				
				
				plots.push(plot);
			</xsl:for-each>
			
			$(document).ready(function()
   			{
   				for (var i = 0; i &lt; plots.length; i++)
   					plots[i].render("#placeholder"+(i+1));
   				tabView = new YAHOO.widget.TabView("demo");	
   			});
		</script>
	</xsl:template>
	
	<xsl:template match="*" mode="adData">
		var AD = new Object();
		AD.intervals = [
			<xsl:for-each select="interval">
				<xsl:value-of select="."/><xsl:if test="position() != last()">, </xsl:if>
			</xsl:for-each>
			];
		AD.errors = [
			<xsl:for-each select="error">
				<xsl:value-of select="."/><xsl:if test="position() != last()">, </xsl:if>
			</xsl:for-each>
			];
	</xsl:template>
	
</xsl:stylesheet>