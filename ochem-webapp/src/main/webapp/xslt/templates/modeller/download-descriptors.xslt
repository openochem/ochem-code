<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />

	<xsl:template name="content">
		<title>Model builder - Descriptors</title>
		<style type="text/css">
			.warning {background-color: #FCC;}
		</style>
		
		<table width="100%">
			<tr><td class="itunes-up">
				<h1>Download your descriptors</h1>
			</td></tr>
			<tr><td class="itunes-right">
				<xsl:if test = "param[@key='calculate'] = 'false'">
				<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
				<script language="javascript" src="js/commons/browser.js"></script>
				</xsl:if>
				Descriptors have been calculated for the set <a href="epbrowser/show.do?basket-select={model/training-set/@id}" tab="Basket compounds"><xsl:value-of select="model/training-set/@name"/></a>
				<a href="basket/edit.do?id={model/training-set/@id}" tab="Basket details">[details]</a><br/>
				Calculated descriptors: <xsl:value-of select="//param[@key='descriptors']"/> 
				<br/><br/>
				<a href="descriptormapper/descriptors.do?type=csv&amp;model={model/@id}"><img src="img/icons/csv_icon.gif"/>Download descriptors as CSV</a><br/>
				<a href="descriptormapper/descriptors.do?type=xls&amp;model={model/@id}"><img src="img/icons/xls.gif"/>Download descriptors as XLS</a><br/>
				<a href="descriptormapper/descriptors.do?type=sdf&amp;model={model/@id}">Download descriptors as SDF.GZ</a><br/>
				</td>
			</tr>
		</table>
	</xsl:template>
</xsl:stylesheet>