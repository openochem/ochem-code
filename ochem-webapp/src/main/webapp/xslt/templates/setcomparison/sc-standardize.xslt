<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />
	<xsl:include href="../model/inc/data-preprocessing-commons.xslt" />
	
	<xsl:template name="content">
		<table height="100%" width="100%">
		<tr>
			<td class="itunes-up">
			<img src="img/icons/compare.png"/>
			<h1><span class="setcompare">SetCompare</span>: Select molecule standardization options</h1>
			</td>
		</tr>
		<tr>
			<td class="itunes-right">
			<form method="post" enctype="multipart/form-data" action="setcomparison/standardizeSubmit.do">
				<xsl:call-template name="data-preprocessing"/>
				<br style="clear: both;"/>
				<input type="submit" name="submit" value="Next &gt;&gt;" class="button"/><br/><br/>
			</form>
			</td>
		</tr>
		</table>
	</xsl:template>
	
</xsl:stylesheet>