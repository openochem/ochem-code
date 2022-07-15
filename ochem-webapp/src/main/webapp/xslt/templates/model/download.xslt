<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />

	<xsl:template name="content">
		<script language="javascript" src="js/commons/midnighter.js"></script>
		<table width="100%">
			<tr><td class="itunes-up">
				<h1>Download descriptors and model statistics</h1>
			</td></tr>
			<tr><td class="itunes-right">
				<form method="post">
					<input type="hidden" name="id" value="{model/@id}"/>
					<input type="checkbox" name="download.hyperlinks"/><label>Create hyperlinks</label><br/>
					<input type="checkbox" name="download.names"/><label>Molecule names</label><br/>
					<input type="checkbox" name="download.values" checked="checked"/><label>Observed and predicted values</label><br/>
					<input type="checkbox" name="download.ad" checked="checked"/><label>Applicability domain data</label><br/>
					<input type="checkbox" name="download.descriptors"/><label>Molecular descriptors</label><br/><br/>
					<input type="submit" name="submit" value="Download XLS file"/>
				</form>				
			</td></tr>
		</table>
		
		
	</xsl:template>
</xsl:stylesheet>