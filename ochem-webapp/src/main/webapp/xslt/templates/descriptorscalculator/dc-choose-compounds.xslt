<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	
	<xsl:include href="../../helper.xslt" />
	<xsl:include href="../model/inc/select-compounds.xslt" />
	
	<xsl:template name="content">
	
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<table height="100%" width="100%">
			<tr><td class="itunes-up silver">
				<img src="img/icons/mol-descriptors.png"/>
				<h1>DescriptorsCalculator: Provide the molecules</h1>
				Provide the molecules, choose the descriptor types and download the descriptor values as an Excel, CSV or SDF file
				</td></tr>
			<tr>
				<td class="itunes-right">
				
				<xsl:if test="//message">
					<div class="warning">
						<xsl:value-of select="//message/message"/>
					</div>
				</xsl:if>
				
				
				<form method="post" enctype="multipart/form-data" action="descriptorscalculator/submitCompounds.do">
					<b>Provide the compounds to calculate molecular descriptors for:</b><br/><br/>
					<xsl:call-template name="select-compounds"/>
					<div class="formsubmit">
						<input type="submit" name="next" value="Next&gt;&gt;"/>
					</div>	
				</form>
				</td>
			</tr>
		</table>
		
		</xsl:template>
</xsl:stylesheet>
