<?xml version="1.0" encoding="UTF-8"?>


<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	
	<xsl:include href="../../helper.xslt" />
	<xsl:include href="../model/inc/structure-optimisation-commons.xslt" />
	
	<xsl:template name="content">
		<table height="100%" width="100%">
			<tr><td class="itunes-up silver">
				<img src="img/icons/mol-descriptors.png"/>
				<h1>DescriptorsCalculator: Structure optimisation</h1>
				Select the structure optimisation options
				</td></tr>
			<tr>
				<td class="itunes-right">
				
				<form method="post" enctype="multipart/form-data" action="descriptorscalculator/optimisationSubmit.do">
					<xsl:call-template name="structure-optimisation"/>
					<div class="formsubmit">
						<input type="button" name="prev" value="&lt;&lt;Back" onclick="history.go(-1);"/>
						<input type="submit" name="next" value="Next&gt;&gt;"/>
					</div>	
				</form>
				</td>
			</tr>
		</table>
		
		</xsl:template>
</xsl:stylesheet>