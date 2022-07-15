<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	
	<xsl:include href="../../helper.xslt" />
	<xsl:include href="../modeller/inc/descriptor-blocks.xslt" />
	
	<xsl:template name="content">
	
		<table height="100%" width="100%">
			<tr><td class="itunes-up">
				<img src="img/icons/mol-descriptors.png"/>
				<h1>DescriptorsCalculator: Select the descriptors for comparison</h1>
				</td></tr>
			<tr>
				<td class="itunes-right">
				Please, select the descriptor types that you would like to calculate:
				
				<form method="post" enctype="multipart/form-data" action="descriptorscalculator/descriptorsSubmit.do">
					<xsl:call-template name="descriptor-blocks"/>
					<br/>
					
					<div class="formsubmit">
						<input type="submit" name="next" value="Next&gt;&gt;"/>
					</div>	
					
				</form>
				
				</td>
			</tr>
		</table>
		
	</xsl:template>
</xsl:stylesheet>
