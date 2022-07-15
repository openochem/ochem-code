<?xml version="1.0" encoding="UTF-8"?>


<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	
	<xsl:include href="../../helper.xslt" />
	<xsl:include href="../model/inc/data-preprocessing-commons.xslt" />
	
	<xsl:template name="content">
	
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<table height="100%" width="100%">
			<tr><td class="itunes-up silver">
				<img src="img/icons/mol-descriptors.png"/>
				<h1>DescriptorsCalculator: Structure preprocessing</h1>
				Select the structure preprocessing options
				</td></tr>
			<tr>
				<td class="itunes-right">
				
				<form method="post" enctype="multipart/form-data" action="descriptorscalculator/dataPreprocessingSubmit.do">
					<xsl:call-template name="data-preprocessing"/>
					<div class="formsubmit">
						<input type="submit" name="next" value="Next&gt;&gt;"/>
					</div>	
				</form>
				</td>
			</tr>
		</table>
		
		</xsl:template>
</xsl:stylesheet>