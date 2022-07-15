<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	
	<xsl:include href="../../helper.xslt" />
	<xsl:include href="../modeller/inc/descriptor-blocks.xslt" />
	
	<xsl:template name="content">
	
		<table height="100%" width="100%">
			<tr><td class="itunes-up">
				<img src="img/icons/compare.png"/>
				<h1><span class="setcompare">SetCompare</span>: Select the descriptors for comparison</h1>
				</td></tr>
			<tr>
				<td class="itunes-right">
				The SetCompare utility allows you to compare two sets of molecules based on their structural features.<br/>
				Please, provide the descriptors which will be used to compare the two sets.<br/><br/>
				<form method="post" enctype="multipart/form-data" action="setcomparison/descriptorsSubmit.do">
					<xsl:call-template name="descriptor-blocks"/>
					<input type="submit" name="submit" value="Next &gt;&gt;" class="button"/><br/><br/>
				</form>
				
				</td>
			</tr>
		</table>
	
		
	</xsl:template>
</xsl:stylesheet>
