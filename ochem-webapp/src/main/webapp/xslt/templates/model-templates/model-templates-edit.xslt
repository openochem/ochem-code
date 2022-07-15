<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />

	<xsl:template name="content">
		<style type="text/css">
			
		</style>
		<title>Model templates</title>
		
		<table height="100%" width="100%">
			<tr><td class="itunes-up">
				<h1>Model template editor</h1>
				Advanced users only! Requires knowledge of OCHEM XML format
				</td></tr>
			<tr>
				<td class="itunes-right">
					<div class="invisible info" id="saved">
						Your template has been successfully saved.
					</div>
					<form action="modeltemplate/save.do" method="post">
						<input type="hidden" name="id" value="{/model/model-template/@id}"/>
						Model template name: <input type="text" value="{/model/model-template/@name}" name="name" size="40"/><br/><br/>
						<textarea style="width: 800px; height: 600px;" name="xml"><xsl:value-of select="/model/model-template/full-xml"/></textarea>
						<br/><input type="submit" value="Save"/>
					</form>
				</td>
			</tr>
		</table>
		<script language="javascript">
			if (getParams["saved"])
				$("#saved").removeClass("invisible");
		</script>
	</xsl:template>
</xsl:stylesheet>