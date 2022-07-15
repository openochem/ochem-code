<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />

	<xsl:template name="content">
		<script language="javascript" src="js/commons/midnighter.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<table width="100%">
			<tr><td class="itunes-up">
				<h1>Model applier</h1>
			</td></tr>
			<tr><td class="itunes-right">
				<form action="modelapplier/provideCenters.do" method="post">
					<div class="formsubmit">
						<input type="button" name="submit" value="&lt;&lt;Back" onclick="location.href='model/select.do';"/>
						<input type="submit" id="next" name="next" value="Start calculations&gt;&gt;"/>
					</div>	
				</form>
			</td></tr>
		</table>
	</xsl:template>
</xsl:stylesheet>