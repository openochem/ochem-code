<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />

	<xsl:template name="content">
		<script language="javascript" src="js/commons/midnighter.js"></script>
		<table width="100%">
			<tr><td class="itunes-up">
				<h1>Model creator</h1>
			</td></tr>
			<tr><td class="itunes-right">
				<h1>Your model has been saved</h1>
				Thank you for your cooperation.
				<br/><br/>
				Your next possible actions are:
				<ul>
					<li><a href="modelapplier/apply.do?render-mode=popup&amp;model={model/@id}">Apply your model</a></li>
					<li><a href="model/profile.do?id={model/@id}">View your model's properties</a></li>
				</ul>
			</td></tr>
		</table>
	</xsl:template>
</xsl:stylesheet>