<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />
	<xsl:template name="content">
	<style type="text/css">
	</style>
	<title>Import the model configuration</title>
		<table width="100%">
			<tr><td class="itunes-up">
				<h1><div class="not-uplodable">Import model configuration</div></h1>
				<p class="not-uploadable">
				Create a new model by loading the template from an XML file
				</p>
			</td></tr>
			<tr><td class="itunes-right">
				<form action="modelconfigurator/importTemplateSubmit.do" method="post" enctype="multipart/form-data">
					<b>OCHEM model configuration files</b> (.ochem.xml) contain information about the training and validation sets,<br/>
					the training method, the modelled property, molecular descriptors and their configurations.<br/><br/>
					If you have such a file, you can use it to create a new model.<br/><br/>
					Select the configuration file: <input type="file" name="file"/><br/><br/>
					<input type="submit" name="submit" value="Import the template and create the model" class="button"/>
				</form>
			</td></tr>
		</table>
	</xsl:template>
</xsl:stylesheet>