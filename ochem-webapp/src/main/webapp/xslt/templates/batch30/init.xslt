<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="2.0">
	<xsl:include href="../../helper.xslt" />
	<xsl:template name="content">
		<title>Batch upload (reloaded)</title>
		<link rel="stylesheet" type="text/css" href="css/batch.css" />
		<link rel="stylesheet" type="text/css" href="css/smoothness/jquery-ui-1.10.3.custom.min.css" />
		<script language="javascript" src="js/blocks/batch-upload.js"></script>
		<script language="javascript">
			$(document).ready(function(){setValidationHook()});
		</script>
		<table width="100%">
		<tr><td class="itunes-up silver">
			<img src="img/icons/batchupload.png"/>
			<h1><doc term="Batch+data+upload">Batch Upload 3.0 - File selection</doc></h1>
			Select the CSV, SDF or Excel file to upload multiple records to the database.
		</td></tr>
		<tr><td class="itunes-right big-padding ui-widget">
			<h1>Instructions</h1>
			<p>Here you have the possibility to upload data from an <b>excel sheet, sdf or csv </b>.<br/>Backslash \ is used as stereochemistry in cvs, which should not contain '\uffff' characters.<br/>The format of these data is strict, and can be viewed in <a href="documents/batch-sample-sci.xls">this sample</a> (scientific format) and <a href="documents/batch-sample-tech.xls">this sample</a> (technical format).</p>
			For more information, consider the wiki page that you can access by clicking on the wiki icon next to the title ("Batch upload browser").<br/>
			If you have difficulties uploading your data, feel free to drop us an e-mail at <a href="mailto:info@ochem.eu">info@ochem.eu.</a><br/><br/>
			<form action="batchupload30/init_submit.do?render-mode=popup" method="post" enctype="multipart/form-data">
			<fieldset>
				<legend>Select a file to upload</legend>
				Upload file <input type="file" name="file"/><br/>
				<xsl:if test="//bu/@state = 'error'">
					<b class="error">
						<xsl:value-of select="//bu/status"/>
					</b>				
				</xsl:if>
			</fieldset><br/>
			<fieldset>
				<legend>Settings</legend>
				<input type="checkbox" name="allow_pubchem" id="allow_pubchem">
					<xsl:if test="//uploadContext/@allowMoleculePubchemSearch = 'true'">
						<xsl:attribute name="checked">true</xsl:attribute>
					</xsl:if>
				</input><label for="allow_pubchem">Allow molecule lookup by name on PubChem</label><br/>
				<input type="checkbox" name="allow_pubmed" id="allow_pubmed">
					<xsl:if test="//uploadContext/@allowArticlePubmedSearch = 'true'">
						<xsl:attribute name="checked">true</xsl:attribute>
					</xsl:if>
				</input><label for="allow_pubmed">Allow article lookup by PMID on PubMed</label><br/>
				<input type="checkbox" name="hidden" id="hidden">
					<xsl:if test="//uploadContext/@hiddenByDefault = 'true'">
						<xsl:attribute name="checked">true</xsl:attribute>
					</xsl:if>
				</input><label for="hidden">Make the uploaded records hidden</label><br/>
			</fieldset><br/>
			<input type="submit" value="Upload" disabled="disabled"/>
			</form><br/><br/>
		</td></tr>
	</table>
	<div id="selenium-batch-upload-page-1" class="invisible"/>
	</xsl:template>
</xsl:stylesheet>