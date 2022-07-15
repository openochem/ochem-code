<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	
	<xsl:template name="content">
		<style type="text/css">
			
			SMALL {font-size: 70%;}
			#TagsBrowser-property, #TagsBrowser-molecule {float: left;}
			
		</style>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/browsers/tags-browser.js"></script>
		<script language="javascript">
			include.plugins('view');
		</script>
		<title>Area of interest</title>
		<table style="width: 100%;"><tr><td class="itunes-up">
		<img src="img/icons/areaofI.png"/>		
		<h1>Area of interest</h1>This page allows you to edit your area of interest and to filter all output according to it.</td></tr>
			<tr><td class="itunes-right">
			Following tags reflect your area of interest. Please provide adequate combination of tags to filter all the data.<br/>
			You will see only data, that matches ALL of the tags specified.
			<br/>
			<br/><br/><small class="scope-property">Tags related to properties: <a action="add">[add tag]</a></small>
			<div class="yellow scope-property">
				<div id="TagsBrowser-property">
				</div>&#160;
			</div>
			<small class="scope-molecule">Tags related to compounds: <a action="add">[add tag]</a></small>
			<div class="yellow scope-molecule">
				<div id="TagsBrowser-molecule">
				</div>&#160;
			</div>
		</td></tr>
		</table>
	</xsl:template>
</xsl:stylesheet>