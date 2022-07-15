<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	
	<xsl:template name="content">
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/blocks/journal-edit.js"></script>
		<style type="text/css">
			.outer TD {}
			.inner TD {padding: 2px 2px 2px 2px;}
			TD.inner {border: 1px solid #999; padding: 10px;}
			INPUT, TEXTAREA {border: 1px solid black; padding: 2px 2px 2px 2px; margin: 1px 1px 1px 1px;}
			.outer IMG {border: 1px solid black;}
			.inner TR TD:first-child {width: 100px;}
			.hidden {display: none;}
		</style>
		<title>Edit journal</title>
		<h1 class="popup-header">Journal editor
			<i>Add new journal or edit exiting journal</i></h1>
		<table class="outer">
			<tr>
				<td class="inner formscope">
					<input type="hidden" name="id" value="{journal/@id}" send="1"/>
					<table>
						<tr><td><b>Title</b></td><td><input type="text" name="n-title" send="1" value="{journal/title}" class="w600"/></td></tr>
						<tr class="journalrow"><td><b>Abbreviation</b></td><td><input type="text" name="n-abbreviation" send="1" value="{journal/abbreviation}" class="w600"/></td></tr>
					</table>			
				</td>
			</tr>
			<tr>
				<td class="inner formscope">
				<table>
					<tr><td><b>Publisher:</b></td><td><input type="text" name="n-publisher" send="1" value="{journal/publisher}" class="w600"/></td></tr>
					<tr id="issn">
						<td><b>ISSN:</b></td><td><input type="text" name="n-issn" n-disable="1" send="1" value="{journal/issn}" class="w100"/>
						<input type="button" action="reload" value="reload"/>
						<br/><div style="font-size: 11px; text-align: left;" n-hidden="1" >(XXXX-XXXX)</div>
						</td>
						</tr>
					<tr><td><b>WebSite:</b></td><td><input type="text" name="n-link" send="1" value="{journal/link}" class="w600"/></td></tr>
				</table>
				</td>
			</tr>
		</table>
		<div class="formscope popup-footer">
			<a action="edit" edit="journal">save</a>
			<a href="javascript:window.closeTab();">cancel</a>
		</div>
	</xsl:template>
</xsl:stylesheet>