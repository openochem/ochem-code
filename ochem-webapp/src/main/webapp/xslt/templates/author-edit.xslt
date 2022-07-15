<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	
	<xsl:template name="content">
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/blocks/author-edit.js"></script>
		<style type="text/css">
			.outer TD {}
			.inner TD {padding: 2px 2px 2px 2px;}
			TD.inner {border: 1px solid #999; padding: 10px;}
			INPUT, TEXTAREA {border: 1px solid black; padding: 2px 2px 2px 2px; margin: 1px 1px 1px 1px;}
			.outer IMG {border: 1px solid black;}
			.inner TR TD:first-child {width: 100px;}
			.hidden {display: none;}
			.right {align: right;}
		</style>
		<h1 class="popup-header">Author editor
			<i>Add new author or edit exiting author</i></h1>
		<table class="outer">
			<tr>
				<td class="inner formscope">
					<input type="hidden" name="id" value="{author/@id}" send="1"/>
					
					<table>
						<tr><td colspan="2"><b>Introducing Author</b></td></tr>
						<tr><td><b>FirstName</b></td><td><input type="text" name="n-firstname" send="1" value="{author/first_name}" class="w600"/></td></tr>
						<tr><td><b>Initials</b></td><td><input type="text" name="n-initials" send="1" value="{author/initials}" class="w50"/></td></tr>
						<tr><td><b>LastName</b></td><td><input type="text" name="n-lastname" send="1" value="{author/last_name}" class="w600"/></td></tr>
					</table>			
				</td>
			</tr>
		</table>
		<div class="formscope right"><input type="button" action="edit" value="save"/></div>
	</xsl:template>
</xsl:stylesheet>