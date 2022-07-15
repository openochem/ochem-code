<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />

	<xsl:template name="content">
		<title>Tag editor</title>
		<style type="text/css">
			INPUT, TEXTAREA {border: 1px solid black; padding: 2px;}
			TD {padding: 2px 2px 2px 2px;}
			TD{border: 1px solid #999; padding: 10px;}
			INPUT, TEXTAREA {border: 1px solid black; padding: 2px 2px 2px 2px; margin: 1px 1px 1px 1px;}
			.right {align:right;}
			.hidden {display: none;}
			.message {font-size: 11px; font-style: italic;}
		</style>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/blocks/tag-edit.js"></script>
		<h1 class="popup-header">Tag editor
		<i>Add new tag or edit exiting tag</i></h1>
		<table width="720">
		<tr class="EditForm"><td>
			<input type="hidden" name="id" value="{tag/@id}" send="1" filter="1"/>
		Name:<input type="text" onkeyup="return checkTagName(this)" name="name"  value="{tag/@name}" send="1"/>
		<div id="tagName" class="message">(min. 2 characters and max. 40 characters)</div>
		</td></tr>
		<tr class="EditForm" id="Tags"><td>
			<div>
				Tag type:
				<select name="type" send="1" disabled="disabled">
					<option value="property">
						<xsl:if test="tag/@type = 'property'">
							<xsl:attribute name="selected">selected</xsl:attribute>
						</xsl:if>
						Properties
					</option>
					<option value="molecule">
						<xsl:if test="tag/@type = 'molecule'">
							<xsl:attribute name="selected">selected</xsl:attribute>
						</xsl:if>
						Molecules
					</option>
				</select>
					
				<br/>
				<input type="checkbox" name="isPublic" send="1">
					<xsl:if test="//tag/@isPublic = 'true'">
						<xsl:attribute name="checked">checked</xsl:attribute>
					</xsl:if>	
				</input>This tag is public (accessible by everybody)
				<br/><input type="checkbox" name="showInBrowser" send="1">
					<xsl:if test="//tag/@showInBrowser = 'true'">
						<xsl:attribute name="checked">checked</xsl:attribute>
					</xsl:if>
				</input>Show this tag in browser of experimental properties
			</div>
		</td></tr>
		<xsl:if test="tag/@id">
		<tr class="invisible" id="upload">
			<td>
				You can assign this tag to a set of molecules by uploading an SDF file:<br/>
				<form action="tags/action.do" id="uploadform" target="uploadframe" enctype="multipart/form-data" method="post">
					<input type="hidden" name="id" value="{tag/@id}"/>
					<input type="hidden" name="action" value="uploadfile"/>
					<span id="progress" class="invisible"><img src="img/roller_small.gif"/></span><input type="file" id="file" name="file"/>
				</form>
				<iframe name="uploadframe" id="uploadframe" width="1" height="1" frameborder="0" onload="ajaxForm.iframeLoaded();"></iframe>
			</td>
		</tr>
		</xsl:if>
		<tr class="EditForm"><td>
			Description:<br/>
			<textarea name="description" onkeyup="return checkTagDesc(this)" send="1"
				style="width: 100%; height: 80px;">
				<xsl:value-of select="tag/description" />
			</textarea>
			<div id="tagDesc" class="message">(min. 50 characters)</div>
		</td></tr>
		</table>
		<div class="EditForm popup-footer">
			<a action="edit" restrict="form">save</a>
			<a href="javascript:window.closeTab();">cancel</a>
		</div>
		
		<div id="waitingDialog"> 
	    <div class="hd">Please wait</div> 
	    <div class="bd" style="text-align: center;"> 
	        Your file with molecules is currently being processed<br/>
	        It may take a while.<br/>
	        <img src="img/roller_small.gif"/> 
	    </div> 
	</div>
		
	</xsl:template>
</xsl:stylesheet>