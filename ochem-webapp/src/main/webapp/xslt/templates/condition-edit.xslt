<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />

	<xsl:template name="content">
		<title>Condition editor</title>
		<style type="text/css">
			INPUT, TEXTAREA {border: 1px solid black; padding: 2px;}
			TD {padding: 2px 2px 2px 2px;}
			TD{border: 1px solid #999; padding: 10px;}
			INPUT, TEXTAREA {border: 1px solid black; padding: 2px 2px 2px 2px; margin: 1px 1px 1px 1px;}
			.right{align:right;}
			.hidden {display: none;}
		</style>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script language="javascript" src="js/blocks/condition-edit.js"></script>
		<title>Condition editor</title>
		<h1 class="popup-header">Condition editor
			<i>Add new Condition or edit exiting Condition</i></h1>
		<table>
		<tr class="EditForm"><td>
		<input type="hidden" name="id" value="{condition/@id}" send="1" filter="1"/>
		Name:<input id="name" type="text" name="name" value="{condition/@name}" send="1"/>
		<input type="checkbox" name="qualitive" send="1" onlynew="1">
			<xsl:if test="condition/@qualitive='true'">
				<xsl:attribute name="checked">checked</xsl:attribute>
			</xsl:if>
		</input>
		Is a qualitive condition
		</td></tr>
		<tr class="EditForm"><td>
			Description:<br/>
			<textarea name="description" send="1"
				style="width: 100%; height: 80px;">
				<xsl:value-of select="condition/description" />
			</textarea>
		</td></tr>
		<tr class="EditForm"><td>
			<div id="Units">
				Unit system:
					<select name="category" send="1" onlynew="1">
						<xsl:for-each select="others/unitcategory">
							<option value="{@id}">
								<xsl:if test="//condition/unitCategory/@id=@id">
									<xsl:attribute name="selected">selected</xsl:attribute>
								</xsl:if>
								<xsl:value-of select="@name"/>
							</option>
						</xsl:for-each>
					</select>
				<label>Default unit:</label> 
					<select name="unit" send="1">
						<xsl:for-each select="condition/unitCategory/unit">
							<option value="{@id}">
							<xsl:if test="@id=//condition/defaultUnit/@id">
								<xsl:attribute name="selected">selected</xsl:attribute>
							</xsl:if>
							<xsl:value-of select="@name"/>
							</option>
						</xsl:for-each>
					</select>
			</div>
		</td></tr>
		<tr id="Options">
		<td>
			<a action="add">[add class]</a>
			<div id="Browser">
			</div>
		</td>
		</tr>
		</table>
		<div class="EditForm popup-footer">
			<a action="myedit">save</a>
			<a href="javascript:window.closeTab();">cancel</a>
		</div>
	</xsl:template>
	
</xsl:stylesheet>