<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />

	<xsl:template name="content">
		<title>Unit editor</title>
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
		<script language="javascript" src="js/blocks/unit-edit.js"></script>
		<h1 class="popup-header">Unit editor
		<i>Add new unit or edit exiting unit</i></h1>
		<table width="780">
		<tr class="EditForm"><td>
			<input type="hidden" name="id" value="{unit/@id}" send="1" filter="1"/>
		Name:<input type="text" onkeyup="return checkUnitName(this)" name="name"  value="{unit/@name}" send="1"/>
		<div id="unitName" class="message">(min. 2 characters and max. 40 characters)</div>
		<xsl:if test="unit/unit-property">
					<br/>Following properties have been used with this unit: 
					<xsl:for-each select="unit/unit-property">
						<a href="properties/edit.do?id={@id}" tab="Unit profile">
							<xsl:if test="@name = ''">
								empty unit
							</xsl:if>
							<xsl:value-of select="@name"/></a> (<a href="epbrowser/show.do?unit={../@id}&amp;property={@id}" tab="Records of property+unit"><xsl:value-of select="@timesUsed"/></a>), 
					</xsl:for-each>
				</xsl:if>
		</td></tr>
		<tr class="EditForm" id="Units"><td>
			<div>
				System of units:
					<select name="category" send="1" onlynew="1">
						<xsl:for-each select="others/unitcategory">
							<option value="{@id}">
								<xsl:if test="//unit/unitCategory/@id=@id">
									<xsl:attribute name="selected">selected</xsl:attribute>
								</xsl:if>
								<xsl:attribute name="defaultUnit"><xsl:value-of select="@default-unit"/></xsl:attribute>
								<xsl:value-of select="@name"/>
							</option>
						</xsl:for-each>
					</select>
					<span id="def-unit"></span>
			</div>
		</td></tr>
		<tr class="EditForm">
			<td>
				Conversion:<br/>
				To default: <input send="1" name="to-default-conversion" value="{unit/toDefaultConversion}"/>
				From default: <input send="1" name="from-default-conversion" value="{unit/fromDefaultConversion}"/>
				<a action="verify" restrict="form">[verify]</a>
			</td>
		</tr>
		<tr class="EditForm"><td>
			Description:<br/>
			<textarea name="description" onkeyup="return checkUnitDesc(this)" send="1"
				style="width: 100%; height: 80px;">
				<xsl:value-of select="unit/description" />
			</textarea>
			<div id="unitDesc" class="message">(min. 50 characters)</div>
		</td></tr>
		</table>
		<div class="EditForm popup-footer">
			<a action="edit" restrict="form">save</a>
			<a href="javascript:window.closeTab();">cancel</a>
		</div>
	</xsl:template>
</xsl:stylesheet>