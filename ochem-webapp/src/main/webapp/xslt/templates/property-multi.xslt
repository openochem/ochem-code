<?xml version="1.0" encoding="UTF-8"?>
<!-- Universal template to track the status of long-running server side operations -->
<xsl:stylesheet
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	<xsl:template name="content">

		<style type="text/css">
			DIV.comp-block-login
			{
			padding: 30px;
			width: 770px;
			}
			DIV.comp-block-login DIV {padding: 20px; border-bottom: 1px solid #999; margin: 5px;}
			DIV.comp-block-login TABLE TD {padding: 5px 15px 5px 0px;}
			TABLE.comp-block-login
			{
			border: outset 1px #AAA;
			border-spacing: 30pt 10pt;
			background-color:#FFF;
			margin: 10px 10px 10px 10px;
			position: relative;

			}
			.button-link {display: block; margin-top: 15px;}
			.comp-block-login INPUT {border: 1px solid #555; width: 200px;}
			TABLE.comp-block-login TD {padding: 15px 10px 15px 10px;}
			.comp-block-login H1 {font-size: 200%; margin-bottom: 13px; border-bottom: 1px solid
			black;}
			DIV.comp-block-login H2 {font-weight: bold; color: #000; margin-bottom: 6px; font-size:
			130%;}
		</style>

		<title>Create multiple properties</title>
		<h1>Create mutiple properties with high/low options:</h1>

		<form action="properties/addmany.do" method="post"
			name="mainform" target="_parent" accept-charset="UTF-8,ISO-8859-1">

			<input type="hidden" name="render-mode" value="redirect" />
			<table>
				<tr>
					<td>
						<td class="inner formscope" colspan="2">
							Properties:
							<br />
							<textarea name="ids" send="1"
								style="width: 100%; height: 80px;">
							</textarea>
							<div class="message">(add new properties, one per line)</div>
						</td>
					</td>
				</tr>
				<tr>
					<td colspan="1" height="30" valign="bottom">
						<a onclick="document.forms[0].submit(); return false;" href="#"
							class="button-link">submit</a>
						<input type="image" value="submitname"
							src="submit-button.gif" width="1" height="1" border="0"
							alt="SUBMIT!" name="image" class="invisible" />
					</td>
				</tr>
			</table>
		</form>

	</xsl:template>
</xsl:stylesheet>
