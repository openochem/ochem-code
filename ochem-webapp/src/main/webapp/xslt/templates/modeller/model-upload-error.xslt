<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />
	<xsl:template name="content">
		<style type="text/css">
			.message {margin: 40px;}
			.message P {margin-bottom: 10px;}
			.message PRE {padding: 20px; background-color: #EEF;}
			.message IMG {float: left; margin-right: 20px;}
			.message {
				align:center;
				width: 70%;
			}
			.errorlist
			{
				background-color: #EEEEEE;
			}
			
			.errorlist li
			{
				font-family: Courier, sans-serif;
				padding: 2px;
			}
		</style>
		<script language="javascript">
		</script>		
			<table width="100%">
			<tr><td class="itunes-up silver">
				<h1><doc term="Uploading+a+stub+QSAR+model">Upload a stub model</doc></h1>
				<p>Select the training and validation sets, and upload the predicted values from an external file</p>
			</td></tr>
			<tr><td class="itunes-right">
				<div class="message">
					<p>
						Unfortunately, an <b>error</b> occurred during your model upload process.<br/>
						Please address the errors below and <a href="/modelupload/show.do">reattempt your model upload</a>.
					</p>
					<ul class="errorlist">
					<xsl:for-each select="//modelUploadResult/errors">
						<li><xsl:value-of select="."/></li>
					</xsl:for-each>
					</ul>
				</div>
			</td></tr>
			</table>
	</xsl:template>
	
</xsl:stylesheet>
