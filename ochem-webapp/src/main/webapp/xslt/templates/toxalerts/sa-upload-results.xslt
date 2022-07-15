<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />
	
	<xsl:template name="content">
		<style type="text/css">
			.highlights {padding: 10px 10px; background-color: #EEF; border: 1px solid #77A; font-size: 130%; margin-top: 10px;} 
		</style>
		<table height="100%" width="100%">
			<tr><td class="itunes-up">
				<img src="img/icons/sa.png"/>
				<h1>Structure alert upload</h1>
				</td></tr>
			<tr>
				<td class="itunes-right">
					Your file has been uploaded.
					<br/>
					Highlights of the event are:<br/>
					<div class="highlights">
						<pre><xsl:value-of select="//message/message"/></pre>
					</div>
					<br/><br/>
					<a class="fb-button" href="alerts/show.do?render-mode=popup">Go back to the alerts browser</a>
				</td>
			</tr>
		</table>
		</xsl:template>
</xsl:stylesheet>
