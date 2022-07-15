<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	<xsl:template name="content">
		<style type="text/css">
			#content {margin: 40px;}
			#header H1 {font-size: 25px; font-family: Georgia; color: #666; padding: 30px 10px; text-shadow: 1px 1px #AAA;}
			#header { padding-left: 40px; background-color: #EEE;}
			#content P {margin-bottom: 10px;}
			#content PRE {padding: 20px; background-color: #EEF;}
			#content IMG {float: left; margin-right: 20px;}
			#friendly-message {text-align: center; font-size: 200%;}
		</style>
		
		<div id="header">
			<h1>We got an issue</h1>
		</div>
		<div id="content">
			<xsl:choose>
			<xsl:when test="//message/@context = 'friendly'">
				<div class="friendly-message">
					<xsl:value-of select="//message/@title"/>
				</div>
			</xsl:when>
			<xsl:otherwise>
				<p>Dear user!</p>
				<p>
					Unfortunately, an error occurred while serving your request.<br/>
					The error may be temporary, so you can try to <a href="javascript: window.location.reload()">retry your request</a>.</p>
				<p>
					If the problem persists, please feel free to report the case to <a href="http://ochem.eu/support/" target="_blank">our bug tracking system</a>.
					<br/>We are dedicated to constantly improve OCHEM and we appreciate any feedback!
				</p>
				
				<p><br/>Thank you!</p>
				<p>Sincerely,<br/> the OCHEM team</p>
				
				<br/><br style="clear: both"/>P.S.: If you want to report a bug, please describe steps that you have performed and include the following information about the failure:<br/>
				<pre>
<b>Time of the failure:</b>  <xsl:value-of select="//message/time"/><br/>
<b>Exception stacktrace:</b>
	
				<xsl:value-of select="//message/message"/>
				</pre>
			</xsl:otherwise>
			</xsl:choose>
		</div>
	</xsl:template>
	
</xsl:stylesheet>
