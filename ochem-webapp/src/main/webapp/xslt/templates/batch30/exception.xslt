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
		</style>
		<script language="javascript">
			$(document).ready(function(){
				$("#toggle").click(function(){
					var link = $("#toggle");
					var content = $("#strace");
					if (link.text() == "[&gt;&gt;]")
					{
						link.html("[&lt;&lt;]")
						content.removeClass("invisible");
					} else
					{
						link.html("[&gt;&gt;]")
						content.addClass("invisible");
					}
				});
			});
		</script>		
			<table width="100%">
			<tr><td class="itunes-up silver">
				<img src="img/icons/batchupload.png"/>
				<h1>Batch upload 3.0 - error</h1>
				An unexpected rrror has occurred during your batch upload
			</td></tr>
			<tr><td class="itunes-right big-padding ui-widget">
				<div class="message">
					<p>
						Unfortunately, an <b>error</b> occurred during your batch upload.<br/>
						The error may be temporary, so you can try to <a href="/batchupload30/init.do">reattempt your upload attempt</a>.<br/>
						Try to use specify explicitly names of all keywords to avoid any mappings by OCHEM. Delete redundant columns (shown in red color), redundant sheets in the file, etc. 
					</p>
					<p>
						If the problem persists, please feel free to report the case to <a href="http://ochem.eu/support/" target="_blank">our bug tracking system</a>.
						<br/>We are dedicated to constantly improve OCHEM and we appreciate any feedback!
					</p>
					
					<p><br/>Thank you!</p>
					<p>Sincerely,<br/> the OCHEM team</p>
					
					<br/><br style="clear: both"/>P.S.: If you want to report a bug, please include the following detailed information about the failure <a id="toggle">[&gt;&gt;]</a><br/>
					<pre id="strace" class="invisible">
<b>Time of the failure:</b><xsl:value-of select="//message/time"/><br/>
<b>Exception stacktrace:</b><br/><xsl:value-of select="//message/message"/>
					</pre>
				</div>
			</td></tr>
			</table>
	</xsl:template>
	
</xsl:stylesheet>
