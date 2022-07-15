<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	<xsl:include href="inc/license.xslt" />
	
	<xsl:template name="content">
		<style type="text/css">
		 .comp-block-login
		{
			border: outset 1px;
		  	border-spacing: 30pt 10pt;
			background-color:#FFF;
			margin: 10px 10px 10px 10px;
			position: relative;
			width: 770px;
			align: center;
		}
		.comp-block-login INPUT {border: 1px solid #555; width: 200px;}
		.comp-block-login TD {padding: 5px 10px 5px 10px;}
		.comp-block-login TEXTAREA {}
	</style>
	
	<table class="comp-block-login" align="center">
		<tr>
				<td align="center">
					<b>Terms of Service:</b>
				</td>
		</tr>
		<tr>
			<td align="justify">
							<textarea rows="5" cols="80" style="width: 100%; text-align: justify;" readonly="readonly" onfocus="this.rows=10">
							<xsl:call-template name="license"/>
							</textarea>
							<br/><br/><b>By clicking on 'I accept' below I acknowledge that I have read and fully understand the foregoing information and agree to abide by 
							<a href="license_agreement.htm" tab="License agreement">License agreement</a> above and the 
							<a href="Privacy_Policy.htm" tab="Privacy policy">Privacy Policy</a>.</b>
					<br/><br/>
					<div align="center">
		<a href="login/acceptLicenseAgreement.do?render-mode=redirect&amp;accepted=1" class="button-link" target="_parent">I accept</a> <a href="login/acceptLicenseAgreement.do?render-mode=redirect&amp;denied=1" class="button-link">I decline</a>
	</div>
			</td>	
		</tr>
		</table>
	<xsl:if test="not(/model/session/user)">
	
	</xsl:if>
	
	
	
	
						
		</xsl:template>
</xsl:stylesheet>
