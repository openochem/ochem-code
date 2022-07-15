<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../../helper.xslt" />
	<xsl:template name="content">
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		
		<table width="100%">
			<tr><td class="itunes-up">
				<h1><div class="uploadable"><doc term="Upload+a+QSAR+model">Upload a model</doc></div></h1>
				<h1><div class="not-uploadable">Model creator</div></h1>
				<p id="uploadtext">Loading...</p>				

			</td></tr>
			<tr><td class="itunes-right">
					<form action="modelconfigurator/submit.do" id="wizardform" method="post">
					<xsl:if test="//param[@key='multipart'] = 'true'">
						<xsl:attribute name="enctype">multipart/form-data</xsl:attribute>
					</xsl:if>			
					<xsl:call-template name="configuration-content"/>
					<div class="formsubmit">
						<input type="hidden" name="page" value="{//param[@key='page']}"/>
						<xsl:choose>
							<xsl:when test="//param[@key='page'] = 'save'">
								<input type="submit" name="next" value="Save"/>
								<input type="submit" name="discard" value="Discard" onclick="return QSPR.discard = true;"/>
							</xsl:when>
							<xsl:when test="//param[@key='page'] = 'start'">
								<input type="button" name="prev" value="&lt;&lt;Back" onclick="history.go(-1);"/>
								<input type="submit" name="next" value="Start calculation&gt;&gt;" style="background-color:#F99"/>
								<input type="submit" name="discard" value="Discard" onclick="return QSPR.discard = true;"/>
							</xsl:when>
							<xsl:otherwise>
								<input type="button" name="prev" value="&lt;&lt;Back" onclick="history.go(-1);"/>
								<input type="submit" name="next" value="Next&gt;&gt;"/>
							</xsl:otherwise>
						</xsl:choose>
					</div>	
				</form>
			</td></tr>
		</table>
		<script language="javascript">			
			
				if (getParams["upload"])
				{
					$("#uploadtext").html('Select model template, training set and descriptor-coefficient-sheet');			
					$(".not-uploadable").remove(); 
				}
				else
				{
					$("#uploadtext").html('Select model template and training set');			
					$(".uploadable").remove(); 
				}
		</script>
	</xsl:template>
</xsl:stylesheet>