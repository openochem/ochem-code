<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />
	<xsl:include href="./inc/select-compounds.xslt" />

	<xsl:template name="content">
		<script language="javascript" src="js/commons/actionable.js?a=1"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		
		<title>Record editor</title>
		<style type="text/css">
			INPUT {border: 1px solid black; padding: 2px 2px 2px 2px;}
			.options TD {padding-right: 5px; padding-bottom: 20px;}
		</style>
		<table width="100%">
			<tr><td class="itunes-up">
				<h1>Model Applier</h1>
			</td></tr>
			<tr><td class="itunes-right">
				<xsl:if test="//message/message">
						<div class="warning">
							<xsl:value-of select="//message/message"/>
						</div>
					</xsl:if>
				<h1>Provide the compound(s) to predict</h1>

				Please provide compounds for which you want to predict the target property<br/>Several options are available:<br/><br/>
				
				<form name="modeldata" id="modeldata" action="modelapplier/apply.do" method="post" enctype="multipart/form-data">
				<div>
					<xsl:call-template name="select-compounds"/>
					<h1>Additional options</h1>
					<input type="checkbox" name="disable-cache"/> Disable prediction cache <br/>
					<xsl:if test="(//session/user/@login = 'itetko') or (session/@max-priority &gt;= 10)">
						<br/><input type="checkbox" id="force_cache" name="force_cache"/><label for="force_cache">Force recalculation of descriptors in cache </label><br/>
						<br/>Preferred calculation server: <input type="text" name="preferred-server"/> (developers only)
					</xsl:if>
				</div>
				<div class="formsubmit">
					<input type="button" action="submit" value="Next&gt;&gt;"/>
				</div>
				</form>			
			</td></tr>
		</table>
	</xsl:template>
</xsl:stylesheet>