<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Model builder - Start calculation</title>
		<style type="text/css">
			.warning {background-color: #FCC;}
		</style>
		<h1>Start calculation of the model</h1>
		Now we are ready to start calculation.<br/>
		Please provide the name for your model:
		<input type="text" name="name" value="{model/@name}" size="100"/><br/><br/>
		<input type="checkbox" name="save-models" checked="checked"/> Save models<br/><br/>

		Task priority:<br/>
		<xsl:if test="session/@max-priority &gt;= 10"> 
			<input type="radio" name="priority" value="10"/> Extra-high priority (please, use for fast tasks only)<br/>
		</xsl:if>
		<input type="radio" name="priority" value="2"></input> High priority (please, use for fast tasks only)<br/>
		<input type="radio" name="priority" value="1" checked="true"/> Normal priority<br/>
		<input type="radio" name="priority" value="0"/> Low priority (for long tasks)<br/>		
		<xsl:if test="(//session/user/@login = 'itetko') or (session/@max-priority &gt;= 10)">
			<br/><br/>Preferred calculation server: <input type="text" name="preferred-server"/> (developers only)
		</xsl:if>
	</xsl:template>
	
</xsl:stylesheet>