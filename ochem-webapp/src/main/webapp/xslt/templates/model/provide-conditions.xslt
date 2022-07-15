<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />

	<xsl:template name="content">
		<script language="javascript" src="js/commons/midnighter.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<table width="100%">
			<tr><td class="itunes-up">
				<h1>Model applier</h1>
			</td></tr>
			<tr><td class="itunes-right">
				<h1>Provide the default conditions (to be used in case if the respective condition is missed)</h1>
				The models you selected require to specify the following conditions:<br/><br/>
				<form action="modelapplier/provideConditions.do" method="post">
				<xsl:for-each select="//external-descriptor">
					<div>
						<input type="hidden" name="condition-id" value="{id}"/>
						<xsl:variable name="pos"><xsl:value-of select="position()+1-1"/></xsl:variable>
						<span><xsl:value-of select="//others/property[$pos+1-1]/@name"/></span>
						<xsl:if test="//others/property[$pos+1-1]/@type = 0">
							<input type="text" name="condition-value" value="{defaultValue}"/>
							<select id="sel-{position()}" name="condition-unit">
								<xsl:for-each select="//others/property[$pos+1-1]/unitCategory/unit">
									<option value="{@id}">
										<xsl:value-of select="@name"/>
									</option>
								</xsl:for-each> 
							</select>
							<script language="javascript">
								$("#sel-<xsl:value-of select="$pos+1-1"/>").val(<xsl:value-of select="unitId"/>);
							</script>
						</xsl:if>
						<xsl:if test="//others/property[$pos+1-1]/@type = 1">
							<input name="condition-unit" value="0" type="hidden"/>
							<select id="sell-{position()}" name="condition-value">
								<xsl:for-each select="//others/property[$pos+1-1]/option">
									<option value="{@id}"><xsl:value-of select="@name"/></option>
								</xsl:for-each>
							</select>
							<script language="javascript">
								$("#sell-<xsl:value-of select="$pos+1-1"/>").val(<xsl:value-of select="defaultValue"/>);
							</script>
						</xsl:if>
					</div>
				</xsl:for-each>
				
					<div class="formsubmit">
						<input type="button" name="submit" value="&lt;&lt;Back" onclick="location.href='model/select.do';"/>
						<input type="submit" id="next" name="next" value="Start calculations&gt;&gt;"/>
					</div>	
				</form>
			</td></tr>
		</table>
	</xsl:template>
</xsl:stylesheet>