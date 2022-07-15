<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	
	<xsl:include href="../../helper.xslt" />
	<xsl:include href="../model/inc/select-compounds.xslt" />
	
	<xsl:template name="content">
	
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<style type="text/css">
			.step {border: 1px solid #CCC; width: 600px; padding: 10px; margin-bottom: 5px; float: left; margin-right: 5px;}
			.step TABLE TD {padding-bottom: 15px;}
			.letter {font-size: 300%; color: #999; font-family: Arial; font-weight: bold; padding-right: 7px; top: 12px; position: relative;}
			.button {}
			.step HR {color: #999;}
		</style>
		<table height="100%" width="100%">
			<tr><td class="itunes-up">
				<img src="img/icons/compare.png"/>
				<h1><span class="setcompare">SetCompare</span>: Select the sets to compare</h1>
				</td></tr>
			<tr>
				<td class="itunes-right">
				The SetCompare utility is experimental. It allows you to compare two sets of molecules based on their structural features.<br/>
				Please, provide the two sets  available options below.<br/><br/>
				<xsl:if test="//message">
					<div class="warning">
						<xsl:value-of select="//message/message"/>
					</div>
				</xsl:if>
				
				
				<form method="post" enctype="multipart/form-data" action="setcomparison/selectSubmit.do">
					<div class="step">
						<span class="letter">1</span>Select the compounds in the <b>first set:</b><br/><hr/><br/>
						<xsl:call-template name="select-compounds">
							<xsl:with-param name="scope">set1</xsl:with-param>
						</xsl:call-template>
					</div>
					<div class="step">
						<span class="letter">2</span>Select the compounds in the <b>second set:</b><br/><hr/><br/>
						<xsl:call-template name="select-compounds">
							<xsl:with-param name="scope">set2</xsl:with-param>
						</xsl:call-template>
					</div>
					<br style="clear: both;"/>
					<input type="submit" name="submit" value="Next &gt;&gt;" class="button"/><br/><br/>
				</form>
				</td>
			</tr>
		</table>
		
		<script language="javascript">
				
		</script>
		
		</xsl:template>
</xsl:stylesheet>
