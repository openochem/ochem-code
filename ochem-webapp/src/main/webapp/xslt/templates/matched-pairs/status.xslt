<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="2.0">
	
	<xsl:include href="../../helper.xslt" />
	
	<xsl:template name="content">
		<title>MMPA module status browser</title>
		<style type="text/css">
			TABLE.values TD {background-color: #F0F0F0; padding: 5px; border: 1px solid white;}
			TABLE.values {border-spacing: 1px; border-collapse: collapse;}
			TR.header {font-weight: bold;}
			
			.commands {margin-left: 30px;}
		</style>
		<table height="100%" width="100%">
			<tr><td class="itunes-up">
				<h1><span class="setcompare">MatchedPairs</span>: Indexing information</h1>
				</td></tr>
			<tr>
				<td class="itunes-right">
				This page shows technical information regarding the MatchedPairs indexing status.<br/><br/>
					Molecules indexing info:<br/>
					<table class="values" cellspacing="1">
						<tr class="header">
							<td>Indexing status</td>
							<td>Number of molecules</td>
						</tr>
						<xsl:for-each select="//molCount/entry">
						<tr>
							<td><xsl:value-of select="key"/></td>
							<td align="right">
								<xsl:value-of select="format-number(value, '#,###')"/> molecules
							</td>
						</tr>
						</xsl:for-each>
					</table>	
					<br/>
					Database information:<br/>
					<table class="values">
						<tr class="header">
							<td>Table</td>
							<td>Rows</td>
							<td>Data size</td>
							<td>Index size</td>
							<td>Total size</td>
						</tr>
						<xsl:for-each select="//tableInfo/entry">
						<tr>
							<td><xsl:value-of select="key"/></td>
							<td align="right"><xsl:value-of select="format-number(value/count, '#,###')"/> rows</td>
							<td align="right"><xsl:value-of select="format-number(value/dataSize, '#,###')"/> MB</td>
							<td align="right"><xsl:value-of select="format-number(value/indexSize, '#,###')"/> MB</td>
							<td align="right"><xsl:value-of select="format-number(value/size, '#,###')"/> MB</td>
						</tr>
						</xsl:for-each>
						<tr>
							<td colspan="4">Totals</td>
							<td align="right"><xsl:value-of select="format-number(//totalSize, '#,###')"/> MB</td>
						</tr>
					</table>		
				</td>
			</tr>
		</table>
		
		<br/>
		<div class="commands">
			<a class="fb-button" href="matchedpairs/status.do?recalculateCounts=1&amp;render-mode=popup">Recalculate pair counts</a>
			<a class="fb-button" href="matchedpairs/status.do?deleteInfrequent=5&amp;render-mode=popup">Delete infrequent transformations</a>
			<a class="fb-button" href="matchedpairs/status.do?deleteDublicatePairs=1&amp;render-mode=popup">Delete dublicate pairs</a>
		</div>
		<br/><br/>
		<div class="commands">
			<a class="fb-button" href="matchedpairs/status.do?resubmit=stuck&amp;render-mode=popup">Resubmit stuck molecules</a>
			<a class="fb-button" href="matchedpairs/status.do?resubmit=failed&amp;render-mode=popup">Resubmit failed molecules</a>
		</div>
		
	</xsl:template>
</xsl:stylesheet>
