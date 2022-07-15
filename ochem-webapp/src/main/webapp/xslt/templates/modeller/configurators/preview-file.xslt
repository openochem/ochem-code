<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Model file preview</title>
		<style type="text/css">
			.preview
			{
						
			}
			
			.previewtable
			{
				border-collapse: separate;
				border-spacing: 1px;
				min-width: 30%;
			}
			
			.previewtd
			{
				padding: 7px 7px 7px 7px;			
				background-color: #CCCCFF;		
			}
			
			.previewth
			{
				padding: 2px 7px 2px 7px;			
				background-color: #AAAAEE;	
				vertical-align: top;	
			}
			
			.green
			{
				background-color: #AAEEAA !important;	
			}	

			.dgreen
			{
				background-color: #66AA66 !important;	
			}	
			
			.red
			{
				background-color: #EEAAAA !important;	
			}			
			
			.button-link
			{
				background-color:#EEEEFF;
				border:1px solid #666688;
				color:#222288;
			}
			
			.scroll
			{
				width: 100%;
				overflow: auto;
			}
				
		</style>
		<h1>Model file preview</h1>
		<br></br>
<!--    <input type="checkbox" name="xls"/> -->
        <table class="previewtable">
        <tr><th class="previewth green">Descriptor</th><th class="previewth green">Value</th></tr>
         <xsl:for-each select="list/list[1]/list">
         	<tr>
         	<xsl:for-each select="string">
         		<td class="previewtd">	<xsl:value-of select="."/></td>
         	</xsl:for-each>
         	</tr>
         </xsl:for-each>
        </table>
		<br></br>
		<table class="previewtable">
			<xsl:for-each select="list/list[2]/list[1]">
				<tr>
					<th class="previewth green">Descriptor</th>
					<th class="previewth green">Value</th>
					<th class="previewth green">Remap Factor</th>
					<th class="previewth green">Remap Constant</th>
					<th class="previewth green">New Descriptor</th>
					<th class="previewth green">New Value</th>
					<th class="previewth green">Offset Shift</th>
				</tr>
			</xsl:for-each>
			<xsl:for-each select="list/list[2]/list">
				<tr>
					<xsl:for-each select="string">
						<td class="previewtd">
							<xsl:value-of select="." />
						</td>
					</xsl:for-each>
				</tr>
			</xsl:for-each>
		</table>
		<br></br>
	</xsl:template>
</xsl:stylesheet>

<!-- Toxicity/xslt/templates/modeller/configurators/upload-file.xslt tee 20091102 -->
