<?xml version="1.0" encoding="UTF-8"?>
<!-- Universal template to track the status of long-running server side operations -->
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	<xsl:template name="content">
		<style type="text/css">
				.summary {width: 800px; font-size: 120%; margin: 10px;}
				.red {background-color:#FEE; border-top: 1px solid black;}
				.lgreen {background-color:#AFA}
				.dgreen {background-color:#8F8}
				.summary TD {padding: 10px;}
		</style>
		<br/><br/>
		<center>
		<table class="summary">
			<tr><td><b>File name:</b></td><td><xsl:value-of select="//export-action/fileName"/>.<xsl:value-of select="//export-action/format"/></td></tr>
		</table><br/>
		<xsl:choose>
			<xsl:when  test="//export-action/fileDownloaded = 'true'">
				
				<table class="summary">
					<tr><td colspan="2">Your file is ready and the download will start in a few seconds. If the download does not start, you can use the direct link below.</td></tr>
					<tr><td><a href="{//export-action/publicURL}" class="fancy-button" style="width: 300px; text-align: center;">Download</a></td><td><a id="close" href="javascript:void(0)" class="fancy-button-red" style="width: 300px; text-align: center;">close this page</a></td></tr>
				</table>
				
				<script language="javascript">
					$(document).ready(function(){
						var url = "<xsl:value-of select="//export-action/publicURL"/>";
						window.location.replace(url);
					});
				</script>
			</xsl:when>
			<xsl:otherwise>
				<table class="summary">
					
					<tr class="lgreen">
					<td width="300px">
						Available free bonus points <a class="infolink" title="You receive free recoverable bonus points every week. In case they are not suddicient for your demands, you can earn points by contributing high quality data to OCHEM users."></a>
					</td><td class="right"><xsl:value-of select="//export-action/freeBonusPointsAvailable"/></td></tr>
					<tr class="lgreen"><td width="300px">Available earned bonus points</td><td class="right"><xsl:value-of select="//export-action/earnedBonusPointsAvailable"/></td></tr>
					<tr class="dgreen"><td width="300px">Available total bonus points</td><td class="right"><b><xsl:value-of select="//export-action/totalBonusPointsAvailable"/></b></td></tr>
					<tr class="red"><td width="300px"> Export cost in bonus points</td><td class="right">(<xsl:value-of select="//export-action/totalBonusPoints"/>)</td></tr>
				</table><br/>
				
				
					<xsl:choose>
				
						<xsl:when  test="//export-action/enoughBonuses = 'true'">
							<table class="summary">
								<tr><td><a href="export/confirm.do?id={//export-action/@id}" class="fancy-button" style="width: 300px; text-align: center;">Spend <xsl:value-of select="//export-action/totalBonusPoints"/> bonus points and download</a></td><td><a id="close" href="javascript:void(0)" class="fancy-button-red" style="width: 300px; text-align: center;">close this page</a></td></tr>
							</table>
						</xsl:when>
						<xsl:otherwise>
							<table class="summary">
								<tr><td>Unfortunately, you can not download this file due to insufficient bonus points. Please, wait for the free bonus points to restore or earn more bonus points.</td></tr>
								<tr><td colspan="2"><a id="close" href="javascript:void(0)" class="fancy-button-red" style="width: 300px; text-align: center;">close this page</a></td></tr>
							</table>
						</xsl:otherwise>
					</xsl:choose>
			</xsl:otherwise>
		</xsl:choose>
		</center>
		<script language="javascript">
			$("#close").click(function(){ 
				window.closeActiveTab(); 
			});
		</script>		
	</xsl:template>
</xsl:stylesheet>
