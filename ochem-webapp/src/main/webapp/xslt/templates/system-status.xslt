<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	<xsl:template name="content">
		<title>System status - periodic results</title>
		<style type="text/css">
			tr.header td {font-size: 80%; text-align: center; padding: 5px 10px 5px 10px;}
			table#results tr td {border-bottom: 1px solid black; padding-left: 15px; padding-right: 15px;}
			table#results td.first {padding: 3px 15px 3px 15px; }
			table#results span {height: 25px; width: 5px; margin: 1px 1px 1px 0px; display: block; float: left;}
			span.red {background-color: #E00;}
			span.green {background-color: #5A5;}
			span.gray {background-color: #EEEEF3;}
			
			.dura {
				font-size: 70%; text-align: right; position: relative; top: 5px;
				width: 50px !important;
				padding-left: 5px; padding-right: 3px; color: #777;
			}
			
			table#results tr:hover > TD {background-color: #DED;} 
			
			
			.recently-failed {color: red;}
		</style>
		<table width="100%">
	  		<tr>
	  			<td class="itunes-up">
	  				<h1>System status</h1>
	  				Results of automatic hourly tests
	  			</td>
	  		</tr>
	  		<tr><td class="itunes-right">
	  			<xsl:if test="//param[@key='running']">
	  				<img src="img/roller_small.gif"/> 
	  				<font color="#223">The tests are currently running!
	  				<xsl:if test="//param[@key='running']"> 
	  					(<xsl:value-of select="//param[@key='running']"/>)
	  				</xsl:if>
	  				</font>
	  			</xsl:if>
	  			<xsl:if test="not(//param[@key='running'])">
	  				<a href="systemstatus/forceRun.do?type={model/others/testresult[1]/@testType}" class="fb-button">Force the full test run now</a>
	  			</xsl:if>
	  			<a class="fb-button" tab="System monitors" href="systemstatus/monitors.do">Interactive monitors</a>
				<table id="results" cellspacing="1">
					<tr class="header"><td></td></tr>
				</table>
			</td>
			</tr>
		</table>
		<xsl:for-each select="model/others/testresult">
			<xsl:if test="detailedStatus">
				<div class="invisible" id="details{@id}">
					<xsl:value-of select="detailedStatus"/>	
				</div>
			</xsl:if>
		</xsl:for-each>
		
		<script type="text/javascript" src="js/blocks/system-status.js"/>
		<script language="javascript">	
		
		<xsl:for-each select="model/others/testresult">
			<xsl:if test="@testType='selenium'">
				sysstat.getSeleniumTestArray('<xsl:value-of select="@className"/>', '<xsl:value-of select="@methodName"/>');
			</xsl:if>
		</xsl:for-each>
		
		<xsl:for-each select="model/others/testresult">
		<xsl:choose>
			<xsl:when test="@testType='selenium'">
				sysstat.addSeleniumTest(<xsl:value-of select="@id"/>, '<xsl:value-of select="@className"/>', '<xsl:value-of select="@methodName"/>', <xsl:value-of select="@succeeded"/>, '<xsl:value-of select="@time"/>', '<xsl:value-of select="@duration"/>');
				<!-- sysstat.addSeleniumTestResult('<xsl:value-of select="@id"/>', '<xsl:value-of select="@className"/>', '<xsl:value-of select="@methodName"/>', <xsl:value-of select="@succeeded"/>, '<xsl:value-of select="@time"/>', '<xsl:value-of select="@duration"/>');-->			
			</xsl:when>
			<xsl:otherwise>
				sysstat.addTestResult(<xsl:value-of select="@id"/>, '<xsl:value-of select="@name"/>', <xsl:value-of select="@succeeded"/>, '<xsl:value-of select="@time"/>', '<xsl:value-of select="@duration"/>', '<xsl:value-of select="@methodName"/>');
			</xsl:otherwise>
		</xsl:choose>
		</xsl:for-each>
		
		sysstat.highlightRecentlyFailedTests();
		//table.find("span").tooltip();
		</script>
	</xsl:template>
	
</xsl:stylesheet>
