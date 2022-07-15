<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xs="http://www.w3.org/2001/XMLSchema" version="2.0">
	<xsl:include href="../helper.xslt" />
	<xsl:template name="content">
	<script type="text/javascript" src="js/commons/actionable.js" />
	<script type="text/javascript" src="js/commons/browser.js" />
		<title>Bonus points overview</title>
		<style>
			.report TD, .report TH {padding-right: 25px; vertical-align: top; }
			.report TH {font-weight: bold;}
			TH.right {text-align: right !important;}
			.transactions TD {padding: 4px 25px;}
			.transactions TR {border-bottom: 1px solid gray; }
			.sum {font-weight: bold; border-top: 1px solid black;}
		</style>
		
		<table height="100%" width="100%">
			<tr>
				<td class="itunes-up" colspan="2">
					<img src="img/icons/award-48.png"/>
					<h1>Bonuses summary</h1>
					Overview of the bonus points
				</td>
			</tr>
			<tr>
				<td class="itunes-right">
					OCHEM is a free platform and most of the data can be exported without restrictions.<br/>
					However, to stimulate give-and-take, export of some sensitive data (e.g., molecular structures) requires bonus points, which can be earned by contributing and publishing new data to OCHEM.
					<a href="https://docs.ochem.eu/display/MAN/Bonus+points+system" target="_blank">Learn more about bonus points.</a>
					<br/><br/>
					<h2>Your bonus points overview</h2>
					<table class="report">
						<tr>
							<td>
								Free bonus points
								<a class="infolink" title="You receive free recoverable bonus points every week. In case they are not suddicient for your demands, you can earn points by contributing high quality data to OCHEM users."></a>
							</td>
							<td class="right"><xsl:value-of select="/model/bonuses-report/freeBonusesAllowance"/></td>
						</tr>
						<tr>
							<td>Spent during the last 7 days<br/><br/></td>
							<td class="right">(<xsl:value-of select="/model/bonuses-report/freeBonusesSpentLastWeek"/>)</td>
						</tr>
						<tr>
							<td>Regular bonus points earned</td>
							<td class="right"><xsl:value-of select="/model/bonuses-report/regularBonusesEarned"/></td>
						</tr>
						<tr>
							<td>Regular bonus points spent<br/><br/></td>
							<td class="right">(<xsl:value-of select="/model/bonuses-report/regularBonusesSpent"/>)</td>
						</tr>
						<tr class="sum">
							<td>Total bonus points balance</td>
							<td class="right"><xsl:value-of select="/model/bonuses-report/totalBonusesSaldo"/></td>
						</tr>
					</table>
					
					<br/><br/>
					<h2>Bonus points earnings/spendings</h2><br/>
					<table class="report transactions">
						<tr>
							<th>Date</th>
							<th>Description</th>
							<th class="right">Decrease</th>
							<th class="right">Increase</th>
						</tr>
						<tbody id="Browser">
						</tbody>
					</table>
				</td>
			</tr>
		</table>
		
		<script language="javascript">
			include.plugins('view');
			
			$(function(){
				var browser = new Browser();
				browser.url = "bonuses/list.do";
				browser.view = new View({element: "transaction-template"});
				browser.useTableRows = true;
				browser.itemElement = "bonusTransaction";
				browser.initialize();
			});
		</script>
		
		<script type="text/template" id="transaction-template">
			<td>[%=data.time %]</td>
			<td>[%=data.description %]</td>
			<td class="right">[% if (data.amount &lt; 0) { %]([%=data.amount %])[% } %]</td>
			<td class="right">[% if (data.amount &gt; 0) { %][%=data.amount %][% } %]</td>
		</script>
		
	</xsl:template>
	
	
</xsl:stylesheet>