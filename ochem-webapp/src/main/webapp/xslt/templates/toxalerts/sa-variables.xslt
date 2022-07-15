<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />
	
	<xsl:template name="content">
		<style type="text/css">
			TABLE.torefactor TD, TABLE.torefactor TH {background-color: #FFD; padding: 10px; text-align: left; border-bottom: 1px solid #AAA;}
			TABLE.torefactor TH {font-style: bold; border-bottom: 1px solid black; background-color: #FFC;}
			.smarts {font-family: Courier; font-weight: bold;}
		</style>
		<table height="100%" width="100%">
			<tr><td class="itunes-up">
				<h1><span class="toxalerts">ToxAlerts</span>: Alert substitution variables</h1>
				The list of variables you can use in your SMARTS patterns
				
				</td></tr>
			<tr>
				<td class="itunes-right">
					ToxAlerts allows to substitute complex sub-patterns in SMARTS with simpler and more interpretable variables. <br/>
					For example, instead of <span class="smarts">[CH1]([F,Cl,Br,I])([F,Cl,Br,I])[CH2][CH3]</span> you could use <span class="smarts">[CH1]([$Hal])([$Hal])[CH2][CH3]</span><br/>
					<h1>Supported substitution variables</h1>
					Currently, the list of available variables is managed by the ToxAlerts moderation team. If you wish to suggest more variables, please feel free to <a tab="Write a message" href="/mail/action.do?render-mode=popup&amp;action=write&amp;receiver=midnighter&amp;subject=ToxAlerts">drop us a message</a>.
				
					<table class="torefactor">
						<tr>
							<th>Variable</th>
							<th>Title</th>
							<th>Description</th>
							<th>Substitution</th>
						</tr>
						<xsl:for-each select="//alert-variable">
							<tr>
								<td><b>$<xsl:value-of select="variableName"/></b></td>
								<td><xsl:value-of select="name"/></td>
								<td><xsl:value-of select="description"/></td>
								<td><xsl:value-of select="substitution"/></td>
							</tr>
						</xsl:for-each>
					</table>
				</td>
			</tr>
		</table>
		</xsl:template>
</xsl:stylesheet>
