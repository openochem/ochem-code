<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="2.0">
	<xsl:include href="../helper.xslt" />
	<xsl:include href="inc/data-report.xslt" />
	<xsl:template name="content">
	<script language="javascript" src="js/commons/actionable.js"></script>
	<style>
	h1 {
		color: #333;
		display:block;
		font-family: Helvetica;
		font-size:24px;
		margin-bottom: 10px;
	}
	
	.acknowledgements td {
		padding: 10px;
		width: 700px;
	}
	
	body.popup {
		padding: 40px !important;
		
	}
	
	.user {padding-top: 20px; padding-bottom: 20px;}
	.avatar {margin-right: 20px;}
	.user {border-bottom: 1px solid gray; border-top: 1px solid gray;}
	.user TABLE TD {padding-right: 10px; padding-bottom: 3px;}
	
	TD.contributions {vertical-align: top; padding: 20px; width: 50%; border: 10px solid white; background-color: #FFE;}
	
	.avatar {border: 1px solid #333;}
	
	</style>
	<title>Add New Application</title>
	
	<table height="100%" width="100%">
			<tr>
				<td>
				<h1>Existing Applications</h1>
				</td>
			</tr>
			
			<tr>
				<td class="compact-item">
						<b>Application ID	</b>
				</td>
				<td class="compact-item">
					<b>Application Secret</b>
				</td>
				<td class="compact-item">
					<b>Domain List</b>
				</td>
				<td class="compact-item">
					<b>Remove App</b>
				</td>
			</tr>
			
			
			<xsl:for-each select="/model/others/app">
				<tr>
					<td class="compact-item">
						<xsl:value-of select="clientID"/>
					</td>
					<td class="compact-item">
						<xsl:value-of select="secret"/>
					</td>
					<td class="compact-item">
						<xsl:value-of select="domainList"/>
					</td>
					<td class="compact-item">
						<form method="post" target="_parent" accept-charset="UTF-8,ISO-8859-1">
					<xsl:attribute name="action">oauth/delete.do?clientID=<xsl:value-of select="clientID"/></xsl:attribute>
					<xsl:attribute name="name">deleteform_<xsl:value-of select="@id"/></xsl:attribute>  
					<input type="hidden" name="render-mode" value="redirect"/>
					<input type="hidden" name="create-new" value="false"/>
						<table>
							<tr>
								<td colspan="2" height="30"	valign="bottom">
									<a href="#" class="button-link"><xsl:attribute name="onclick">document.deleteform_<xsl:value-of select="@id"/>.submit(); return false;</xsl:attribute>Delete</a>
									<input type="image" value="submitnew" src="submit-button.gif" width="1" height="1" border="0" alt="SUBMIT!" name="image" class="invisible"/>
								</td>
							</tr>
						</table>
					</form>
					</td>
				</tr>
			</xsl:for-each>

		</table>	
		
		<form method="post" name="mainform" target="_parent" accept-charset="UTF-8,ISO-8859-1">
					<xsl:attribute name="action">oauth/create.do?user=<xsl:value-of select="/model/session/user/@login"/></xsl:attribute> 
					<input type="hidden" name="render-mode" value="redirect"/>
					<input type="hidden" name="create-new" value="true"/>
						<table>
							<tr>
								<td colspan="2" height="30"	valign="bottom">
									List of allowed domains (comma-separated, localhost is always included) <input type="text" name="domain-list" palceholder="localhost is always included"/>
								</td>
								<td colspan="2" height="30"	valign="bottom">
									<a onclick="document.mainform.submit(); return false;" href="#" class="button-link">Create New</a>
									<input type="image" value="submitnew" src="submit-button.gif" width="1" height="1" border="0" alt="SUBMIT!" name="image" class="invisible"/>
								</td>
							</tr>							
						</table>
					</form>
	
	</xsl:template>
	
</xsl:stylesheet>
