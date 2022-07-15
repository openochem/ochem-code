<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />

	<xsl:template name="content">
		<xsl:variable name="entity">unit</xsl:variable>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/browsers/unit-system-browser.js"></script>
		<script language="javascript">
			include.plugins('view');
			var sampleBrowser = new UnitSystemBrowser("unit-system", "unit-system");
			$(document).ready(function() {sampleBrowser.initialize();});
		</script>
		<title>Unit system browser</title>
		<table height="100%" width="100%">
			<tr><td class="itunes-up">
				<xsl:call-template name="area-of-interest"/>
					<img src="img/icons/property.png"/>
					<h1>Systems of units</h1>
				</td></tr>
			<tr>
				<td class="itunes-right">
					<a action="addnewdialog" title="Add new system of units">
						[Create new system]
					</a>
					<div id="addnew" class="invisible yellow">
						Create new system of units<br/>
						Name <input type="text" filter="1" name="newname"/>
						Default unit name <input type="text" filter="1" name="defunitname"/>
						<a action="addnewsystem">[Okay]</a><a action="cancel">[Cancel]</a>
					</div>
					
					<div class="pager-strip">
						<span><b class="showed">none</b> of <b class="total">none</b></span>
						<div id="pager" class="pgr">
						</div>
					</div>
					
					<div id="Browser">
					</div>
					
					<div class="pager-strip">
						<span><b class="showed">none</b> of <b class="total">none</b></span>
						<div id="pager" class="pgr">
						</div>
					</div>
					
				</td>
			</tr>
		</table>
	</xsl:template>
</xsl:stylesheet>