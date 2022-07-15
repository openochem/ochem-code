<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />	
	<xsl:template name="content">
		<script language="javascript" src="js/blocks/exim.js"></script>
		<script language="javascript">	
		$(document).ready(function(){
			$("form").submit(onFormSubmit);
			$("#iframe").load(onFrameLoad);
			$("#basketlink").click(onBasketSelect);
		});
	</script>

	<style type="text/css">
		body iframe {border: 0px;}
		.error {color: #771111;};
	</style>
	
	<title>Data structure import / export</title>
	<table width="100%">
		<tr><td class="itunes-up">
			<h1>Data structure import / export</h1>
			<p>Here you can import or export a property/unit/article structure for your basket</p>
			<iframe name="iframe" id="iframe" width="1" height="1"></iframe>
		</td></tr>
		<tr><td class="itunes-right big-padding">
			<table>
			<tr><td>
				Export:
			</td><td>	
				<form enctype="multipart/form-data" name="exportform" action="structureexim/prepfile.do?out=json" method="post" target="iframe">
					<a id="basketlink" href="javasctipt:void(0)">[...]</a>
					<input type="hidden" name="basketid" value="-1"/>
					<input type="submit" name="submit" disabled="disabled" value="Export"/>
				</form>
			</td></tr>
			<tr><td>
				Import:
			</td><td>	
				<form enctype="multipart/form-data" name="importform" action="structureexim/putfile.do?out=json" method="post" target="iframe">
					<input type="file" name="file" size="25" />
					<input type="submit" name="submit" value="Import"/>
				</form>
			</td></tr>
			<tr><td colspan="2"><b id="status"></b></td></tr>
			<tr><td colspan="2"><div id="log"></div></td></tr>
			</table>
		</td></tr>		
	</table>
	</xsl:template>
</xsl:stylesheet>