<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />

	<xsl:template name="content">
		<style type="text/css">
			#area {margin: 10px;}
			.big TD {padding-right: 5px; padding-bottom: 5px;}
			.big INPUT {font-size: 15pt; width: 600px;}	
			#status {padding-left: 30px; font-size: 20pt;}
			#depiction IMG {width: 100px; height: 100px; border: 1px solid gray;}
		</style>
		<title>SMART matcher</title>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<div id="area">
		<h1>Simple SMARTS matcher</h1>
		This is a simple utility to check whether a SMILES molecule matches against a SMARTS pattern.<br/>Please, enter a SMILES of your molecule and the SMARTS pattern (extended syntax and <a href="alerts/variables.do" tab="Substitution variables">substitution variables</a> are supported).<br/><br/>
		
		<table class="big">
			<tr><td></td><td><div id="depiction"></div></td></tr>
			<tr><td>SMILES</td><td><input type="text" send="1" name="smiles" value="CO"/></td><td rowspan="2" id="status"></td></tr>
			<tr><td>Extended SMARTS</td><td><input type="text" name="smarts" send="1" value="NOT [!C;!H]"/></td></tr>
		</table>
		<br/><a action="match" class="button-link">Check matching</a>
		
		</div>
		<script language="javascript">
			var MyActionable = function(){
				AjaxForm.call(this);
				this.actionURL = webRoot + "/alerts/match.do";
				this.onMatchSuccess = function(response)
				{
					var smi = a.getValue('smiles');
					$("#depiction").html('<img src="depiction.jsp?mol='+URLEncode(smi)+'"/>');
					$("#status").html(response.message.message == "1" ? "Matched!" : "Not matched");
				}
			};
			var a = new MyActionable();
		</script>
	</xsl:template>
</xsl:stylesheet>