<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />

	<xsl:template name="content">
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript">
			var form = new AjaxForm();
			var self = form;
			
			form.doEditmolecule = function()
			{
				var molWin = openTab("Edit molecule", webRoot+"molecule/show.do?id="+self.getValue('n-molecule'));
				molWin.callback = function(newId)
				{
					$("[value=molecule]").attr("checked", "checked");
					self.setValue("n-molecule",newId,newId);
					$('#depiction').attr('src', 'depiction.jsp?id='+newId);
					molWin.closeTab();
				}
			}
					
		</script>
		<title>Record editor</title>
		<style type="text/css">
			INPUT {border: 1px solid black; padding: 2px 2px 2px 2px;}
			.options TD {padding-right: 5px; padding-bottom: 20px;}
		</style>
		<table width="100%">
			<tr><td class="itunes-up">
				<h1>Apply the model</h1>
			</td></tr>
			<tr><td class="itunes-right">
				<h1>Provide the secondary compound and molar fraction for the mixture</h1>
				<form name="mixturedata" id="mixturedata" action="modelapplier/provideMixture.do" method="post">
				<div>
				<input type="hidden" name="n-molecule" send="1" value="23427"/>
				<table class="options">
				<tr>			
						<td>Draw secondary compound<br/><i>(click on depiction to the right to draw)</i></td>
						<td><a action="editmolecule"><img id="depiction" style="border: 1px solid gray;" src="depiction.jsp?id=23427" width="100" height="100"/></a></td>
				</tr>
				<tr>
						<td colspan="2">Molar fraction of the primary compound: <input type="text" name="n-molfrac" value="0.9"/></td>
				</tr>				
				</table>				
				</div>
				<div class="formsubmit">
					<input type="button" name="back" value="&lt;&lt;Back" onclick="location.href='model/select.do';"/>
					<input type="submit" value="Next&gt;&gt;"/>
				</div>
				</form>			
			</td></tr>
		</table>
	</xsl:template>
</xsl:stylesheet>