<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />

	<xsl:template name="content">
		
		<table width="100%" height="100%" cellspacing="0">
			<tr>
				<td class="itunes-up" colspan="2">
					<h1>Basket splitter</h1>
				</td>
			</tr>
			<tr>
				<td class="itunes-right" >
					<form method="post" id="mainform">
						<input type="hidden" name="id" value="{//basket/@id}"/>
						You are going to split the basket <b><xsl:value-of select="//basket/@name"/></b> into two new baskets.
						
						<h3>Provide the basket names</h3>
						Basket 1: <input type="text" value="{//basket/@name} (training)" name="basket-name-1" size="70"/><br/>
						Basket 2: <input type="text" value="{//basket/@name} (test)" name="basket-name-2" size="70"/><br/>
						
						<h3>Select the splitting method</h3>
						<input type="radio" name="splitting-method" value="random" checked="checked"/> Random splitting<br/>
						<div class="params random">
							Size of the Basket 2, in percentages: <input type="text" value="20" name="validation-set-percentage"/> %
						</div>
						<!--
						<input type="radio" name="splitting-method" value="y-based" disabled="disabled"/> Y-based splitting (not implemented yet)
						<br/><br/>
						-->
						Your original basket will be preserved.<br/><br/>
						<input type="submit" value="Split the basket" onclick="start(); return false;"/>
					</form>
				</td>
			</tr>
		</table>
		
		<script language="javascript">
			function start()
			{
				new LongOperation({
					url: "basket/split.do",
					formScope: "#mainform",
					finished: function(status)
					{
						window.alert("The basket has been splitted successfully!");
						window.location.href = webRoot + "/basket/show.do?render-mode=popup";
					}
				}).start();
			}
		</script>
		
	</xsl:template>
</xsl:stylesheet>