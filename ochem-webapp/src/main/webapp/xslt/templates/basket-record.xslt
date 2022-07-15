<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />

	<xsl:template name="content">
		<script language="javascript" src="js/commons/midnighter.js"></script>
		<style type="text/css">
			input[type=radio] {margin-right: 5px;}
			#type-file, #type-basket {margin-left: 20px; background-color: #FFA; padding: 20px;}
		</style>
		<table width="100%">
			<tr><td class="itunes-up">
				<h1>Upload records in Basket</h1>
			</td></tr>
			<tr><td class="itunes-right">
			<form action="basket/addRecords.do" method="post" enctype="multipart/form-data">
				Select the desired action:<br/><br/>
					The following actions affect the basket <b><xsl:value-of select="basket/@name"/></b>:<br/>
					<input type="radio" name="action" value="add" checked="true"/>Add records to the basket<br/>
					<input type="radio" name="action" value="delete"/>Delete records from the basket<br/>
					<br/>
					The following options affect only the models based on the basket <b><xsl:value-of select="basket/@name"/></b>:<br/>
					<input type="radio" name="action" value="exclude"/>Exclude records from the basket <br/>
					<input type="radio" name="action" value="include"/>Include the records to the basket<br/>
					<input type="radio" name="action" value="includeall"/>Include all the previously excluded records to the basket<br/>
				
			<br/><br/>
			
			<div id="records">
			<b>Which records should be processed?</b><br/><br/>
			<input type="radio" name="type" value="file" checked="true"/>Records listed in a file<br/>
			<div id="type-file">
				The file must contain "RECORDID" column, which should list internal OCHEM record identifiers.<br/>
				<label>Excel file: </label><input type="file" name="xls-file"/>
			</div><br/>
			<input type="radio" name="type" value="basket"/>Records from another basket<br/>
			<div id="type-basket">
				Select a basket: <a href="#" onclick="selectBasket(this); return false;">[...]</a>
				<input type="hidden" name="another-basket-id" value="{basket/@id}" class="invisible"/>
			</div>
			<br/>
			</div>
					<br/><br/>
					<input type="hidden" name="id" value="{basket/@id}" class="invisible"/>
					<input type="submit" value="Submit"/>			
			</form>
			</td>
			</tr>
		</table>
		<script language="javascript">
			function updateVisibility()
			{
				$("#records").setClass("invisible", $("input[name='action']:checked").val() == "includeall");
				$("#type-file").setClass("invisible", $("input[name='type']:checked").val() != "file");
				$("#type-basket").setClass("invisible", $("input[name='type']:checked").val() != "basket");
			}
			
			$("input[type=radio]").click(function(){
				updateVisibility();
			});
			
			// Duplication - refactor this kind of functions
			var selectBasket = function(link)
			{
				var targetLink = link;
				var baskWin = openTab("Select compound set", webRoot+"basket/show.do?render-mode=popup");
				baskWin.callback = function(basket)
				{
					$(targetLink).html(basket.name);
					$('input[name="another-basket-id"]').val(basket.id);
					
					baskWin.closeTab();
				}
			}
			
			updateVisibility();
		</script>
	</xsl:template>
</xsl:stylesheet>