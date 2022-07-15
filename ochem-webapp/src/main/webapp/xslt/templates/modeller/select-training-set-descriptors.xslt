<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />
	<xsl:template name="content">
	<style type="text/css">
		div.validation {margin-top: 5px; margin-left: 5px;}
	</style>
	<title>Model Builder</title>
		<table width="100%">
			<tr><td class="itunes-up">
				<h1>Calculation of descriptors</h1>
				<p>
				Select the set to calculate molecular descriptors for
				</p>
			</td></tr>
			<tr><td class="itunes-right">
				<form name="modelform" id="modelform" action="modelconfigurator/choosesubmit.do" method="post" enctype="multipart/form-data">
					Select the set:
					<a action="selectset" bindto="trainingsetid" storein="tset-title" title="Click to change">[...]</a>
					<span id="trainingsetid-more" class="invisible">   <a action="moreonset" set="trainingsetid">[details]</a></span><br/>
					<br/>	
					<div class="formsubmit"> 
					<input type="button" action="submit" value="Next&gt;&gt;"/>
					</div>
					<input type="hidden" name="trainingsetid"/>
					<input type="hidden" name="tset-title"/>
					<input type="hidden" name="descriptors" value="1"/>
				</form>
			</td></tr>
		</table>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script language="javascript">
			var form = new AjaxForm();
			var self = form;
			
			form.doSelectset = function(link)
			{
				var targetLink = link;
				var baskWin = openTab("Select compound set", webRoot+"basket/show.do?render-mode=popup");
				baskWin.callback = function(basket)
				{
					form.onSetSelected(targetLink, basket);
					baskWin.closeTab();
				}
			}
			
			form.onSetSelected = function(targetLink, basket)
			{
				var properties;
				self.setValue($(targetLink).attr('bindto'), basket.id, basket.name);
				$("span#" + $(targetLink).attr('bindto') + "-more").removeClass("invisible");
				$(document).trigger("DOM_updated", $(document));
			}
			
			var setSetById = function(id)
			{
				$.ajax({
					url: 'basket/list.do?id='+id+'&amp;out=json&amp;public=1',
					dataType: "json",
					success: function(response)
					{
						form.onSetSelected($("a[action=selectset]").get(0), response.list.basket);
					}
				});	
			}
			
			if (getParams["basketId"])
				setSetById(getParams["basketId"]);
			
			form.doSubmit = function()
			{
				var ts = self.getValue("trainingsetid");
				if (ts == undefined || ts == "") 
					return window.alert("Please select a valid trainingset for the model");
				document.modelform.submit();
			}
		</script>
	</xsl:template>
</xsl:stylesheet>