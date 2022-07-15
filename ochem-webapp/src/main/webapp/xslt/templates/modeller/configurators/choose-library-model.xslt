<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<script language="javascript" src="actionable.js"></script>
		<script language="javascript" src="ajax-form.js"></script>
		<p>Please, select the base model. The new model will use the base model and enrich its performance using the new compounds from the training set.</p><br/>
		The base model:
		<a action="selectmodel" bindto="basemodel-id" storein="basemodel-title" title="Click to change">[...]</a>
		<div class="invisible" id="multi-property">
			Use predictions for the following property:
			<select name="property-num">
				
			</select>
		</div>		
		<input type="hidden" name="basemodel-title"/>
		<input type="hidden" name="basemodel-id"/>
		<script language="javascript">
		
			var form = new EditForm();
			form.doSelectmodel = function(link)
			{
				var targetLink = link;
				var win = openTab("Select the base model", webRoot+"model/select.do?render-mode=popup&amp;template=ASNN&amp;single=true");
				win.callback = function(entity)
				{
					$("#multi-property select").html();
					form.setValue("basemodel-id", entity.id, entity.name);
					if (entity.modelMappings.length)
					{
						$("#multi-property").removeClass("invisible");
						for (var i = 0; i &lt;  entity.modelMappings.length; i++)
							$("#multi-property select").append($("<option/>").attr("value", i).html(entity.modelMappings[i].property.name));
					}
					updateVisibility();
					win.closeTab();
				}
			}
			
			function updateVisibility()
			{
				if (form.getValue("basemodel-id"))
					$("input[name=next]").removeAttr("disabled");
				else
					$("input[name=next]").attr("disabled", "1");
			}
			
			$(document).ready(function(){
				form.initialize();	
				updateVisibility();
			});
			
		</script>
		
	</xsl:template>
</xsl:stylesheet>