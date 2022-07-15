<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />
	
	<xsl:template name="content">
		<style type="text/css">
		</style>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<table height="100%" width="100%">
			<tr><td class="itunes-up silver">
				<h1><img src="img/icons/desc-storage.png"></img>Descriptors upload success</h1>
				</td></tr>
			<tr>
				
				<td class="itunes-right">
					Your descriptors have been successfully uploaded to the descriptor storage!<br/>
					Now, you can use them for the development of QSAR models.
					
					<xsl:if test="//param[name='external-mols']">
						<br/><br/>
						<b>It has not escaped our notice that...</b><br/> you have used external identifiers of molecules (EXTERNAL_ID column).<br/>
						If you want to use such molecules for modeling, do not forget to upload these external identifiers to your experimental data. 
					</xsl:if>
					
					<br/><br/>
					<a href="descriptorsstorage/show.do?render-mode=popup" class="fb-button">Back to the descriptor storage overview</a>
					
				</td>
			</tr>
		</table>
		<script language="javascript">
			var a = new Actionable();
			
			a.doSubmit = function()
			{
				if (!$("[name='desc-type']").val())
				{
					window.alert("Please, provide the descriptors type name");
					return false;
				}
				$(".upload form").submit();
			}
			
			a.doUpload = function()
			{
				$(".upload").removeClass("invisible");
			}
			
			a.doDelete = function()
			{
				var block = this.currentBlock;
				a.ajax.send({
					url: 'descriptorsstorage/delete.do',
					data: 'config-id=' + a.currentEntityId,
					success: function()
					{
						block.remove();
					}
				});
			}
			
			a.setPosition = function(element)
			{
				this.currentBlock = $(element).parents('[rec-id]');
				this.currentEntityId = this.currentBlock.attr('rec-id');
			}	
			
			a.doCancelupload = function()
			{
				$(".upload").addClass("invisible");
			}
			
			function submitForm()
			{
				$(".upload form").submit();
			}
		</script>
	</xsl:template>
</xsl:stylesheet>