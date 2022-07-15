<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />
	
	<xsl:template name="content">
		<style type="text/css">
			.upload {padding: 20px; background-color: #FFE; border: 1px solid #AAA; margin-top: 10px;}
			.upload TABLE TD {padding: 10px 20px 10px 0px;}
			.top IMG, .upload IMG {margin-right: 3px;}
			
			.descriptors-summary TD, .descriptors-summary TH {padding: 5px;; border-bottom: 1px solid white; background-color: #FAFAFA;}
			.descriptors-summary TH {font-weight: bold;}
		</style>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<table height="100%" width="100%">
			<tr><td class="itunes-up silver">
				<h1><img src="img/icons/desc-storage.png"></img>Descriptors storage overview (under construction!)</h1>
				An overview of the cached and externally uploaded descriptor values
				</td></tr>
			<tr>
				
				<td class="itunes-right">
					<div class="top">
						<a href="#" action="upload"><img src="img/icons/xls.gif"/>Upload precalculated descriptors</a>
					</div>
					<div class="upload invisible">
						<form action="descriptorsstorage/uploadSubmit.do" method="POST" enctype="multipart/form-data">
							OCHEM allows you to use your own precalculated descriptor values for the development of QSAR models.<br/>
							You can upload descriptors as an Excel or csv file with columns "MOLECULE (or EXTERNALID for molecules without structures), DESC_NAME1, DESC_NAME2, ...", where MOLECULE corresponds to OCHEM identifier of the molecule (e.g. M12345) or SMILES string.
							N.B.!: do not use identifiers starting with "M" as EXTERNALID. They will be erroneously considered as OCHEM identifiers.
							<br/>
							A couple of examples with and without molecules: <a href="documents/custom-descriptors-sample.xls"><img src="img/icons/xls.gif"/>Exemplary file 1</a>
							<a href="documents/custom-descriptors-sample-no-mols.xls"><img src="img/icons/xls.gif"/>Exemplary file 2</a>
							<table>
								<tr><td>Descriptor type name<br/><small>e.g., Dragon6, MyCustomDescriptors, E-State</small></td><td><input type="text" name="desc-type"/></td></tr>
<!--								<tr><td>Descriptor configuration XML<br/><small>can be empty for your custom descriptors</small></td><td><textarea type="text" name="desc-conf-xml"/></td></tr>
-->								<tr><td>Excel file with descriptors<br/><small>first line lists the descriptor names as specified above</small></td><td><input type="file" name="desc-file"/></td></tr>
							</table> 
							<br/>
							<a class="fb-button" href="#" action="submit" style="color: white;">Start the upload</a>
							<a class="fb-button" href="#" action="cancelupload">Cancel</a>
						</form>
					</div>
					<br/>
					<xsl:choose>
						<xsl:when test="//desc-config-entry[user]">
							Privately stored descriptors
							<table class="descriptors-summary">
							<tr>
								<th>Descriptor type</th>
								<th>User</th>
								<th>Entries</th>
								<th>Storage size</th>
								<th>Operations</th>
							</tr>
							<xsl:for-each select="//desc-config-entry[user]">
								<tr rec-id="{objectID}" user-id ="{user}">
									<td class="desc-type"><span><xsl:value-of select="type"/></span>
										<div class="invisible">
											<xsl:value-of select="description"/>
										</div>
									</td>
									<td align="right"><xsl:value-of select="user"/></td>
									<td align="right"><xsl:value-of select="entriesCount"/></td>
									<td align="right"><xsl:value-of select="entriesSize"/> bytes</td>
									<td align="right"><a action="delete"><img src="img/icons/delete.gif"/></a>
										<a tab="Export precalculated descriptors" href="descriptorsstorage/export.do?config-id={objectID}" title="Export the descriptors as a file"><img src="img/icons/xls.gif"/></a>
									</td>
								</tr>
							</xsl:for-each>
							</table>
						</xsl:when>
						<xsl:otherwise>
							There no descriptors in the private storage. Try uploading some.
						</xsl:otherwise>
					</xsl:choose>
					<br/><br/>
					
					<xsl:if test="//desc-config-entry[not(user)]">
					Descriptors stored in the public cache
					<table class="descriptors-summary">
							<tr>
								<th>Descriptor type</th>
								<th>Entries</th>
								<th>Storage size</th>
								<th>Operations</th>
							</tr>
							<xsl:for-each select="//desc-config-entry[not(user)]">
								<tr rec-id="{objectID}" user-id="">
									<td class="desc-type"><span><xsl:value-of select="type"/></span>
										<div class="invisible">
											<xsl:value-of select="description"/>
										</div>
									</td>
									<td align="right"><xsl:value-of select="entriesCount"/></td>
									<td align="right"><xsl:value-of select="entriesSize"/> bytes</td>
									<td align="right">
										<xsl:if test="/model/session/user/rank &gt;= 10">
										<a action="delete"><img src="img/icons/delete.gif"/></a>
										<a tab="Export precalculated descriptors" href="descriptorsstorage/export.do?config-id={objectID}&amp;public=1" title="Export the descriptors as a file"><img src="img/icons/xls.gif"/></a>
										</xsl:if>
									</td>
								</tr>
							</xsl:for-each>
							</table>
						</xsl:if>
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
			
			a.doDelete = function(link)
			{
				var block = this.currentBlock;
				a.ajax.send({
					url: 'descriptorsstorage/delete.do',
					data: 'config-id=' + a.currentEntityId + "&amp;user=" + a.user,
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
				this.user = this.currentBlock.attr('user-id');
			}	
			
			a.doCancelupload = function()
			{
				$(".upload").addClass("invisible");
			}
			
			function submitForm()
			{
				$(".upload form").submit();
			}
			
			$(".desc-type").each(function(){
				$(this).find("span").attr("title", "<b>Configuration:</b>" + "<br/>" +  $(this).find("div").html());
			});
		</script>
	</xsl:template>
</xsl:stylesheet>