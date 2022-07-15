<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:template name="recalculate">
		<script language="javascript" src="js/commons/actionable.js"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<xsl:if test="/model/model/@id">
			<div class="popup-footer">
				<xsl:if test="/model/model/@id and /model/model/owner = 'true'">
					<xsl:if test="model/published='false' and model/session/user">
						<a href="javascript:publish()">Publish the model</a>
					</xsl:if>
					<xsl:if test="(/model/model/template/@name != 'Consensus') and (/model/model/template/@name != 'Uploaded model')">
						<a href="javascript:start()">Recalculate model and statistics</a>
					</xsl:if>
				</xsl:if>
				<xsl:if test="/model/model/template/@name != 'Uploaded model'">
					<a href="modelapplier/apply.do?model={/model/model/@id}" tab="Apply a model">Apply the model to new compounds</a><br/>
				</xsl:if>
				<xsl:if test="/model/model/template/@name = 'Uploaded model'">
					<div>This is an uploaded model, which can not be recalculated or applied to new compounds</div>
				</xsl:if>
			</div>
			<div id="progress" style='font-size: 9pt;'><span id="status"></span></div>
			<div id="publish" class="yellow invisible">
				You are about to <b>publish your model</b>. Please, only publish <b>final models</b> that are formally accepted in a journal and, if possible, <b>only one</b> per basket (e.g., consensus).
				Availability of multiple model will confuse the users. Also, you can always share unpublished model by providing its Temporal Public ID (see at the top). 
				Users (e.g., reviewers, your colleagues) who received Temporal link can access the unpublished model and its data. <br></br>
				Once the model is published, it will be visible to all the users together with it's training and validation sets. 
				In the article provide a link to the respective OCHEM article id https://ochem.eu/article/xxxx and not to individual models.
				To complete publication, please <a href="j">select article</a> where the model has been published.
				<br/>
				<div class="popup-footer">
				<a href="javascript:selectArticle()">Select article</a><a href="javascript:cancel()">Cancel</a>
				</div>
			</div>
		</xsl:if>
				
		<script language="javascript">
			var ajax = new QSPR.Ajax();
			var timer = 0;
			
			function publish()
			{
				$("#publish").removeClass("invisible");
			}
			
			function cancel()
			{
				$("#publish").addClass("invisible");
			}
			
			function selectArticle()
			{
				
				var win = openTab("Select article to publish the model", webRoot+"article/show.do?render-mode=popup");
				win.callback = function(article)
				{
					location.href = 'model/action.do?id=<xsl:value-of select="model/@id"/>&amp;action=publish&amp;article=' + article.id;
					win.closeTab();
				}
			}
			
			function checkStatus()
			{
				ajax.url = 'pendingtasks/recalculatestatus.do';
				ajax.send({
					data:'',
					success: function(json)
					{
						var status = json.message.message;
						$("#status").html(status);
						if (status == "Finished")
						{
							window.location.reload(false);	
							clearInterval(timer);					
						}
						if (status.indexOf("Error") != -1)
							clearInterval(timer);
					},
					error: function(msg)
					{
						$("#status").html('QSPR server is not available');
					}
				});
			}
			
			function start()
			{
				$("#status").html('Initializing ...');
				ajax.url = 'pendingtasks/recalculate.do?id=<xsl:value-of select="model/@id"/>';
				ajax.send({
					success: function()
					{
						$("#status").html('Requesting status...');
						timer = setInterval("checkStatus()", 3000);
					}
				});
			}
			
			function rename()
			{
			
				var span = $("#model-name");
				var newName = prompt("Enter the model name", span.html());
				if (newName != null &amp;&amp; newName != '')
				{
					ajax.url = 'model/action.do?id=<xsl:value-of select="model/@id"/>&amp;action=rename&amp;name='+URLEncode(newName);
					ajax.send(
					{
						success: function()
						{
							span.html(newName);
						}
					});
				}
			}
		</script>
		
		<script language="javascript">
		 var ModelApprove = function() {
		 	var ajax = new QSPR.Ajax();
		 	
		 	var form = new AjaxForm();
		 	
		 	$(function(){
		 		$("#approve-dialog").dialog({
					height: 500,
					width: 700,
					modal: true,
					autoOpen: false,
					buttons: {
						"Approve the model": function() {
							ajax.send({
								url: "model/approve.do",
								data: "id=" + <xsl:value-of select="model/@id"/> + form.getActionQuery(),
								success: function(response) {
									window.alert(response.message.message);
									window.location.reload();
								}
							});
						},
						"Cancel": function() {
							$( this ).dialog( "close" );
						0}
					}
				});
		 	});
		 	
			this.dialog = function() {
				$("#approve-dialog").dialog("open");		
			}
		};
		
		var modelApprove = new ModelApprove();
			
		</script>
		<div id="approve-dialog" title="Model approval">
			<p><small>Before approving the model, please make sure that the model is functional (try <a href="modelapplier/apply.do?model={model/@id}" tab="Apply the model">applying this model</a>).
			Please, rate the quality of the model quality considering the predictive performance, the publication status and the value of this model for the community.</small><br/><br/>
			<input type="checkbox" name="publishedAndCited" send="1"/> The model has been published and cited in a scientific publication<br/><br/>
			Model quality assessment:
			<select name="qualityGrade" send="1">
				<option value="3">Excellent</option>
				<option value="2">Good</option>
				<option value="1">Acceptable</option>
			</select></p><br/><br/>
		</div>
	</xsl:template>
	
</xsl:stylesheet>