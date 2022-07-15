<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />

	<xsl:template name="content">
		<title>Property editor</title>
		<style type="text/css">
			.form TD {padding: 10px 10px; vertical-align: top;}
			.form INPUT[type=text], .form TEXTAREA {width: 500px; background-color: #FFFFEE; border: 1px solid #444;}
			
			.popup-content
			{
				margin: 20px;
			}
			
			.failure-info {
				border: 1px solid gray;
				padding: 10px;
				margin: 5px 0px;
			}
		</style>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		
		<h1 class="popup-header">Task profile
		<i>A profile for a the calculation task</i></h1>
		
		
		
		<div class="popup-content">
		<xsl:if test="pending-task/published = 'true'">
			<img src="img/icons/public.png" title="This task is publicly available."/>This was published by <a tab="User profile" title="Click to show user profile" href="user/profile.do?login={pending-task/session/user/@login}"><xsl:value-of select="pending-task/session/user/@login"/></a> and is available to all OCHEM users<br/><br/>
		</xsl:if>
		<input type="hidden" name="id" value="{pending-task/id}" send="1"/>
		<table class="form">
			<tr>
				<td>Task type</td>
				<td>
					<xsl:value-of select="pending-task/type"/>
				</td>
			</tr>
			<tr>
				<td>Status</td>
				<td>
					<xsl:choose>
						<xsl:when test="pending-task/status = 'ready'">Task has been successfully completed<br/><a tab="Task results" href="pendingtasks/fetchnew.do?id={pending-task/id}">View task results</a></xsl:when>
						<xsl:when test="pending-task/status = 'error'">This task has failed
							<div class="failure-info">
								<xsl:value-of select="pending-task/detailedStatus"/>
							</div>
						</xsl:when>
					</xsl:choose>	
				</td>
			</tr>
			<tr>
				<td>Name</td>
				<td><textarea name="name" send="1" style="height: 2.7em;"><xsl:value-of select="pending-task/name"/></textarea></td>
			</tr>
			<tr>
				<td>Description</td>
				<td>
					<textarea name="description" style="height: 3.7em;" send="1"><xsl:value-of select="pending-task/description"/></textarea>
				</td>
			</tr>
			<tr>
				<td>Associated publication</td>
				<td>
					<a href="" action="article" bindto="article" title="Click to select another publication">
						<xsl:choose>
							<xsl:when test="pending-task/article">
								<xsl:value-of select="pending-task/article/@title"/>
							</xsl:when>
							<xsl:otherwise>
								[...]
							</xsl:otherwise>
						</xsl:choose>
					</a>
					<div class="invisible article-info">
						<a action="article_details">[publication details]</a>
					</div>
					<input type="hidden" name="article"  value="{pending-task/article/@id}" send="1"/>
				</td>
			</tr>
			<xsl:if test="pending-task/model">
				<tr>
					<td>Model</td>
					<td><xsl:value-of select="pending-task/model/@name"/><br/>
						<a href="model/profile.do?id={pending-task/model/@id}" tab="Modle profile">[model profile]</a>	
					</td>
				</tr>
			</xsl:if>
			<tr>
				<td colspan="2"><a action="edit" class="fb-button">Save changes</a></td>
			</tr>
		</table>
		
		</div>
		
		<script language="javascript">
			var form = new EditForm();
			form.actionURL = "pendingtasks/editProfile.do";
			
			form.doArticle = function()
			{
				var win = openTab("Select article to publish the task", webRoot+"article/show.do?render-mode=popup");
				win.callback = function(article)
				{
					form.setValue("article", article.id, article.title);
					win.closeTab();
				}
			}
			
			form.onEditSuccess = form.doCancel = function()
			{
				closeTab();
			}
			
			form.doArticle_details = form.doCancel = function()
			{
				openTab("Article profile", webRoot+"article/profile.do?id=" + form.getValue("article"));
			}
			
			$(document).ready(function(){
				if (form.getValue("article"))
					$(".article-info").removeClass("invisible");
				
			});
			
		</script>
	</xsl:template>
	
</xsl:stylesheet>