<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />

	<xsl:template name="content">
		<style type="text/css">
			.itunes-right H2 {font-weight: bold; font-family: Arial; margin-top: 20px;}
			.itunes-right H1 { font-family: Georgia; border-bottom: 1px solid; color: #666; }
			
			.itunes-right TABLE TD {padding-right: 20px; padding-bottom: 10px; vertical-align: top;}
			SPAN.author {padding-right: 10px;}
			
			.models {margin-left: 10px;}
			TABLE.models TD {padding-bottom: 5px;}
			.models IMG {margin-right: 2px;}
			A IMG {margin-right: 2px;}
			
			.attach {
				border: 1px solid gray;
				padding: 10px;
				background: #FFE;
			}
			
		</style>
		<title>Article profile</title>
		<table height="100%" width="100%">
			<tr><td class="itunes-up silver">
				<h1><img src="img/icons/document2.png"></img>Article profile</h1>
				A brief overview of the article details
				</td>
			</tr>
			<tr>
				<td class="itunes-right">
				<h1>General information</h1>
					<table>
						<tr><td>Title: </td><td><b><xsl:value-of select="/model/article/@title"/></b>
						
						<xsl:if test="/model/article/@pdf-available = 'true'">
							<br/><a href="pdf/show.do?id={article/@id}&amp;type=pdf"><img src="img/icons/pdf16.gif"/>PDF file available</a>						 
						</xsl:if>
						<xsl:if test="string-length(/model/article/url) &gt; 0">
							<br/><a href="{/model/article/url}" target="_blank"><img src="img/icons/weblink16.gif"/>External web link</a>						
						</xsl:if>
						
						</td></tr>
						<tr><td>Authors: </td><td>
							<xsl:for-each select="//article/authors/author">
								<span class="author"><xsl:value-of select="@printed-name"/>;</span>
							</xsl:for-each>
							</td></tr>
						<tr><td>Journal reference: </td><td><xsl:value-of select="//journal/title"/>, <b><xsl:value-of select="//publication-date/year"/></b>; 
							<xsl:value-of select="/model/article/volume"/> (<xsl:value-of select="/model/article/issue"/>);
							<xsl:value-of select="/model/article/pageNumbers"/>
							</td></tr>
						<xsl:if test="/model/article/pmid">
						<tr><td>PubMed reference: </td><td><a href="http://www.ncbi.nlm.nih.gov/pubmed/{/model/article/pmid}" target="_blank"><xsl:value-of select="/model/article/pmid"/></a></td></tr>
						</xsl:if>
						<tr><td>Internal identifier: </td><td name="article_id">A<xsl:value-of select="/model/article/@id"/></td></tr>
					</table>

				<h1>Data and models</h1>	
					<xsl:if test="article/@property-record &gt; 0">
					This article is referenced from
						<xsl:if test="article/@property-record-nondummy &gt; 0">
							<a href="epbrowser/show.do?article={article/@id}&amp;approval-status=all&amp;hidedummy=true" tab="Article compounds"><xsl:value-of select="article/@property-record-nondummy"/> experimental records</a> 
						</xsl:if>
						
						<xsl:if test="(article/@property-record-nondummy &gt; 0) and (article/@property-record-dummy &gt; 0)">
						 and 
						</xsl:if>
						 
						<xsl:if test="article/@property-record-dummy &gt; 0">
							<a href="epbrowser/show.do?article={article/@id}&amp;approval-status=all&amp;dummy=true" tab="Article compounds"><xsl:value-of select="article/@property-record-dummy"/> molecules without property</a> 
						</xsl:if>
					</xsl:if>

					<xsl:if test="//article/@structural-alerts &gt; 0">
					This article is referenced from <a href="alerts/show.do?article={article/@id}&amp;approval-status=all" tab="Article alerts"><xsl:value-of select="article/@structural-alerts"/> structural alerts</a>
					<br/><br/>
					</xsl:if>
					
					<xsl:if test="//article/mmp-set">
					<h2>Molecular matched pairs (MMPs)</h2>
					This article is associated with MMP studies.<br/>
					Below is the list of the connected MMP transformation sets:<br/><br/>
					<table class="models">
						<xsl:for-each select="//article/mmp-set">
							<tr>
								<td><xsl:value-of select="name"/></td>
							</tr>
							<xsl:for-each select="properties">
							<tr>
								<td>&#160;for <i><xsl:value-of select="@name"/></i></td>
								<td><xsl:value-of select="transformationsCount[1] + transformationsCount[2]"/> transformations</td>
								<td><a href="matchedpairs/annotatedTransformations.do?set={../@id}&amp;property={@id}" tab="MMP transformations for {@name}">view transformations</a></td>
							</tr>
							</xsl:for-each>
						</xsl:for-each>
						</table>
					</xsl:if>
										
					<xsl:if test="//article/@tasks-count &gt; 0">
					<br/><br/>
					This article is related to <a href="pendingtasks/published.do?article-id={article/@id}" tab="Article-related tasks"><xsl:value-of select="article/@tasks-count"/> calculation tasks </a>:
						<table class="models">
						<xsl:for-each select="//article/pending-task">
							<tr><td><xsl:value-of select="type"/> task</td>
							<td><xsl:value-of select="name"/></td>
							<td><a href="pendingtasks/profile.do?id={id}" tab="Task profile"><img src="img/icons/edit.gif"/>view task profile</a></td></tr>
						</xsl:for-each>
						</table>
					</xsl:if>	
					
					<xsl:if test="//article/@model-list &gt; 0">
							This article is connected to <b><xsl:value-of select="article/@model-list"/> predictive model(s)</b>:
							<div class="models">
							<table>
							<xsl:for-each select="//article/models">
									<tr>
										<td valign="top"><a href="model/profile.do?public_id={publicId}" tab="Model profile" title="Open model profile"><xsl:value-of select="@name"/></a></td>
										<td>
											trained using the dataset <a href="basket/edit.do?id={training-set/@id}" tab="Training set profile"><xsl:value-of select="training-set/@name"/></a> 
											(<a href="basket/edit.do?id={training-set/@id}" tab="Training set profile"><img src="img/icons/edit.gif"/>view dataset profile</a> or <a href="basket/exportBasket.do?id={training-set/@id}" tab="Export the dataset"><img src="img/icons/xls.gif"/>export the dataset</a>)
											<xsl:for-each select="validation-sets/validation-set">
											<br/>validated using dataset <a href="basket/edit.do?id={@id}" tab="Validation set profile"><xsl:value-of select="@name"/></a> 
											(<a href="basket/edit.do?id={@id}" tab="Validation set profile"><img src="img/icons/edit.gif"/>view dataset profile</a> or <a href="basket/exportBasket.do?id={@id}" tab="Export the dataset"><img src="img/icons/xls.gif"/>export the dataset</a>)
											</xsl:for-each>
										</td>
									</tr>
							</xsl:for-each>
							</table>
							</div> 
							<br/>						
					</xsl:if>
					
					
					
					<xsl:if test="//session/user">
					<br/>
						<a class="fb-button" href="article/edit.do?id={//article/@id}" tab="Article editor">Open article editor</a>
						<xsl:if test="//session/user and //session/user/rank &gt;= 1">
						<a class="fb-button" onclick="$('.attach').removeClass('invisible'); $(this).addClass('invisible'); return false;" href="#">Attach files to this article</a>
						<br/><br/>
						<div class="invisible attach">
						It is only possible to attach PDF or Excel files, one file type per article.<br/><br/>
						<form action="article/managePdfs.do" enctype="multipart/form-data" method="post">
							<input type="hidden" name="id" value="{//article/@id}"/>
							<input type="hidden" name="action" value="attach"/>
							<input type="file" name="attachment"/>
							<input type="submit" value="Attach the file"/>
						</form>
						</div>
					</xsl:if>
					</xsl:if>
				</td>
			</tr>
		</table>
		
	</xsl:template>
</xsl:stylesheet>