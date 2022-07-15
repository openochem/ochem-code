<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	
	<xsl:template name="content">
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script language="javascript" src="js/blocks/article-edit.js"></script>
		<script language="javascript">
		function TestDataCheck() 
		{
			var title = $("[name='n-title']").val();
			var journal = $("[name='n-journal-title']").val();
			var journalId = $("[name='n-journal-key']").val();
			var retVal;
			if (title == "" || (journal == "" &amp;&amp; journalId == "")){
				alert("Title and Journal must not be empty!");
				return false;
			} else {
				return true;
			}
		}
		</script>
		<style type="text/css">
			.outer TD {}
			.inner TD {padding: 2px 2px 2px 2px;}
			TD.inner {border: 1px solid #999; padding: 10px;}
			BODY IFRAME {border: 0px;}
			INPUT, TEXTAREA {border: 1px solid black; padding: 2px 2px 2px 2px; margin: 1px 1px 1px 1px;}
			.outer IMG {border: 1px solid black;}
			.inner TR TD:first-child {width: 100px;}
			.floating {float: left; clear:none;}
			.right{align:right;}
			.hidden {display: none;}
			.message {font-size: 11px; font-style: italic;}
			#journal-info {font-size: 8pt; font-color: #224; background-color: #AAE; padding: 3px;}
		</style>
		<title>Article profile</title>
		<h1 class="popup-header">
		Article/Book editor
		<a help="editor-help" class="infolink"></a>
		<div id="editor-help" class="invisible">
		 	Editor for articles and Books.
		</div>
		<i>Edit and explore information about an article/book</i></h1><!-- wiki_info -->
		<iframe name="iframe" id="iframe" width="1" height="1" onload="articleForm.onIframeLoaded();"></iframe>
		<div id="demo" class="yui-navset"> 
	    <ul class="yui-nav"> 
	    	<li class="selected"><a href="#tab1"><em>Manual submit</em></a></li>
	    	<xsl:if test="article/mediaType='article'">
		    	<xsl:if test="article/@id &lt; 0">
		        	<li><a href="#tab2"><em>PubMed submit</em></a></li>
		        </xsl:if>
		    	 <xsl:if test="article/@id &lt; 0">
		        	<li><a href="#tab3"><em>External format upload</em></a></li>
		        </xsl:if> 
	        </xsl:if>
	        <xsl:if test="article/mediaType='book'">
	        	 <xsl:if test="article/@id &lt; 0">
		        	<li><a href="#tab4"><em>ISBN-DB submit</em></a></li>
		        </xsl:if>
	        </xsl:if>
	    </ul> 
		<div class="yui-content">
			<div>
				<form enctype="multipart/form-data" onSubmit="return TestDataCheck()" name="articleEditor" action="article/action.do?out=json&amp;action=edit" method="post" target="iframe" accept-charset="utf-8">
				<table class="outer">
					<tr>
						<td class="inner formscope">
							<input type="hidden" name="id" value="{article/@id}" send="1"/>
							<input type="hidden" name="media-type" value="{article/mediaType}" send="1"/>
							<table>
								<tr><td><b>Title</b></td><td><input type="text" name="n-title" send="1" value="{article/@title}" class="w600" restricted="1"/></td></tr>
								<xsl:if test="article/mediaType='article'">
									<tr><td><b>Abstract</b></td><td><textarea name="n-abstract" send="1" class="w600" restricted="1"><xsl:value-of select="article/articleAbstract"/></textarea></td></tr>
								</xsl:if>
								<xsl:if test="article/mediaType='book'">
									<tr><td></td><td>
										<input type="hidden" name="parent-book-id" send="1" value="{article/parent/@id}"/>
										<input type="checkbox" name="is-chapter" send="1">
											<xsl:if test="article/isChapter='true'">
												<xsl:attribute name="checked">checked</xsl:attribute>
											</xsl:if>
										</input>This is a chapter of a book 
										<a action="selectparentbook" id="selectparentbook" bindto="parent-book-id">
											<xsl:choose>
												<xsl:when test="article/parent">[<xsl:value-of select="article/parent/@title"/>]</xsl:when>
												<xsl:otherwise>[...]</xsl:otherwise>
											</xsl:choose>
										</a>
									</td></tr>
								</xsl:if>
							</table>			
						</td>
					</tr>
					<tr>
						<td class="inner authorsscope">
						<input type="hidden" name="id" value="{article/@id}" filter="1"/>
						<table width="100%">
							<tr>
								<td width="80">Authors:</td><td><a action="new">[add]</a>
								<xsl:if test="article/@id &gt; 0">
								<a action="swap">Swap</a>
								<a help="swap-help" class="infolink"></a>
								<div id="swap-help" class="invisible">
								    swap last name and first name
								</div>
								</xsl:if>
								<small>eg. <i>LastName, initial.</i></small></td>
							</tr>
							<tr>
							<td></td>
							<td>
								<div id="ArticleAuthorsBrowser"></div>
							</td>	
							</tr>				
						</table>
						</td>
					</tr>
					
					<tr>				
						<td class="inner formscope">
						<table>
						<xsl:if test="article/mediaType='article'">
							<tr>
								<td>
									Journal/Book: 
								</td>
								<td>
								
									<div id="ac-source" style="height: 25px; width: 320px;"></div>
									<input type="hidden" name="journal-id" value="{article/journal/@id}" />
									<input type="hidden" name="journal-title" value="{article/journal/title}" />
									<a action="editjournal" name="journal-link">
									    [select]
									</a>
									<div id="journal-info"></div>
								</td>
							</tr>
						</xsl:if>
						<tr class="book-only">
							<td>Publication date</td>
							<td><input type="text" name="n-date" class="w100" send="1" value="{article/publication-date/@printed-name}"/></td>
						</tr>
						<tr class="book-only">					
							<td>Volume</td>
							<td><input type="text" name="n-volume" class="w100" send="1" value="{article/volume}"/></td>	
						</tr>
						<xsl:if test="article/mediaType='article'">
						<tr>					
							<td>Issue</td>
							<td><input type="text" name="n-issue" class="w100" send="1" value="{article/issue}"/></td>
						</tr>
						</xsl:if>
						<tr>					
							<td>Page numbers</td>
							<td><input type="text" name="n-pages" class="w100" send="1" value="{article/pageNumbers}"/></td>													
						</tr>
						</table>											
						</td>
					</tr>
					<tr><td class="inner formscope">
					<table width="100%">
						<xsl:if test="article/mediaType='article'">
						<tr>					
							<xsl:choose>
							<xsl:when test="article/pmid">
							
								<input type="hidden" send="1" name="n-pubmed" value="{article/pmid}"/>
								<td>PubMed ID:</td>
								<td><b><a target="_blank" href="http://www.ncbi.nlm.nih.gov/pubmed/{article/pmid}"><xsl:value-of select="article/pmid"/></a></b><input type="button" action="reloadpm" value="reload"/> <a help="Article-reload-hint" class="infolink"> </a> </td>
							</xsl:when>
							<xsl:otherwise>
							<td colspan="2">
								No PubMed ID for this article
								<a id="inlink" action="introducepmid">[Introduce]</a>
								<div id="inblock" class="invisible">
								<input type="hidden" send="1" name="art-id" value="{article/@id}"/>
								<input type="text" send="1" name="n-pubmed" value="{article/pmid}"/>
								<input type="button" action="reloadpm" value="reload"/> <a help="Article-reload-hint" class="infolink"> </a>
								</div>
							</td>
							</xsl:otherwise>
							</xsl:choose>
							<span class="invisible" id="Article-reload-hint">
								With the reload button information about the article is fetched from PubMed.<br/> 
								Existing information will be overwritten!<br/>
							</span>
						</tr>
						</xsl:if>
						<tr>
							<td>OCHEM ID:</td>		
							<td>A<xsl:value-of select="article/@id"/></td>
						</tr>
					</table>
					</td></tr>
					<tr><td class="inner formscope">
					<table width="100%">
						<xsl:if test="article/mediaType='article'">
						<tr>					
							<td>Affiliation</td>
							<td><input type="text" name="n-affiliation" value="{article/affiliation}" class="w560" send="1"/></td>													
						</tr>
						</xsl:if>
						<xsl:if test="article/mediaType='book'">
							<tr class="book-only">
								<xsl:choose>
							<xsl:when test="article/isbn13">
								<input type="hidden" send="1" name="n-isbn13" value="{article/isbn13}"/>
								<td>ISBN-13:</td>
								<td><b><a target="_blank" href="http://isbndb.com/search-all.html?kw={article/isbn13}"><xsl:value-of select="article/isbn13"/></a></b></td>
							</xsl:when>
							<xsl:otherwise>
							<td>ISBN-13</td>
							<td><input type="text" name="n-isbn13" value="{article/isbn13}" class="w100" send="1"/></td>
							</xsl:otherwise>
							</xsl:choose>
							</tr>
							<tr class="book-only">					
								<td>ISBN-10:</td>
								<td><input type="text" name="n-isbn" value="{article/isbn}" class="w100" send="1"/></td>													
							</tr>
							<tr class="book-only">					
							<td>Publisher</td>
							<td><input type="text" name="n-publisher" value="{article/publisher}" class="w560" send="1"/></td>													
						</tr>
						</xsl:if>
						<tr>					
							<td>
								<xsl:choose>
									<xsl:when test="string-length(article/url) > 0">
										<a target="_blank" href="{article/url}">Preprint/Public URL</a>
									</xsl:when>
									<xsl:otherwise>
										Preprint/Public URL
									</xsl:otherwise>
								</xsl:choose>
							</td>
							<td><input type="text" name="n-url" value="{article/url}"  class="w560" send="1" n-disable="1"/></td>													
						</tr>
						<tr class="book-only">					
							<td>
							<xsl:choose>
							<xsl:when test="article/doi">
							<a target="_blank" href="http://dx.doi.org/{article/doi}">DOI</a>
							</xsl:when>
							<xsl:otherwise>
							DOI
							</xsl:otherwise>
							</xsl:choose>
							</td>
							<td>http://dx.doi.org/<input type="text" name="n-doi" value="{article/doi}" class="w470" send="1" n-disable="1"/></td>													
						</tr>
						<xsl:if test="//model/session/user/@id">
							<xsl:choose>
							<xsl:when test="article/@pdf-available ='true'">
								<input type="hidden" name="n-delpdf" value=""/>
								<tr id="pdfavail">
									<td colspan="2"><b>PDF Available:</b><a href="pdf/show.do?id={article/@id}&amp;type=pdf">[show PDF]</a><a action="deletepdf" class="delete-link">[x]</a></td>						
								</tr>					
								<tr id="nopdfavail" class="hidden">
									<td><b>Upload PDF</b></td><td><input type="file" name="file" width="100%"/></td> 
								</tr>					
							</xsl:when>
							<xsl:otherwise>
								<tr>
									<td><b>Upload PDF</b></td><td><input type="file" name="file" width="100%"/></td> 
								</tr>					
							</xsl:otherwise>
							</xsl:choose>	
							<xsl:if test="article/@batch-file ='true'">
								<input type="hidden" name="n-batch" value=""/>
								<tr id="batchavail">
									<td colspan="2"><b>Batch File Available:</b><a href="pdf/show.do?id={article/@id}&amp;type=excel">[show Batch File]</a><a action="deleteexcel" class="delete-link">[x]</a></td>						
								</tr>					
							</xsl:if>
						</xsl:if>
					</table>
					</td></tr>
					<tr>					
						<td class="inner formscope" colspan="2">
							Comment:<br/>
							<textarea maxlength="255" onkeyup="return ismaxlength(this)" name="n-comment" send="1"
								style="width: 100%; height: 80px;" n-disable="1">
								<xsl:value-of select="article/comment" />
							</textarea>
							<div class="message">(max 255 characters)</div>
						</td>
					</tr>
					
					<tr>
						<td class="inner">
							<xsl:if test="article/@property-record &gt; 0">
							This article is referenced from <a href="epbrowser/show.do?article={article/@id}&amp;approval-status=all" tab="Article compounds"><xsl:value-of select="article/@property-record"/> experimental records </a>
							</xsl:if> 
							
							<xsl:if test="article/@tasks-count &gt; 0">
							<br/>
							This article is related to <a href="pendingtasks/published.do?article-id={article/@id}" tab="Article-related tasks"><xsl:value-of select="article/@tasks-count"/> published calculation tasks </a>
							</xsl:if>						
						</td>					
					</tr>
					<xsl:if test="article/@model-list &gt; 0">
						<tr>
							<td class="inner">
								Browse models: <a href="model/select.do?article={article/@id}" tab="Article models"><xsl:value-of select="article/@model-list"/> models</a> 						
							</td>					
						</tr>
					</xsl:if>	
				</table>
				<div class="formscope right"><input type="button" action="cancel" value="cancel"/><input type="submit" value="save"/></div>
				</form>
			</div>
			<xsl:if test="article/mediaType='article'">
				<xsl:if test="article/@id &lt; 0">
					<div>
						<form enctype="multipart/form-data" name="articleEditor" action="article/action.do?out=json&amp;action=edit&amp;pubmed=true" method="post" target="iframe">
						<input type="hidden" name="id" value="-1"/>
						<input type="hidden" name="media-type" value="{article/mediaType}" send="1"/>
						<table class="outer">
							<tr><td class="inner formscope">
								<table width="100%">
									<tr>
										<td><b>Preload from PubMed</b></td><td><input type="text" send="1" name="n-pubmed" width="100%"/><input type="submit" value="load"/></td> 
									</tr>					
								</table>						    
							</td></tr>			
						</table>
						</form>
					</div>
				</xsl:if>
				<xsl:if test="article/@id &lt; 0">
					<div>
						<form enctype="multipart/form-data" name="articleEditor" action="article/upload.do?out=json" method="post" target="iframe">
						<input type="hidden" name="id" value="-1"/>
						<input type="hidden" name="media-type" value="{article/mediaType}" send="1"/>
						<select name="format">
							<option value="ris">RIS</option>
							<option value="en">EndNote - plain text</option>
							<option value="ex">EndNote - xml</option>
<!-- 							<option value="bt">BibTex</option> -->
							<option value="isi">ISI</option>
						</select>
						<table class="outer">
							<tr><td class="inner formscope">
								<table width="100%">
									<tr>
										<td><b>Load from file:</b></td><td><input type="file" name="file" width="100%"/><input type="submit" value="load"/></td> 
									</tr>	
									<tr>
										<td colspan="2" style="font-size: 12px; color: red"><input type="checkbox" name="overwrite">check will allow you to overwrite existing article. </input></td>
									</tr>				
								</table>						    
							</td></tr>
							<tr><td style="font-size: 12px; color: red">
								*Its highly recommended to verify data manually after uploading article from external source.
							</td></tr>		
						</table>
						</form>
					</div>
				</xsl:if>	
			</xsl:if>	
			<xsl:if test="article/mediaType='book'">
				<xsl:if test="article/@id &lt; 0">
					<div>
						<input type="hidden" name="id" value="-1"/>
						<input type="hidden" name="media-type" value="{article/mediaType}" send="1"/>
						<table class="outer">
							<tr><td class="inner formscope">
								<table width="100%">
									<tr>
										<td><b>Preload from ISBNDB</b></td><td><input id="load-isbn" type="text" send="1" name="load-isbn" width="100%"/><input type="button" action="load_isbn" value="load"/></td> 
									</tr>					
								</table>						    
							</td></tr>
							<tr><td style="font-size: 12px; color: red">
								*Its highly recommended to verify data manually after uploading article from external source.
							</td></tr>		
						</table>
					</div>
				</xsl:if>
				</xsl:if>	
		</div>
		</div>
		
		<div id="articleDialog"> 
		    <div class="hd">Message</div> 
	   		 <div class="bd"> 
		        <form> 
		        	<input type="hidden" name="art-id" value=""></input> 
		            <div id="dupArticle">duplicate article</div>
		        </form> 
		    </div> 
		</div>
	</xsl:template>
</xsl:stylesheet>