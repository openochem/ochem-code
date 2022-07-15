<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />

	<xsl:template name="content">
		<style type="text/css">
			.actions A
			{
				margin-right: 10px;
			}
		</style>
		<title>Article browser</title>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/browsers/article-browser.js"></script>
		<script language="javascript">
			include.plugins('view');
			var sampleBrowser = new ArticleBrowser();
			sampleBrowser.options.editwindow_scrollbars = 1;
			sampleBrowser.options.editwindow_resizable = 'no';
			$(document).ready(function() {sampleBrowser.initialize();});
		</script>
		<title>Article browser</title>
		<table height="100%">
			<tr>
				<td class="itunes-up" colspan="2">
				<img src="img/icons/article.png"/>
					<xsl:call-template name="area-of-interest"/>
					<h1><doc term="Browser+of+scientific+publications#Browserofscientificpublications-articleBrowserBrowserofarticles">Articles/books browser</doc></h1>
					Search for articles, books and other sources for scientific numerical data.
				</td>
			</tr>
			<tr>
				<td class="itunes-left">
					<div class="openable opened">
					<br/>
						Media type:
						<select name="media-type" filter="1">
							<option value="all">All sources</option>
							<option value="article">Articles</option>
							<option value="book">Books and chapters</option>
							<option value="temporal">Unpublished</option>
						</select>
						<h1>Basic filters</h1>
						<div class="openable-content">
							<div class="book-only">
								ISBN/QID
								<input type="text" name="identifier" prompt="Type to filter by ISBN number" filter="1"/>
							</div>
							Pubmed ID / OCHEM ID
							<input type="text" name="article-identifier" prompt="Type to filter by OCHEM (AXXX) PubMed ID (XXX)" filter="1"/>
							Title
							<input type="text" name="title" prompt="Type any part of title" filter="1"/>
							<br />
							Authors
							<input type="text" name="name" prompt="Type any part of author name" filter="1"/>
							Year
							<input typ="text" name="year" filter="1"/>
						</div>
					</div>
					
					<div class="book-only">
						<input type="checkbox" name="hide-chapters" filter="1"/>Don't show chapters
					</div>
					<div class="openable opened article-only">
						<h1>Journal filters</h1>
						<div class="openable-content">
							Journal
							<input typ="text" name="journal" filter="1"/>
							<br />
							Volume
							<input typ="text" name="volume" filter="1"/>
							<br/>
						</div>
					</div>
					<div class="openable opened">
						<input type="checkbox" name="models" value="true" filter="1"/>With models				
					</div>
					<div class="openable">
						<h1>Additional filters</h1>
						<div class="openable-content">
							Introducer / modifier user name
							<input type="text" name="username" prompt="Type exact user name" filter="1"/>
						</div>
					</div><br/>
					
					<a action="refresh" class="button-link" >refresh</a>
				</td>
				<td class="itunes-right">
					<div class="upper-command-panel">
						<a title="Create new" action="edit">
							[Create new <img src="img/icons/new.gif"/>]
						</a>
					</div><br/>
					
					<div class="pager-strip">
						<span><b class="showed">none</b> of <b class="total">none</b></span>
						<div id="pager" class="pgr">
						</div>
					</div>
					
					<div id="Browser">
					</div>
					
					<div class="pager-strip">
						<span><b class="showed">none</b> of <b class="total">none</b></span>
						<div id="pager" class="pgr">
						</div>
					</div>
				</td>
			</tr>
		</table>
	</xsl:template>
</xsl:stylesheet>