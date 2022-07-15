<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	
	<xsl:template name="content">
	<style type="text/css">
	HTML, BODY, #home {height: 100%;}
	#home TD
	{
		padding: 15px; vertical-align: top;
	}
	
	.bottom {background-color: #EEF; vertical-align: bottom; text-align: center; height: 10px;}
	.bottom A {padding-left: 15px; padding-right: 15px; font-family: Arial;} 
	
	
	#ochem-tweet
	{
	}
	
	.tweet_header
	{
		 text-align: center;
		 font-size: 100%;
		 min-height: 20px;
		 font-family: Georgia;
	}
	
	.tweet,
	.query {
		font: 100% Georgia, serif;
		color: #085258;
	}
	
		.tweet_list, .active-users, .yellow-box, .blue-box {
			-webkit-border-radius: 0.5em;
			-moz-border-radius: 0.5em;
			border-radius: 0.5em;
			list-style: none;
			margin: 0;
			padding: 0;
			overflow-y: hidden;
			background-color: #EEF;
			
		}
		
		.blue-box {padding: 10px;}
		.yellow-box {padding: 10px; text-align: center; background-color: #faf6d7 !important;}
		
			.tweet_list .awesome,
			.tweet_list .epic {
				text-transform: uppercase;
			}
			
			.tweet_list li {
				overflow-y: auto;
				overflow-x: hidden;
				padding: 0.5em;
			}
			
				.tweet_list li a {
					color: #0C717A;
				}
			
			.tweet_list .tweet_even {
				background-color: #91E5E7;
			}
			
			.tweet_list .tweet_avatar {
				padding-right: .5em; float: left;
			}
			
				.tweet_list .tweet_avatar img {
					vertical-align: middle;
				}
				
		#tagCloud { width:290px; background-color:#575454; text-align:center; padding:5px; overflow:auto; font-size:70%; font-family:arial; }
		#tagCloud h2 {font-size:2.5em; margin:0 0 10px 0; background:url(images/cloud.gif) no-repeat 0; padding:15px 0 15px 80px; }
		#tagList { margin:0; padding:0; }
		#tagList li { list-style-type:none; float:left; margin:0 10px; height:35px; }
		#tagList li a { text-decoration:none;}
		#tagList li a:hover { text-decoration:underline; }
		A.reset {font-size: inherit !important; display: inline !important; font-family: inherit !important;}
	</style>
	<script language="javascript" src="js/lib/jquery.tweet.js"></script>
	<script language="javascript">
		jQuery(function($){
	      $("#tweet").tweet({
	         avatar_size: 32,
        	count: 4,
	        query: "#ochem_eu",
	        loading_text: "loading tweets..."
	      });
	    });
	</script>
	<title>OCHEM home page</title>
	<table width="100%" height="100%" id="home">
		<tr>
			<td width="33%">
				<div class="tweet_header">Welcome to OCHEM! Your possible actions</div>
					<ul class="options blue-box">
					<li>
						<a href="epbrowser/show.do">Explore Open OCHEM data</a>
						Search chemical and biological data: 
						experimentally measured, published and exposed to public access by our users. You can also <a href="batchupload30/show.do" class="reset">upload your data.</a>
					</li>
					<li>
						<a href="modelconfigurator/choose.do">Create QSAR models</a>
						Build QSAR models for predictions of chemical properties. The models can be based on the experimental data published in our database.
					</li>
					<li>
						<a href="predictor/show.do">Run predictions</a>
						Apply one of the available models to predict property you are interested in for your set of compounds.
					</li>
					<li>
						<a href="alerts/home.do">Screen compounds with ToxAlerts</a>
						Screen your compound libraries against structural alerts for such endpoints as mutagenicity, skin sensitization, aqueous toxicity, etc.
					</li>
					<li>
						<a href="static/tutorials.do">Tutorials</a>
						Check our video tutorials to know more about the OCHEM features.
					</li>
					<li>
						<a href="static/acknowledgements.do">Our acknowledgements</a>
					</li>		
				</ul>
				<br/>
				<div class="tweet_header">
					Feedback and help
				</div>
				<div class="blue-box">
					<ul class="options">
						<li>
							<a target="_blank" href="https://docs.ochem.eu/display/MAN/OCHEM+Introduction">User's manual</a>
							Check an online user's manual
						</li>
					</ul>
				</div>
			</td>
			<td>
				<div class="tweet_header">Check out the properties available on OCHEM</div>
				<div class="yellow-box">
				OCHEM contains <a href="epbrowser/show.do" tab="Compound properties browser">
					<xsl:value-of select="//param[@key='total-records']"/> records</a> 
					for <a tab="Properties browser" href="properties/show.do"> <xsl:value-of select="//param[@key='total-properties']"/> properties</a> (with at least <xsl:value-of select="//param[@key='minimum-count']"/> records) 
					collected from <a href="article/show.do" tab="Browser of articles">
					<xsl:value-of select="//param[@key='total-articles']"/> sources</a>
				</div>
				<div id="property-cloud" class="yellow-box">
					Loading<br/>
				</div>
			</td>
			<td width="400px">
				<div class="tweet_header">Latest active users</div>
				<div class="active-users">
					<ul class="tweet_list">
						<xsl:for-each select="//others/user">
						<li>
							<a href="user/profile.do?login={@login}" tab="User profile" class="tweet_avatar"><img src="img/user_small.png"/></a>
							<span class="tweet_time"><a href="user/profile.do?login={@login}" tab="User profile"><xsl:value-of select="@login"/>:</a></span>
							<span class="tweet_text"><xsl:value-of select="Title"/>&#160;<xsl:value-of select="FirstName"/>&#160;<xsl:value-of select="LastName"/><br/>
							<xsl:value-of select="latestActivity"/></span>
						</li>
						</xsl:for-each>
					</ul>
				</div>
				<br/>
				<div class="tweet_header">Latest published models</div>
				<div class="models">
					<ul class="tweet_list">
						<xsl:for-each select="//others/model">
						<li>
							<a href="model/profile.do?public_id={publicId}" tab="Model profile" class="tweet_avatar"><img src="img/model.png"/></a>
							<a href="model/profile.do?public_id={publicId}" tab="Model profile"><xsl:value-of select="modelMappings[1]/property/@name"/> model</a> published by 
							<a href="user/profile.do?login={session/user/@login}" tab="User profile"><xsl:value-of select="session/user/@login"/></a><br/>
							<xsl:value-of select="last-modified"/>
						</li>
						</xsl:for-each>
					</ul>
				</div>
			</td>
		</tr>
		
		<tr>
			<td colspan="3" class="bottom">
				<a href="static/acknowledgements.do">Our acknowledgements</a> | <a href="static/tutorials.do">Tutorials</a> | <a target="_blank" href="https://docs.ochem.eu/display/MAN.html">User's manual</a>
			</td>
		</tr>
	</table>
	
	<script language="javascript">
		$.getJSON("home/properties.do?out=json", function(data) {
				
					$("#property-cloud").html("");
					//create list for tag links
					$("<ul/>").attr("id", "tagList").appendTo("#property-cloud");
					
					//create tags
					$.each(data.others.property, function(i, val) 
					{
						val.timesUsed = val.timesUsed;
						var li = $("<li/>");
						
						//create link
						$("<a/>").text(val.name).attr({title: "" + val.timesUsed + " records available", href:"epbrowser/show.do?property=" + val.id}).attr("tab", "Property records").appendTo(li);
						
						//set tag size
						var freq = 1.0 * val.timesUsed;
						freq = Math.log(freq);
						var size = 2*(freq - Math.log(10))/Math.log(100000) + 0.5;
						li.children().css("fontSize", size + "em");
						
						//add to list
						li.appendTo("#tagList");
					});
					
					$(document).trigger("DOM_updated", $(document));
				});
	</script>
	
	</xsl:template>
	
</xsl:stylesheet>
