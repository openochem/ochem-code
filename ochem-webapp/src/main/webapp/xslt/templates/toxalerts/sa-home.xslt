<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />
	
	<xsl:template name="content">
		<style type="text/css">
			DIV.toxalerts {padding: 40px; max-width: 800px;}
			DIV.toxalerts H1 {font-size: 200%; margin-bottom: 13px; border-bottom: 1px solid black;}
			DIV.toxalerts P {margin-bottom: 13px; text-indent:50px;}
			DIV.toxalerts TABLE TD A {font-size: 100%; padding-top: 40px; padding-bottom: 40px; margin-right: 20px;}
			
			.button {
				display: inline-block;
				outline: none;
				cursor: pointer;
				text-align: center;
				text-decoration: none;
				font: 14px/100% Arial, Helvetica, sans-serif;
				padding: .5em 2em .55em;
				text-shadow: 0 1px 1px rgba(0,0,0,.3);
				-webkit-border-radius: .5em; 
				-moz-border-radius: .5em;
				border-radius: .5em;
				-webkit-box-shadow: 0 1px 2px rgba(0,0,0,.2);
				-moz-box-shadow: 0 1px 2px rgba(0,0,0,.2);
				box-shadow: 0 1px 2px rgba(0,0,0,.2);
			}
			.button:hover {
				text-decoration: none;
			}
			.button:active {
				position: relative;
				top: 1px;
			}
			
			.orange {
				color: #fef4e9;
				border: solid 1px #da7c0c;
				background: #f78d1d;
				background: -webkit-gradient(linear, left top, left bottom, from(#faa51a), to(#f47a20));
				background: -moz-linear-gradient(top,  #faa51a,  #f47a20);
				filter:  progid:DXImageTransform.Microsoft.gradient(startColorstr='#faa51a', endColorstr='#f47a20');
			}
			.orange:hover {
				background: #f47c20;
				background: -webkit-gradient(linear, left top, left bottom, from(#f88e11), to(#f06015));
				background: -moz-linear-gradient(top,  #f88e11,  #f06015);
				filter:  progid:DXImageTransform.Microsoft.gradient(startColorstr='#f88e11', endColorstr='#f06015');
			}
			.orange:active {
				color: #fcd3a5;
				background: -webkit-gradient(linear, left top, left bottom, from(#f47a20), to(#faa51a));
				background: -moz-linear-gradient(top,  #f47a20,  #faa51a);
				filter:  progid:DXImageTransform.Microsoft.gradient(startColorstr='#f47a20', endColorstr='#faa51a');
			}
			
			.white {
				color: #606060;
				border: solid 1px #b7b7b7;
				background: #fff;
				background: -webkit-gradient(linear, left top, left bottom, from(#fff), to(#ededed));
				background: -moz-linear-gradient(top,  #fff,  #ededed);
				filter:  progid:DXImageTransform.Microsoft.gradient(startColorstr='#ffffff', endColorstr='#ededed');
			}
			.white:hover {
				background: #ededed;
				background: -webkit-gradient(linear, left top, left bottom, from(#fff), to(#dcdcdc));
				background: -moz-linear-gradient(top,  #fff,  #dcdcdc);
				filter:  progid:DXImageTransform.Microsoft.gradient(startColorstr='#ffffff', endColorstr='#dcdcdc');
			}
			.white:active {
				color: #999;
				background: -webkit-gradient(linear, left top, left bottom, from(#ededed), to(#fff));
				background: -moz-linear-gradient(top,  #ededed,  #fff);
				filter:  progid:DXImageTransform.Microsoft.gradient(startColorstr='#ededed', endColorstr='#ffffff');
			}
		</style>
		<title>Model Templates</title>
		<div class="toxalerts">
			<h1>Welcome to <span class="toxalerts">ToxAlerts</span>!</h1>
			<p>
				<i>Structural alerts</i> (also known as "<i>toxicophores</i>") are molecular patterns known to be associated with particular type of toxicity. 
				The studies performed last decade has shown that structural alerts is an efficient technique to detect potentially toxic chemicals. 
				Screening chemical compounds against known structural alerts can be a good practice to complement the QSAR models and to help interpreting their predictions.
			</p>
			<p>
				<span class="toxalerts">ToxAlerts</span> is a platform for screening chemical compounds against structural alerts. The platform allows to search structural alerts, introduce your own alerts and screen chemical libraries for alert-hitting compounds. 
			</p>
			<br/><br/>
			
			<br/>
			<table align="center">
				<tr>
					<td><a href="alerts/show.do?render-mode=popup" class="white button">View available alerts</a></td>
					<td><a href="alerts/upload.do?render-mode=popup" class="white button">Upload new alerts</a></td>
					<td><a href="alerts/screen.do?render-mode=popup" class="white button">Screen your molecules</a></td>
				</tr>
			</table>
			<br/><br/>
				In case of any questions, ideas, or problems with the software, feel free do <a href="mailto:info@ochem.eu" target="_blank">drop us a message</a>. We highly appreciate any feedback from you!
				
		</div>
		</xsl:template>
</xsl:stylesheet>
