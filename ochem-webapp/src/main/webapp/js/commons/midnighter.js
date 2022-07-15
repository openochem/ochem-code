	// Support for HTML style content in new jQuery UI widget
	$.widget("ui.tooltip", $.ui.tooltip, {
	    options: {
	        content: function () {
	            return $(this).prop('title');
	        }
	    }
	});

// List of supported browsers
// Order matters here!
include.opera();
var BSupport = new Object();
BSupport.yes = "<font color='green'>Supported</font>";
BSupport.part = "<font color='#CCCC00'>Partially supported</font>";
BSupport.no = "<font color='red'>Not supported!</font>";
BSupport.NA = "<font color='red'>Not tested</font>";
var supportedBrowsers = [
	["Explorer 8", "", BSupport.part], 
	["Explorer 7", "", BSupport.part], 
	["Explorer 6", "", BSupport.part],
	["Chrome 1", "", BSupport.part],
	["Explorer", "", BSupport.no],
	["Opera", "", BSupport.no], 
	["", "", BSupport.NA]
];

var isBrowserSupported = function()
{
	}

var QSPR = function() {};

// A stub for console in IE
var console = console || {
	log: function(){
		
	}
};

QSPR.Ajax = function(url)
{
	this.url = url;
	this.dataType = "json";
	var self = this;
	
	this.send = function(params)
	{
		if (!params.error)
			params.error = this.onError;
		if (!params.success)
			params.success = this.onSuccess;
		if (!this.owner)
			this.owner = this;
		var url = (params.url != undefined) ? params.url : this.url;
		var self = this;
		
		if (params.tracker)
		{
			params.tracker.setStatus("Querying results..");
			params.tracker.start();
			params.tracker.show();
		}
		
		var wrappedSuccess = function(data)
		{
			if (self.waitingDialog)
				self.waitingDialog.hide();
			
			if (params.tracker)
				params.tracker.stop();
			
			var json = data;
			if (self.dataType == "xml")
			{
				json = xml2json(data, "\t");
				json = eval('('+json+')');
				json = json.model;
			}
			
			if (typeof(json.message) != "undefined" && (json.message.type == "exception" || json.message.type == "error"))
				params.error.call(self.owner, json.message.title, json.message);
			else
			{
				params.success.call(self.owner, json);
				if (params.tracker)
					params.tracker.hide();
			}
			
			if (params.after)
				params.after.call(self.owner, data);
			
			if (self.afterRequest)
				self.afterRequest.call(self.owner);
		}
		
		var wrappedError = function(data)
		{
			if (params.tracker)
			{
				params.tracker.stop();
				params.tracker.error();
			}
			params.error.call(self.owner, data);
			if (params.after)
				params.after.call(self.owner, data);
			if (self.afterRequest)
				self.afterRequest.call(self.owner);
		}
		
		
		if (self.beforeRequest)
			self.beforeRequest.call(self.owner);
		
		this.waitingDialog = this.waitingDialog || false;
		if (self.waitingDialog)
			this.waitingDialog.show();
		
		return $.ajax
		(
			{ 
				type: 'POST', 
				url: url, 
				data: "out="+self.dataType+"&"+params.data,
				dataType: self.dataType, 
				success: wrappedSuccess,
				error: wrappedError
			}
		);
	}
	
	this.onError = function(message, messageEntity)
	{
		if (!message)
			message = "";
		if (message instanceof Object)
			window.alert("Request to server failed!");
		else if (message == "")
			window.alert('Error occured while executing this action');
		else
			self.showError(messageEntity);
	}
	
	//this.showError = function(messageEntity)
	//{
	//	// Default error display behavior / Midnighter
	//	window.alert(messageEntity.title);
	//}
	
	this.showError = function(messageEntity)
	{
		// Extended error display. Makes difference between UserFriendlyExceptions and all others, treated as unexpected / Midnighter
		if (messageEntity.context == "friendly")
			window.alert(messageEntity.title);
		else
			window.alert("Dear user,\n\nunexpected error occured during processing your request. "
				+"Please make a screenshot report a bug on email info@ochem.eu.\n\nAdditional information: \n\n"+messageEntity.message.substring(0, 500)+"\n(...)");
	}
}


$.fn.setClass = function(className, condition)
{
	if (condition)
		this.addClass(className);
	else
		this.removeClass(className);
}

$.fn.setChecked = function(condition)
{
	if (condition)
		this.attr("checked", true);
	else
		this.removeAttr("checked");
	
	return this;
}

$.fn.setDisabled = function(condition)
{
	if (condition)
		this.attr("disabled", true);
	else
		this.removeAttr("disabled");
	
	return this;
}

setCheckbox = function(name, value)
{
	var cb = $("input[name="+name+"]");
	if (value == 'true' || value == '1')
		cb.attr("checked", "checked");
	else
		cb.removeAttr("checked");
	cb.change();
}

setRadio = function(name, value)
{
	$("input[name=" + name + "][value='" + value + "']").attr("checked", "checked");
}

setValue = function(name, value)
{
	if (value == '')
		return;
	var elem = $("[name="+name+"]").val(value);
	if (elem && elem.change)
		elem.change();
}

function include_js(script_filename) 
{
    var html_doc = document.getElementsByTagName('head').item(0);
    var js = document.createElement('script');
    js.setAttribute('language', 'javascript');
    js.setAttribute('type', 'text/javascript');
    js.setAttribute('src', script_filename);
    html_doc.appendChild(js);
    return false;
}

function resizeOuterTo(w,h) 
{
 if (parseInt(navigator.appVersion)>3) {
   if (navigator.appName=="Netscape") {
    top.outerWidth=w;
    top.outerHeight=h;
   }
   else top.resizeTo(w,h);
 }
}

function getWikiUrl(str)
{
	if (str != undefined)
		return str.replace(/[ \.]/g, "_")
}

function getFilter(name)
{
	return $('[filter][name="'+name+'"]').attr('value');
}

	var getParams = new Array(); 
  
  	// Parse GET params
function retrieveGETqs() 
{ 
    var query = window.location.search.substring(1); 
    var parms = query.split('&'); 
    for (var i=0; i<parms.length; i++) { 
        var pos = parms[i].indexOf('='); 
        if (pos > 0) { 
            var key = parms[i].substring(0,pos); 
            var val = parms[i].substring(pos+1); 
            getParams[key] = val; 
        } 
    } 
}
retrieveGETqs();

function Utf8encode(string) 
{
	string = string+"";
	string = string.replace(/\r\n/g,"\n");
	var utftext = "";
	 
	for (var n = 0; n < string.length; n++) 
	{
 		var c = string.charCodeAt(n);
 		if (c < 128) 
 		{
 			utftext += String.fromCharCode(c);
		}
		else if((c > 127) && (c < 2048)) {
			utftext += String.fromCharCode((c >> 6) | 192);
			utftext += String.fromCharCode((c & 63) | 128);
		}
		else {
			utftext += String.fromCharCode((c >> 12) | 224);
			utftext += String.fromCharCode(((c >> 6) & 63) | 128);
			utftext += String.fromCharCode((c & 63) | 128);
		}
 
	}
	return utftext;
}
	
function Utf8decode(utftext) 
{
	var string = "";
	var i = 0;
	var c = c1 = c2 = 0;
 
	while ( i < utftext.length ) 
	{
 		c = utftext.charCodeAt(i);
 
		if (c < 128) {
			string += String.fromCharCode(c);
			i++;
		}
		else if((c > 191) && (c < 224)) {
			c2 = utftext.charCodeAt(i+1);
			string += String.fromCharCode(((c & 31) << 6) | (c2 & 63));
			i += 2;
		}
		else {
			c2 = utftext.charCodeAt(i+1);
			c3 = utftext.charCodeAt(i+2);
			string += String.fromCharCode(((c & 15) << 12) | ((c2 & 63) << 6) | (c3 & 63));
			i += 3;
		} 
	}
	return string;
}
	
function URLEncode(plaintext)
{
	//return escape(Utf8encode(plaintext));
	return encodeURIComponent(plaintext);
};

function URLDecode(encoded)
{
	//return Utf8decode(unescape(encoded));
	return decodeURIComponent(encoded);
};

function clone(myObj)
{
	if(typeof(myObj) != 'object') return myObj;
	if(myObj == null) return myObj;

	var myNewObj = new Object();

	for(var i in myObj)
		myNewObj[i] = clone(myObj[i]);

	return myNewObj;
}

//NoS 11.05.12 A huge change... watcho out for trouble!
//Rob 17.07.12 I encountered a problem with XMLElement List<Long>. If there was only one element in the list, array(idList) --> [3,2,3,4,6,8,2],
//			   i.e the single figures of the ID as array. With lists of other entities there were no problems 
//NoS&Rob 16.05.13  if (! (val instanceof Array)): array not recognized
//					if (! (val.length)): one element array (long) is recognized as char[] and therefore not wrapped (s. above)
// 					if (! (val.constructor.name == "Array")) might work
function array(val)
{
	if (val == undefined)
		return new Array();
	//if (! (val instanceof Array))
	//if (! (val.length))
    if (!(Object.prototype.toString.call(val) == "[object Array]"))
		return new Array(val);
	return val;
}

if (window.focus) 
	window.focus();
	
String.prototype.trim = function() {
	return this.replace(/^\s+|\s+$/g,"");
}
String.prototype.ltrim = function() {
	return this.replace(/^\s+/,"");
}
String.prototype.rtrim = function() {
	return this.replace(/\s+$/,"");
}

function concatObject(obj) {
	  str='';
	  for(prop in obj)
	  {
	    str+=prop + " value :"+ obj[prop]+"\n";
	  }
	  return(str);
	}

function createCookie(name,value,days) {
	if (days) {
		var date = new Date();
		date.setTime(date.getTime()+(days*24*60*60*1000));
		var expires = "; expires="+date.toGMTString();
	}
	else var expires = "";
	document.cookie = name+"="+value+expires+"; path=/";
}

function readCookie(name) {
	var nameEQ = name + "=";
	var ca = document.cookie.split(';');
	for(var i=0;i < ca.length;i++) {
		var c = ca[i];
		while (c.charAt(0)==' ') c = c.substring(1,c.length);
		if (c.indexOf(nameEQ) == 0) return c.substring(nameEQ.length,c.length);
	}
	return null;
}



var fontSize = readCookie("fontsize");
if (!fontSize)
	fontSize = 100;
else
	fontSize = 1*fontSize;

function changeFontSize(delta)
{
	fontSize += delta;
	$('body').css({'font-size':fontSize + '%'});
	createCookie('fontsize', fontSize, 60);
	$(".main-tabs iframe").each(function(){
		if (this.contentWindow)
			if (this.contentWindow.changeFontSize)
				this.contentWindow.changeFontSize(delta);
	});
}

function getBrowserInfo()
{
	var userAgent = navigator.userAgent.toLowerCase();
}

$(document).ready(function()
{
	// After tabs it stopped working correctly
	// Requires refactoring - history of tabs etc.
	//YAHOO.util.History.initialize("yui-history-field", "yui-history-iframe");
	changeFontSize(0);
	
	if ($("#main-tab-container").length > 0)
	{
		addTabs();
		
		var compatibility = getCompatibilityStatus();
		if (compatibility.indexOf("Supported") < 0)
			$("#browser-version").html(BrowserDetect.browser + " " + BrowserDetect.version + " on " + BrowserDetect.OS + " - <a href='https://docs.ochem.eu/display/MAN/Supported+Browsers' target='_blank'>" + compatibility +"</a>");
		else
			$("#browser-version").remove();
	}
	
	$(document).trigger("DOM_updated", $(document));
	$('a[help]').each(function(){
		$(this).attr('title', $("#" + $(this).attr('help')).html());
		$(this).tooltip({
			showURL: false, 
			delay: 0, 
			content: $("#" + $(this).attr('help')).html()
		});
	});
	$(":not(div)[title]").filter(":not([help])").tooltip({showURL: false});
	$(".infolink[href]").filter(":not([title])").attr("title", "Click to read");
});

$(document).bind('DOM_updated', null, function()
{
	$("a[tab]").each(function()
	{
		if (!this.tabHandler)
		{
			this.tabHandler = true;
			$(this).click(function(){
				
				openTab($(this).attr('tab'), $(this).attr('href'), $(this).attr("waitMsg"));
				return false;
			});
		}
	});
	
	$(":not(div)[title]").filter(":not([help])").tooltip({showURL: false});
});


// Here goes experimental part
// "Tabs" feature, suggested by Wladimir
// should it be a success, we will move it to appropriate place
// Midnighter

function addTabs()
{
	QSPR.tvMain = new YAHOO.widget.TabView('main-tab-container');
	var iframe = $("#main-tab-container > div > div:last").find("iframe").get(0);
	var firstTab = QSPR.tvMain.get('tabs')[0];
	$(iframe).load(function(){
		this.firstTab = firstTab;
		tabLoaded(firstTab, this);
	});
	
	resizeTabs();
	
	$(window).bind('resize', function() 
	{
		resizeTabs();
	});
}

function resizeTabs()
{
	// Adjust size of main frame
	var size = Math.round($("#bottom-aligner").position().top - $("TD.content").position().top - $("DIV.main-tabs .yui-content").position().top);
	//window.alert("" + Math.round(""+ $("#bottom-aligner").position().top) + " - " + $("TD.content").position().top + " - " + $("DIV.main-tabs .yui-content").position().top);
	
	if (size > 0)
		$(".main-tabs iframe").css({height: ""+size+"px"});
	
	//console.log($(document).hasScrollBar());
	
	//if ($(document).hasScrollBar())
	//{
	//	size -= 100;
	//	$(".main-tabs iframe").css({height: ""+size+"px"});
	//}
}

(function($) {
    $.fn.hasScrollBar = function() {
        return this.get(0).scrollHeight > this.height();
    }
})(jQuery);

var handleClose = function(e, tab) 
{ 
	YAHOO.util.Event.preventDefault(e); 
	
	// Call the handler of this event, if it is defined in a window
	try
	{
		if (tab.contentWindow)
			if (tab.contentWindow.onTabClose)
				if (!tab.contentWindow.onTabClose())
					return;
	}
	catch (err)
	{
		// 
	}
	
	if (QSPR.tvMain.get('tabs').length <= 2)
	{
		$("div.main-tabs > ul").addClass("invisible");
		resizeTabs();
	}
	
	QSPR.tvMain.removeTab(tab); 
};

function tabLoaded(tab, iframe)
{
	try {
		iframe.contentWindow.iFrame = iframe;
	}
	catch(e) { //ignore permission denied error for wiki pages.
		//IE displays permission denied error for loading cross domain pages which is an expected security behaviour of IE.
		//since the wiki domain is trusted. we can ignore this error.
	}
	
	if (tab.originalLabel == "" || tab.originalLabel == undefined || tab.originalLabel == "undefined" || !tab.originalLabel)
	{
		var title = "OCHEM";
		try {
			title = iframe.contentWindow.$("title").html();
		}
		catch (e) {
		}
		tab.set('label', title+"<span class='close'>X</span>");
		YAHOO.util.Event.on(tab.getElementsByClassName('close')[0], 'click', handleClose, tab);
	}
}

//$(window).bind('beforeunload', function() {
//    if (QSPR.tvMain)
//    {
//    	if (QSPR.tvMain.get('tabs').length >= 2)
//    		return "You have " + QSPR.tvMain.get('tabs').length +" tabs open.\nPushing the BACK button will close the tabs.\n\nAre you sure you would like to continue?";
//    }
//}); 

function openTab(label, url, waitMsg)
{
	if (!QSPR.tvMain)
	{
		var wnd = window.parent.openTab(label, url, waitMsg);
		wnd.openerTab = window;
		return wnd;
	}
	else
	{
		var tab;
		
		// Show tab bar
		$("div.main-tabs > ul").removeClass("invisible");
		
		if (url.indexOf("?") == -1)
			url += "?render-mode=popup";
		else
			url += "&render-mode=popup";
		QSPR.tvMain.addTab(tab = new YAHOO.widget.Tab(
				{ 
					label: label+"<span class='close'>X</span>", 
					content: "<div class='frame-status invisible'><span></span><br/><img src='img/roller_transparent.gif'/></div><iframe frameborder='0' src='" + url + "' class='iframe'></iframe>"
				}));
		tab.originalLabel = label;
		YAHOO.util.Event.on(tab.getElementsByClassName('close')[0], 'click', handleClose, tab);
		
		QSPR.tvMain.set("activeIndex", QSPR.tvMain.get('tabs').length - 1);
		var iframe = $("#main-tab-container > div > div:last").find("iframe").get(0);
		var curWindow = window;
		
		if (waitMsg != undefined)
		{
			$(iframe).parent().find("div.frame-status").removeClass("invisible").find("span").html(waitMsg);
		}
		
		
		$(iframe).load(function(){
			$(this).parent().find("div.frame-status").addClass("invisible");
			tab.contentWindow = this.contentWindow;
			this.tab = tab;
			this.contentWindow.opener = curWindow;
			tabLoaded(tab, this);
			
			try {
				this.contentWindow.tab = tab;
			}
			catch(e) { //ignore permission denied error for wiki pages.
				//IE displays permission denied error for loading cross domain pages which is an expected security behaviour of IE.
				//since the wiki domain is trusted. we can ignore this error.
			}
			
		});
		
		resizeTabs();
		return iframe.contentWindow;
	}
}

// Setup a tracker for long-loading iframes
function frameLoading(iFrame)
{
	if (!iFrame)
	{
		if (!window.iFrame)
			throw "No parent window defined!";
		return window.parent.frameLoading(window.iFrame);
	}
	
	if (!iFrame.loadHookSet)
	{
		iFrame.loadHookSet = true;
		$(iFrame).load(function(){
			this.operation.stopCheckingStatus(); // Once frame is loaded, stop the tracker
		});
	}
	
	$(iFrame).parent().find("div.frame-status").removeClass("invisible").find("span").html("Loading");
	var longOperation = iFrame.operation = new LongOperation({tracker: $(iFrame).parent().find("div.frame-status span")});
	longOperation.operationId = 666; // Testing stage
	longOperation.startCheckingStatus();
}

function closeTab()
{
	// Hide tab bar
	if (window.parent.QSPR.tvMain.get('tabs').length <= 2)
	{
		window.parent.$("div.main-tabs > ul").addClass("invisible");
		window.parent.resizeTabs();
	}
	
	window.parent.QSPR.tvMain.removeTab(window.tab);
}

function closeActiveTab()
{
	// Hide tab bar
	if (window.parent.QSPR.tvMain.get('tabs').length <= 2)
	{
		window.parent.$("div.main-tabs > ul").addClass("invisible");
		window.parent.resizeTabs();
	}
	window.parent.QSPR.tvMain.removeTab(window.parent.QSPR.tvMain.get("activeTab"));
}

/// Browser detection script
// http://www.quirksmode.org/js/detect.html
var BrowserDetect = 
{
		init: function () {
			this.browser = this.searchString(this.dataBrowser) || "An unknown browser";
			this.version = this.searchVersion(navigator.userAgent)
				|| this.searchVersion(navigator.appVersion)
				|| "an unknown version";
			this.OS = this.searchString(this.dataOS) || "an unknown OS";
		},
		searchString: function (data) {
			for (var i=0;i<data.length;i++)	{
				var dataString = data[i].string;
				var dataProp = data[i].prop;
				this.versionSearchString = data[i].versionSearch || data[i].identity;
				if (dataString) {
					if (dataString.indexOf(data[i].subString) != -1)
						return data[i].identity;
				}
				else if (dataProp)
					return data[i].identity;
			}
		},
		searchVersion: function (dataString) {
			var index = dataString.indexOf(this.versionSearchString);
			if (index == -1) return;
			return parseFloat(dataString.substring(index+this.versionSearchString.length+1));
		},
		dataBrowser: [
			{
				string: navigator.userAgent,
				subString: "Chrome",
				identity: "Chrome"
			},
			{ 	string: navigator.userAgent,
				subString: "OmniWeb",
				versionSearch: "OmniWeb/",
				identity: "OmniWeb"
			},
			{
				string: navigator.vendor,
				subString: "Apple",
				identity: "Safari",
				versionSearch: "Version"
			},
			{
				prop: window.opera,
				identity: "Opera"
			},
			{
				string: navigator.vendor,
				subString: "iCab",
				identity: "iCab"
			},
			{
				string: navigator.vendor,
				subString: "KDE",
				identity: "Konqueror"
			},
			{
				string: navigator.userAgent,
				subString: "Firefox",
				identity: "Firefox"
			},
			{
				string: navigator.vendor,
				subString: "Camino",
				identity: "Camino"
			},
			{		// for newer Netscapes (6+)
				string: navigator.userAgent,
				subString: "Netscape",
				identity: "Netscape"
			},
			{
				string: navigator.userAgent,
				subString: "MSIE",
				identity: "Explorer",
				versionSearch: "MSIE"
			},
			{
				string: navigator.userAgent,
				subString: "Gecko",
				identity: "Mozilla",
				versionSearch: "rv"
			},
			{ 		// for older Netscapes (4-)
				string: navigator.userAgent,
				subString: "Mozilla",
				identity: "Netscape",
				versionSearch: "Mozilla"
			}
		],
		dataOS : [
			{
				string: navigator.platform,
				subString: "Win",
				identity: "Windows"
			},
			{
				string: navigator.platform,
				subString: "Mac",
				identity: "Mac"
			},
			{
				   string: navigator.userAgent,
				   subString: "iPhone",
				   identity: "iPhone/iPod"
		    },
			{
				string: navigator.platform,
				subString: "Linux",
				identity: "Linux"
			}
		]
};

BrowserDetect.init();

function checkMessages() {
	if ($("#message-count").length == 0)
		return;
	
	$.ajax({
		url: "user/getSessionData.do?out=json",
		dataType: "json",
		success: function(response) {
			var count = response.session.user.message || 0;
			unreadMessagesCount = count;
			if (count > 0)
				$("#message-count").html("(" + count + " new mails!)");
			else
				$("#message-count").html("");
			setTimeout(checkMessages, 5000);
		},
		error: function() {
			console.log("Error checking messages");
		}
	});
}

$(function(){
	setTimeout(checkMessages, 1000);
});

function getCompatibilityStatus()
{
	var curBrowser = "" + BrowserDetect.browser + " " + BrowserDetect.version;
	
	if (BrowserDetect.browser == "Explorer")
	{
		$(function(){
			$(".announcement").html("You are using Internet Explorer which is not compatible with OCHEM. We recommend Firefox, Chrome or Safari for the best experience.");
			$(".announcement").parent().removeClass("invisible");
		});
		
	}
	
	if (BrowserDetect.browser == "Chrome")
		return BSupport.yes;
	
	if (BrowserDetect.browser == "Firefox")
		if (BrowserDetect.version >= 6)
			return BSupport.yes;
		else if (BrowserDetect.version >= 4)
			return BSupport.part;
		else
			return BSupport.no;
	
	if (BrowserDetect.browser == "Safari")
		if (BrowserDetect.version >= 4)
			return BSupport.yes;
		else
			return BSupport.no;
	
	for (var i = 0; i < supportedBrowsers.length; i++)
	{
		if (curBrowser.indexOf(supportedBrowsers[i][0]) != -1)
			return supportedBrowsers[i][2];
	}
	
	return "This browser has not been fully tested!";
}

function isInteger(value)
{
	var intRegex = /^\d+$/;
	return intRegex.test(value);
}

var tagsToReplace = {
	    '&': '&amp;',
	    '<': '&lt;',
	    '>': '&gt;',
	    '\n': '<br/>'
	};
	
function replaceTag(tag) {
    return tagsToReplace[tag] || tag;
}

function replaceHTMLTags(str) {
	if (!str)
		return "";
	str = str.trim().replace("/\&lt\;br\&gt\;/g", "");
    return str.replace(/[&<>\n]/g, replaceTag);
}

String.prototype.startsWith = function(str)
{return (this.match("^"+str)==str)}

if (!String.prototype.format) {
	  String.prototype.format = function() {
	    var args = arguments;
	    return this.replace(/{(\d+)}/g, function(match, number) { 
	      return typeof args[number] != 'undefined'
	        ? args[number]
	        : match
	      ;
	    });
	  };
	}

var templatesMap = new Object();
var getTemplate = function(id) {
	if (!templatesMap[id])
		return templatesMap[id] = new View({element: id});
	else
		return templatesMap[id];
}

jQuery.fn.hint = function (blurClass) {
	  if (!blurClass) { 
	    blurClass = 'hint-blur';
	  }

	  return this.each(function () {
	    // get jQuery version of 'this'
	    var $input = jQuery(this),

	    // capture the rest of the variable to allow for reuse
	      title = $input.attr('hint'),
	      $form = jQuery(this.form),
	      $win = jQuery(window);

	    function remove() {
	      if ($input.val() === title && $input.hasClass(blurClass)) {
	        $input.val('').removeClass(blurClass);
	      }
	    }

	    // only apply logic if the element has the attribute
	    if (title) { 
	      // on blur, set value to title attr if text is blank
	      $input.blur(function () {
	        if (this.value === '') {
	          $input.val(title).addClass(blurClass);
	        }
	      }).focus(remove).blur(); // now change all inputs to title

	      // clear the pre-defined text when form is submitted
	      $form.submit(remove);
	      $win.unload(remove); // handles Firefox's autocomplete
	    }
	  });
	};
