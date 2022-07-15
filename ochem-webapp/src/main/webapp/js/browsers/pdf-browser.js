include.plugins('view');
var sampleBrowser = new Browser();			
sampleBrowser.itemElement = "row";
sampleBrowser.itemTemplate = "js/templates/pdf.ejs";
sampleBrowser.url = webRoot + "pdf/list.do";
sampleBrowser.actionURL = webRoot + "pdf/action.do";

sampleBrowser.onDeleteSuccess = function()
{
	sampleBrowser.request();
};


$(document).ready(function() {sampleBrowser.initialize();});
