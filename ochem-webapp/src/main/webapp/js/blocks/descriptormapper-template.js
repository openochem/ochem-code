var form = new AjaxForm();
var self = form;
form.doSelectset = function(link)
{
	var targetLink = link;
	var baskWin = openTab("Select compound set", webRoot+"basket/show.do?render-mode=popup");
	baskWin.callback = function(basket)
	{
		var properties;
		self.setValue($(targetLink).attr('bindto'), basket.id, basket.name);
		
		$(document).trigger("DOM_updated", $(document));		
		baskWin.closeTab();
	}
}

form.doSubmit = function()
{
	var ts = self.getValue("trainingsetid");
	var vs = self.getValue("validationsetid");
	if (ts == undefined || ts == "") 
		return window.alert("Please select a valid trainingset for the model");		
	document.modelform.submit();
}

$(document).ready(function()
{
});

