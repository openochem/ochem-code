function DotprofileBrowser()
{
	var self = this;
	
	this.controller = "modeldot";
	this.editController = "eprecord";
	this.scope  = ".dotbrowser";
	
	UnifiedBrowser.call(this);
	
	this.pager.useHistory = false;
	
	this.filters.scope = this.scope;
	this.container = "DotBrowser";
	this.itemTemplate = "js/templates/dotprofile.ejs";
	this.itemElement = "exp-property";
	this.filters.firstTime = true;
	this.minusCounter = -1;
	
	this.onItemDrawn = function()
	{
		this.currentBlock.find('a[title]').tooltip({showURL: false});
	}
	
	this.doZoom = function(link)
	{
		var img = $(link).find('IMG');
		var id = img.attr("id");
		if (img.hasClass('big'))
		{
			img.removeClass('big');
			img.attr("src","depiction.jsp?id="+id);
		} else
		{
			img.attr("src","depiction.jsp?w=300&h=300&id="+id);
			img.addClass('big');
		}
	}
	
	this.onExcludeSuccess = function()
	{
		window.alert("Record has been excluded from training set successfully. Recalculate the model to see updated statistics");
	}
	
	this.onIncludeSuccess = function()
	{
		window.alert("Record has been included in training set successfully. Recalculate the model to see updated statistics");
	}
	
	this.doEdit = function()
	{
		var id = (this.currentRecordId) ? this.currentRecordId : "-1";
		var editQueryString = this.editQueryString || "";
		var editController = this.editController ? this.editController : this.controller;
		var printedItemName = this.printedItemName || this.itemElement;
		
		var win = openTab("Edit "+printedItemName, webRoot + editController + "/edit.do?render-mode=popup&id=" + id + "&" + editQueryString);
		win.callback = function(entity) 
		{
			win.closeTab();
		};
	}
	
	this.listenEvent('items_loaded', function(){
		self.mainDiv.find('div[rec-id]').eq(0).addClass("selected-point");
	});
}

function DescriptorsBrowser()
{
	var self = this;
	
	this.controller = "modeldot";
	this.scope  = ".descriptorsbrowser";
	
	
	UnifiedBrowser.call(this);
	this.pager.useHistory = false;

	this.filters.scope = this.scope;

	this.container = "DescriptorsBrowser";
	this.itemTemplate = "js/templates/dotdescriptors.ejs";
	this.itemElement = "string";
	this.url = "modeldot/descriptors.do";
	this.filters.firstTime = true;
	this.minusCounter = -1;
	
	this.onItemDrawn = function()
	{
		this.currentBlock.find('a[title]').tooltip({showURL: false});
		//this.currentBlock.css({float: "left", width: "100px"});
	}
	
	this.doZoom = function()
	{
	}
}