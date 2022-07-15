function ModelRecordBrowser()
{
	var self = this;
	
	this.controller = "modeldot";
	this.scope  = ".errorbrowser";
	this.type = "";
	Browser.call(this);
	
	this.pager.useHistory = false;
	
	this.filters.scope = this.scope;
	
	this.container = "ErrorBrowser";
	this.itemTemplate = "js/templates/modelrecord.ejs";
	this.itemElement = "exp-property";
	
	this.url = "modeldot/errorlist.do";
	this.filters.firstTime = true;
	this.minusCounter = -1;
	
	this.getAdditionalData = function() //A stub, if a browser needs to add some data
	{
		return "&type="+self.type;
	}
	
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
}