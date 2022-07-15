//ArticleBrowser.Inherits(Browser);
function ArticleBrowser()
{
	this.controller = "article";
	var self = this;
	
	Browser.call(this);
	this.itemElement = "article";	

	this.onDeleteSuccess = function(response)
	{
		this.filters.setValue("deletion-confirmation", "");
		if (response.message)
			// Ask user confirmation bzw. reason for deletion
			this.overviewDeletion(response.message);
		else
			this.request();	
		//this.deleteRecord();
	}	
	
	this.onDeleteError = function(msg)
	{
		this.filters.setValue("deletion-confirmation", "");
		window.alert(msg);
	}
	
	this.overviewDeletion = function(recordsSummary)
	{
		if (window.confirm(recordsSummary.message))
		{
			this.filters.setValue("deletion-confirmation", "OK");
			this.callAction("delete");
		}
	}
	
	this.draw = function(entity)
	{
		return self.views[entity.mediaType].render(entity);
	}
	
	this.onItemSaved = function(entity)
	{	
	    if (entity.oldId < 0)
	    {
	    	self.page = 1;
		    self.filters.setValue("identifier", "Q"+entity.id);
	    }    
	}
	
	this.onMediaTypeChange = function()
	{
		var select = $("[name='media-type']");
		$(".book-only").addClass("invisible");
		$(".article-only").addClass("invisible");
		$(".book-only, .article-only").find("[filter]").attr("filter", "0");
		if (select.val() == "book")
		{
			$(".book-only").removeClass("invisible");
			$(".book-only").find("[filter]").attr('filter', "1")
		}
		if (select.val() == "article")
		{
			$(".article-only").removeClass("invisible");
			$(".article-only").find("[filter]").attr('filter', "1")
		}
		self.editQueryString = "media-type="+select.val();	
	}
	
	this.listenEvent('items_loaded', function(){
		self.onMediaTypeChange();
	});
	
	$(document).ready(function(){
		self.views = new Array();
		self.views['article'] = new View({url: 'js/templates/article.ejs'});
		self.views['book'] = new View({url: 'js/templates/book.ejs'});
		self.views['url'] = new View({url: 'js/templates/url.ejs'});
		$("[name='media-type']").change(function(){
			self.onMediaTypeChange.call(this);	
		});
	});
}
