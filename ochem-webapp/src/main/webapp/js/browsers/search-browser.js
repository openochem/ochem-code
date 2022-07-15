function SearchBrowser()
{
	this.controller = "ncbisearch";
	Browser.call(this);
	this.itemElement = "search-molecule";	
	this.itemTemplate = "js/templates/searchmol.ejs";
	var self = this;
	
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

	
	this.onItemDrawn = function()
	{
		this.currentBlock.find('a').tooltip({showURL: false});
	}
	
	this.filters.onFilterChange = function()
	{
	}
	
	this.onSearchSuccess = function(xml)
	{
		
	}
}