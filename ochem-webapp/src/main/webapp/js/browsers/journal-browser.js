function JournalBrowser()
{
	Browser.call(this);
	
	this.url = "journal/list.do";
	this.itemElement = "journal";
	this.controller = "journal";
	this.actionURL = "journal/action.do";
	this.pager.selectors.pager = ".pgr";
	this.itemsPerPage = 10;
	
	this.beforeDelete = function()
	{
		return window.confirm("Do you really want to delete this journal");
	}
	
	this.onDeleteSuccess = function()
	{
		this.deleteRecord();
	}
	
	this.draw = function(entity)
	{
		return new View({url: 'js/templates/journal.ejs'}).render(entity);
	}
	
}
