function UnitSystemBrowser(controller, item)
{
	this.controller = "unit";
	Browser.call(this);
	this.url = "unit/systemslist.do"
	this.itemElement = "unitcategory";
	this.itemTemplate = 'js/templates/unit-system.ejs';
	this.pager.selectors.pager = ".pgr";
	
	this.doAddnewdialog = function()
	{
		$("#addnew").removeClass("invisible");
	}
	
	this.doCancel = function()
	{
		$("#addnew").addClass("invisible");
	}
	
	this.onAddnewsystemSuccess = function()
	{
		this.request(true);
		$("#addnew").addClass("invisible");
	}
}