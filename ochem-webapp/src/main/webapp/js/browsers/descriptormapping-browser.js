function DescriptorMappingBrowser()
{
	this.controller = "descriptormapper";
	this.url = "descriptormapper/list.do";
	this.actionURL = "descriptormapper/action.do";
	this.itemTemplate = "js/templates/descriptor-mapping.ejs";
	Browser.call(this);
	this.pager.selectors.pager = ".pgr";
	this.itemElement = "descriptorMapping";
	
	this.doEditmapping = function(link)
	{
		$(this.currentBlock).find(".view").addClass("invisible");
		$(this.currentBlock).find(".edit").removeClass("invisible");
	}
	
	this.onSavemappingSuccess = function(item)
	{
		item = item.descriptorMapping;
		var options = {"replaceId":item.id};
		this.drawFromJSON(item, options);
	}
}
 
include.plugins('view');
var sampleBrowser = new DescriptorMappingBrowser();
$(document).ready(function(){
	sampleBrowser.initialize();
});