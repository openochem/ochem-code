function MessageBrowser()
{
	this.controller = "mail";
	Browser.call(this);
	this.itemElement = "message";
	this.itemTemplate = 'js/templates/message.ejs';
	
	this.onReadSuccess = function(xml)
	{
		this.request();	
	}
	
	this.onAllreadSuccess = function(xml)
	{
		this.pageNum = 1;
		this.request();	
	}
	
	this.doShowall = function(entity)
	{
		this.filters.values.type = "showall";
		this.request(true);
	}
	
	this.doUnread = function(entity)
	{
		this.filters.values.type = "unread";
		this.request(true);
	}
	
	this.doSent = function(entity)
	{
		$("#sent").addClass("fontweight");
		$(".upper-command-panel").addClass("invisible");
		$("#reply").addClass("invisible");
		$("#inbox").removeClass("fontweight");
		
		this.filters.values.type = "sent";
		this.request(true);
	}
	
	this.doInbox = function(entity)
	{
		$("#inbox").addClass("fontweight");
		$("#sent").removeClass("fontweight");
		$(".upper-command-panel").removeClass("invisible");
		
		this.filters.values.type = "inbox";
		this.request(true);
	}
}
