var
	dsAuthorsAC = 
		new YAHOO.widget.DS_XHR("author/list.do", ["author", "printed-name", "id"]);
	dsAuthorsAC.scriptQueryParam = "query"; 
	dsAuthorsAC.responseType = YAHOO.widget.DS_XHR.TYPE_XML; 
	dsAuthorsAC.maxCacheEntries = 0; 
	dsAuthorsAC.scriptQueryAppend = "out=xml";
	
var
	dsSourcesAC = 
		new YAHOO.widget.DS_XHR("article/list.do", ["article", "title", "id"]);
	dsSourcesAC.scriptQueryParam = "query"; 
	dsSourcesAC.responseType = YAHOO.widget.DS_XHR.TYPE_XML; 
	dsSourcesAC.maxCacheEntries = 0; 
	dsSourcesAC.scriptQueryAppend = "out=xml";
	
var
	dsJournalsAC = 
		new YAHOO.widget.DS_XHR("journal/list.do", ["journal", "title", "id"]);
	dsJournalsAC.scriptQueryParam = "query"; 
	dsJournalsAC.responseType = YAHOO.widget.DS_XHR.TYPE_XML; 
	dsJournalsAC.maxCacheEntries = 0; 
	dsJournalsAC.scriptQueryAppend = "out=xml";
	
var
	dsBookAC = 
		new YAHOO.widget.DS_XHR("journal/list.do", ["journal", "title", "id"]);
	dsBookAC.scriptQueryParam = "query"; 
	dsBookAC.responseType = YAHOO.widget.DS_XHR.TYPE_XML; 
	dsBookAC.maxCacheEntries = 0; 
	dsBookAC.scriptQueryAppend = "out=xml&book=1";
	
var
	dsUnitsAC = 
		new YAHOO.widget.DS_XHR("unit/list.do", ["unit", "name", "id"]);
	dsUnitsAC.scriptQueryParam = "query"; 
	dsUnitsAC.responseType = YAHOO.widget.DS_XHR.TYPE_XML; 
	dsUnitsAC.maxCacheEntries = 0; 
	dsUnitsAC.scriptQueryAppend = "out=xml";
	
var
	dsProperties = new YAHOO.util.DataSource("properties/list.do?");
	dsProperties.responseType = YAHOO.util.DataSource.TYPE_XML;
	dsProperties.connXhrMode = "queueRequests";
	dsProperties.responseSchema = 
	{
		resultList: "model.list",
		fields:
		[
			"id", "name"
		]
	};
	
function getProperties(query, callBack)
{
	dsProperties.sendRequest('out=xml&query='+query+'');
}	

var
	dsPropertiesAC = 
		new YAHOO.widget.DS_XHR("properties/list.do", ["property", "name", "id"]);
	dsPropertiesAC.scriptQueryParam = "query"; 
	dsPropertiesAC.responseType = YAHOO.widget.DS_XHR.TYPE_XML; 
	dsPropertiesAC.maxCacheEntries = 0; 
	dsPropertiesAC.scriptQueryAppend = "out=xml";
	
var
	dsConditionsAC = 
		new YAHOO.widget.DS_XHR("conditions.list", ["condition", "name", "id"]);
	dsConditionsAC.scriptQueryParam = "query"; 
	dsConditionsAC.responseType = YAHOO.widget.DS_XHR.TYPE_XML; 
	dsConditionsAC.maxCacheEntries = 0; 
	dsConditionsAC.scriptQueryAppend = "out=xml";