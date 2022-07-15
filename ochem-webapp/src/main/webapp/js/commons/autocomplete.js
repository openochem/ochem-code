
// Midnighter
// Autocomplete
// represents reusable class for 
// autocomplete UI block

function Autocomplete(id, ds, name)
{
	this.id = id;
	if (name == undefined)
		this.name = this.id;
	else
		this.name = name;
	this.preffix = "ac-"+id;
	this.datasource = ds;
	Autocomplete.paHash['_'+id] = this;
	this.queryMade = false;
	
	this.getKey = function()
	{
		return document.getElementById(this.preffix+"-key").value;	
	}
	
	this.setKey = function(key)
	{
		var elem = document.getElementById(this.preffix+"-key");
		elem.value = key;
		var oldValue = $(elem).attr('oldValue');
		$('#'+this.preffix+'-odiv').attr('value', key);
		if (oldValue != key)
			this.evtChange.fire(); //Experimental
	}
	
	this.getValue = function()
	{
		return document.getElementById(this.preffix+"-input").value;
	}
	
	this.setValue = function(value)
	{	
		document.getElementById(this.preffix+"-input").value = value;
		$('#'+this.preffix+'-odiv').attr('title', value);	
	}

	this.attachHandlers = function()
	{
		this.autocomplete = new YAHOO.widget.AutoComplete(
    				this.preffix+'-input',
    				this.preffix+'-container', 
    				this.datasource);
    	this.autocomplete.holder = this;
    	this.autocomplete.prehighlightClassName = "yui-ac-prehighlight";
    	this.autocomplete.typeAhead = false;
    	this.autocomplete.useShadow = true;
    	this.autocomplete.forceSelection = false;
    	
    	this.autocomplete.formatResult = 
    		function(oResultItem, sQuery) 
    		{
        		var sMarkup = oResultItem[0];
        		return (sMarkup);
    		};	
    		
    	this.autocomplete.itemSelectEvent.fire = 
    		this.onSelect;
    		
    	this.autocomplete.dataReturnEvent.fire =
    		function (oSelf, sQuery, aResults)
    		{
    			oSelf.holder.queryMade = true;    			
    		} 
    	
    	this.unmatched = 
    		function(ac)
    		{
    			ac.setValue("");
    		}
    		
    	this.autocomplete.unmatchedItemSelectEvent.fire = 
    		function(oSelf)
    		{
    			oSelf.holder.setKey("");
    			if (oSelf.holder.queryMade)
    				oSelf.holder.unmatched(oSelf.holder);
    		}
    		
    	$('#'+this.preffix+'-odiv').get(0).autocomplete = this;
    		
    	promptInputHandler($(this.preffix+"-input"));
	}
		
	this.getJQuery = function(_class)
	{
	    if (!_class)
	    	_class = " ";
	    	
		return this.jQuery = 
			$('<div class="autocomplete'+_class+'" class="nomargin nopadding" id="'+this.preffix+'-odiv" filter="1" name="'+this.id+'"/>')
			.append('<input type="hidden" style="border: none; visibility: hidden; display: none;" name="'+this.name+'-key" id="'+this.preffix+'-key"/>')
			.append('<input class="autocomplete-input'+_class+' nomargin" prompt="Start typing to select" name="'+this.name+'-title" id="'+this.preffix+'-input" type="text"/>')
			.append('<div class="nomargin nopadding" id="'+this.preffix+'-container"/>');
	}
	
	this.onSelect = function(oSelf, elItem, oData)
	{
		var holder = oSelf.holder;
		if (holder.getKey() != oData[1])
		{
    		holder.setKey(oData[1]);
    		holder.evtChange.fire(oSelf.holder.id);
    	}	
	}
	
	this.evtChange = new YAHOO.util.CustomEvent("changed");
}

// "Static" declarations
Autocomplete.getInstance = function(id)
{
	return Autocomplete.paHash['_'+id];
}

Autocomplete.updateFromHTML = function()
{
	for (key in Autocomplete.paHash)
	{
		var ac = Autocomplete.paHash[key];
		var div = $('div#'+ac.preffix+'-odiv');
		ac.setKey(div.attr('value'));
		if (div.attr('title') != undefined)
			ac.setValue(div.attr('title'));
	}
}


Autocomplete.paHash = new Object();
