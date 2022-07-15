include_js("js/commons/actionable.js?ver=1.7.5");

function MoleculeNameBrowser()
{
	var self = this;
	
	this.controller = "eprecord";
	
	UnifiedBrowser.call(this);
	this.itemTemplate = "js/templates/molsynonym.ejs";
	this.url = "eprecord/getsynonyms.do";
	
	this.itemElement = "moleculename";
	this.container = "MolNameBrowser";
	this.scope = this.filters.scope = "#MolNames";
	this.filters.firstTime = false;
	this.minusCounter = -1;
	
	
	this.getAdditionalData = function() //For the browser
	{
		return "&mol-id="+this.parent.oldId;
	}
	
	this.reset = function()
	{
		var container = $('#'+this.container);
		container.html("Information out of date. Please <a action='Update'>refresh</a>");
		this.attachActions(container);
		this.parent.tabView = new YAHOO.widget.TabView("demo");
	}	
	
	this.doUpdate = function()
	{
		return this.parent.callAction("Update");
	}
	
}

function MoleculeDepictionBrowser()
{
	var self = this;
	this.controller = "molecule";
	Browser.call(this);
	this.itemClass = "";
	this.itemElement = "molecule";
	this.container = "MolDepictions";
	this.scope = this.filters.scope = "MolDepictions";
	this.filters.firstTime = false;
	
	var view = false;
	
	this.draw = function(entity)
	{
		if (!view)
			 view = new View({url: 'js/templates/moleculedepiction.ejs'});
		var varentity = {data: entity};
		return view.render(varentity);
	}
	
	this.getAdditionalData = function()
	{
		return "&browser=depictions&mol-id="+this.parent.oldId
	}	
	
	this.reset = function()
	{
		var container = $('#'+this.container);
		container.html("Information out of date. Please <a action='Update'>refresh</a>");
		this.attachActions(container);
		this.parent.tabView = new YAHOO.widget.TabView("demo");
	}
	
	this.doUpdate = function()
	{
		return this.parent.callAction("Update");
	}
	
	this.doGet = function()
	{
		if (this.parent.oldId != -1)
			return this.parent.callAction("Get");
	}
	
	this.setPosition = function(entity)
	{		
		var tmpId = $(entity).attr("id");
		if (tmpId != undefined)
		{
			this.parent.oldId = tmpId;	
		}
	}
	
	this.onItemDrawn = function()
	{
		this.currentBlock.addClass("floating").attr("style","");
	}
	
}

function MoleculePropertyBrowser()
{
	var self = this;
	this.controller = "epbrowser";
	Browser.call(this);	
	this.container = "MolPropertiesBrowser";
	this.itemElement = "exp-property";
	this.scope = this.filters.scope = this.pager.selectors.scope = "#MolProperties";
	this.filters.firstTime = false;
	
	var view = false;
	
	this.draw = function(entity)
	{
		if (!view)
			view = new View({url: 'js/templates/moleculeproperty.ejs'});
		var varentity = {data: entity};
		return view.render(varentity);
	}
	
	this.getAdditionalData = function()
	{
		return "&similarmol="+this.parent.oldId
	}	
	
	this.reset = function()
	{
		var container = $('#'+this.container);
		container.html("Information out of date. Please <a action='Update'>refresh</a>");
		this.attachActions(container);
		this.parent.tabView = new YAHOO.widget.TabView("demo");
	}
	
	this.doUpdate = function()
	{
		return this.parent.callAction("Update");
	}
	
	this.doGet = function()
	{
		return this.parent.callAction("Get");
	}
	
	this.setPosition = function(entity)
	{		
		var tmpId = $(entity).attr("id");
		if (tmpId != undefined)
		{
			this.parent.oldId = tmpId;	
		}
	}
}

function MoleculeEditor()
{
	var self = this;
	var moleditor = this;
	this.editor = new JSME();
	var intervalID = 0;
	var t;
	var entity;
	var molData;
	var oldSmiles = '';
	this.oldId = getParams["id"];
	if ((this.oldId == undefined) || (!this.oldId))
		this.oldId = -1;
		
	this.scope = ".mainactions";	
	Actionable.call(this);
	
	this.mnBrowser = new MoleculeNameBrowser();
	this.mnBrowser.parent = this;
	
	this.mdBrowser = new MoleculeDepictionBrowser();
	this.mdBrowser.parent = this;

	this.mpBrowser = new MoleculePropertyBrowser();
	this.mpBrowser.parent = this;
		
	
	
	this.initialize = function()
	{
		this.mnBrowser.initialize(true);
		this.mdBrowser.initialize(true);
		this.mpBrowser.initialize(true);
	}

	this.checkChange = function()
	{
		try{
			if (oldSmiles != ''+self.editor.getMolecule())
			{
				clearInterval(moleditor.intervalID);
				self.mnBrowser.reset();
				self.mdBrowser.reset();
				self.mpBrowser.reset();
			}
		}catch(err)
		{
			//
		}	
	}
	
	this.getActionQuery = function(action)
	{
		if (action == "Get")
			return "mol-id="+this.oldId;
		else if (action == "Submit" || action=="Update")
		{
			var data = "";
			try{
				var molData = self.editor.getMolecule();
				if (molData){
					data = "data=" + molData;
				}else{
					window.alert("Please check the compound structure. May be it's not correct")
					data = "mol-id="+this.oldId;
				}
				
			}catch(err)
			{
				data = "mol-id="+this.oldId;
				throw err;
			}
			return data;
		}
		else if (action == "Smiles")
			return "data="+URLEncode(document.molViewer.smiles.value);
		else if (action == "Searcher")
			return "Search="+escape(document.molViewer.namesearch.value);
	}
	
	this.getData = function()
	{
		//if(document.JME.isActive())
		//{
			//applet is ready
			try{
				//document.JME.readMolFile(molData);
				self.editor.setMolecule(molData);
				oldSmiles = ''+self.editor.getMolecule();
				
				clearInterval(moleditor.intervalID);
				moleditor.intervalID = setInterval(moleditor.checkChange, 300);	// Start the change-checking
			}catch(err)
			{
				$("JME").addClass("invisible");
				$("#javaexc").html('<img src="depiction.jsp?id='+entity.id+'" width="300px" height="300px"/>');
			}
			clearTimeout(t);
		//}
		//else{
		//	t = setTimeout(moleditor.getData(), 1000);
		//}
	}
	
	this.onGetSuccess = function(xml) //Get single molecule
	{
		entity = xml.molecule;
		
		if (entity != undefined)
		{
			molData = URLDecode(entity.encodedData);
		} else 
		{
			molData = "";
		}
		
		if (molData != "")
		{
			//if(document.applets[0]!= null)
			//{
				moleditor.getData();
			//}
			//else
			//{
			//	$("JME").addClass("invisible");
			//	$("#javaexc").html('<img src="depiction.jsp?id='+entity.id+'" width="300px" height="300px"/>');
			//}
		}
		
		if(entity != undefined && entity.id != undefined)
			this.oldId = entity.id;
		else
			this.oldId = -1;
		try{
			
		}catch(err)
		{
			//
		}
		if ((this.oldId) && (this.oldId != undefined) && (this.oldId != -1))
		{
			this.mnBrowser.request();
			this.mdBrowser.request();
			this.mpBrowser.request();
		}	
								
	}
	
	this.onSubmitSuccess = function(xml) //Save a molecule
	{
		var id = xml.molecule.id;
		
		if (window.callback != undefined)
			window.callback(id);
		
		if (window.callback_extended != undefined)
			window.callback_extended(xml.molecule);
		
		window.closeTab();		
	}
	
	this.onUpdateSuccess = function(xml) //Save a molecule and reload it's info
	{
		this.onGetSuccess(xml);					
	}

	this.onSmilesSuccess = function(xml) //Reload from smiles
	{
		this.onGetSuccess(xml);					
	}
	
	this.doSearcher = function()
	{
		// URLEncode (ie encodeURIComponent()) is not enough, because prime (') in names is not escaped causing broken url 
		var searchstring = escape((document.molViewer.namesearch.value).trim()); //URLEncode(document.molViewer.namesearch.value);

		var molid = "";
		var map1id = ""; 
		var map2id = "";
		
		if (entity != undefined){
			molid = entity.id; // mol id
			map1id = entity.mapping.mapping1.id; 
			map2id = entity.mapping.id; // mapping2 id
		} 
		
		var propWin = openTab("NCBI Search", webRoot + "ncbisearch/searchBrowser.do?render-mode=popup&search-string=" + searchstring
							+ "&existmol-id=" + molid
							+ "&mapping1-id=" + map1id
							+ "&mapping2-id=" + map2id);
		
		propWin.moveBy(50, 50);
		propWin.callback = function(molecule)
		{
			
			document.molViewer.namesearch.value =  molecule.searchedBy.name;
			
			propWin.closeTab();
			//this.checkChange();
			var a = new Object();
			a.molecule = molecule;
			self.onGetSuccess(a);					
			
			setTimeout('editor.mnBrowser.request();', 1);
			setTimeout('editor.mdBrowser.request();', 1);
			setTimeout('editor.mpBrowser.request();', 1);
		}
	}
	
	function updateName(link)
	{
		var oldname = $(link).attr('name');
		var name = window.prompt("Enter new Molecule name:","");
			if (name != undefined && name !='')
			{
				ajax.send
				(
					{ 
						url: 'molecule/action.do', 
						data: 'mol-id='+oldId+'&oldname='+oldname+'&newname='+name,
						success: function(xml)
						{
							drawUpdated(xml);
						}
					}
				);	
			}
			else
			{
				window.alert('Please enter valid name');
			}
	}
	
	function addName()
	{
		var name = window.prompt("Enter new Molecule name:","");
		if (name != undefined && name !='')
		{
			$.ajax
			(
				{ 
					url: 'molecule/action.do', 
					data: 'mol-id='+oldId+'&name='+name,
					success: function(xml)
					{
						drawUpdated(xml);
					}
				}
			);	
		}
		else
		{
			window.alert('Please enter valid name');
		}
	}
	
	this.iframeLoaded = function()
	{	
		var data = $("#iframe").contents();
		if ((data.text() != "null")&&(data.text() != ""))
		{
			json = eval("("+data.text()+")");
			this.onGetSuccess(json); //????	
			document.molViewer.file.value = "";
		}
	}
	
	this.beforeGet = function() {
		console.log(this.oldId);
		return this.oldId != "-1";
	}
	
	function deleteName(link)
	{
		var name = $(link).attr('name');
		$.ajax
		(
			{ 
				url: 'molecule/action.do', 
				data: 'mol-id='+oldId+'&delname='+name,
				dataType: 'xml', 
				success: function(xml)
				{
					drawUpdated(xml);
				}
			}
		);		
	}
}

