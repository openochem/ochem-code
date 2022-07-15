include.plugins('view');


$(document).ready(function()
{
	var ajax = new QSPR.Ajax();
	var nb = new NeighboursBrowser();
	
	$("#main-profile").html("<img src='img/roller_transparent.gif'/>");
	
	if (ep_id != undefined)
	{
		ajax.send({
			url: "modeldot/list.do?mm_id=" + mm_id + "&ep_id=" + ep_id + "&statnum=1&model_id=" + model_id,
			success: function (response) 
			{
				var arr = array(response.list["exp-property"]);
				var profile = new View({url: "js/templates/dotprofile-neighbor.ejs"}).render(arr[0]);
				$("#main-profile").html(profile);
				$("#main-profile").find(".neighbors-link").remove();
			}
		});
	} else
	if (row_num != undefined && task_num != undefined)	
	{
		ajax.send({
			url: "modelapplier/reslist.do?mm_id=" + mm_id + "&row_num=" + row_num + "&task_num=" + task_num + "statnum=1&model_id=" + model_id,
			success: function (response) 
			{
				var arr = array(response.list["predictionRow"]);
				var profile = new View({url: "js/templates/models/modelresult.ejs"}).render(arr[0]);
				$("#main-profile").html(profile);
				$("#main-profile").find(".neighbors-link").remove();
			}
		})		
	}
	
	nb.initialize(false);
});		

function NeighboursBrowser()
{
	var self = this;
	
	this.controller = "modelneighbours";
	this.scope  = "#neighbours";
	
	UnifiedBrowser.call(this);
	
	this.pager.useHistory = false;
	
	this.filters.scope = this.scope;
	this.container = "nearest-neighbours-browser";
	this.itemTemplate = "js/templates/dotprofile-neighbor.ejs";
	this.itemElement = "exp-property";
	this.filters.firstTime = true;
	this.minusCounter = -1;
	
	this.filters.setFromUrl();
	
	this.onItemDrawn = function()
	{
		this.currentBlock.find('a[title]').tooltip({showURL: false});
	}
}

var labelMap = new Array();;
labelMap["correlation"] = "Correl";
labelMap["rank-correlation"] = "RankCorrel";
labelMap["euclidean"] = "Distance";
labelMap["structural-similarity"] = "Similarity";
