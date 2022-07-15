include.plugins('view');
DataCubeReport = function(url) {
	var ajax = new QSPR.Ajax();
	var self = this;
	var tbody = $(".report TBODY");
	var actionable = new Actionable();
	actionable.scope = $(".groupings");
	this.dataDrawn = new Event();
	this.filters = new Filters();
	this.filters.setFromUrl();
	
	this.load = function() {
		
		var subkontos = [];
		var drawLine = function(node, depth, maxDepth, queryArray) {
			var rowDepth = maxDepth == 1 ? depth + 1 : depth;
			var tr = new $("<tr class='depth-"+rowDepth+"'/>");

			var td = $("<td/>");
			if (depth == 0)
				td.append("Totals");
			else
				td.append(self.drawNode(node));
			tr.append(td);
			if (depth > 0)
				queryArray.push("" + subkontos[depth - 1].name + "=" + node.id);
			self.drawMetrics(tr, node.metrics, queryArray.join("&").toLowerCase() + "&" + self.filters.getQueryString(), node);
			tbody.append(tr);
			var children = array(node.child);
			for (var i = 0, len = children.length; i < len; i++)
				drawLine(children[i], depth + 1, maxDepth, queryArray);
			queryArray.pop();
		}
		
		var getGroupings = function() {
			var q = [];
			$(".groupings").find("input:checked").each(function() {q.push($(this).attr("name"));});
			return q.join(",");
		}
		
		$(".progress").removeClass("invisible");
		tbody.html("");
		$(".report THEAD TR TH.subkonto").remove();
		ajax.send({
			url: url,
			data: "groupings=" + getGroupings() + "&" + self.filters.getQueryString(),
			success: function(response) {
				subkontos = array(response.awaitingDataReport.dataTree.subkontos);
				$(".report THEAD TR").prepend("<th class='subkonto'>" + subkontos.map(function(o) {
					return o.name;
				}).join(" / ") + "</th>");
				
				if (response.others)
					self.responseOthers = response.others;
				
				drawLine(response.awaitingDataReport.dataTree, 0, subkontos.length, []);
				$(document).trigger("DOM_updated", $(".report"));
				self.dataDrawn.fire();
			},
			after: function() {
				$(".progress").addClass("invisible");
			}
		});
	}
	
	var groupingsUpdated = function() {
		$(".groupings DIV.grouping A[action]").removeClass("invisible");
		$(".groupings DIV.grouping:first A[action=up]").addClass("invisible");
		$(".groupings DIV.grouping:last A[action=down]").addClass("invisible");
		self.load();
	}
	
	actionable.doUp = function(link) {
		var div = $(link).closest("DIV");
		div.insertBefore(div.prev());
		groupingsUpdated();
	}
	
	actionable.doDown = function(link) {
		var div = $(link).closest("DIV");
		div.insertAfter(div.next());
		groupingsUpdated();
	}
	
	$(function(){
		$(".groupings INPUT").change(groupingsUpdated);
		$(".groupings DIV.grouping").append("<a action='up' title='Move the grouping up'><img src='/img/icons/up.png'/></a><a action='down' title='Move the grouping down'><img src='/img/icons/down.png'/></a>");
		actionable.attachActions($(".groupings"));
		groupingsUpdated();
	});
}

var Event = function(callback) {
	var self = this;
	this.listeners = new Array();
	
	this.fire = function(val) {
		this.listeners.map(function(listener){
			listener(val);
		});
	}
	
	this.register = function(callback) {
		self.listeners.push(callback);
	}
}

var templatesMap = new Object();
var getTemplate = function(id) {
	if (!templatesMap[id])
		return templatesMap[id] = new View({element: id});
	else
		return templatesMap[id];
}

