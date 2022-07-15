function plot1clicked(setNum, pointNum)
{
	var epId = this.dataArray[setNum].data[pointNum][2];
	openTab("Model's compound", webRoot+"modeldot/show.do?render-mode=popup"+mid+"&statnum=0&ep_id=" + epId);
}

function plot2clicked(setNum, pointNum)
{
	var epId = this.dataArray[setNum].data[pointNum][2];
	openTab("Model's compound", webRoot+"modeldot/show.do?render-mode=popup"+mid+"&statnum=1&ep_id=" + epId);
}

function showEstimatedStatistics() 
{

	var appendStats = function(row, stats) {
		console.log(stats);
		if (qualitativeProperty) 
		{
			var cs = stats.classificationSummary;
			row.append("<td>" + cs.accuracyTotal["formatted-value"] + "% &plusmn; " + cs.accuracyTotal["formatted-std"] + "</td>");
			row.append("<td>" + cs.accuracyBalanced["formatted-value"] + "% &plusmn; " + cs.accuracyBalanced["formatted-std"] + "</td>");
			row.append("<td>" + cs.mcc["formatted-value"] + " &plusmn; " + cs.mcc["formatted-std"] + "</td>");
		}
		else
		{
			row.append("<td>" + stats.r2 + " &plusmn; " + stats.r2std + "</td>");
			row.append("<td>" + stats.q2 + " &plusmn; " + stats.q2std + "</td>");
			row.append("<td>" + stats.rmse + " &plusmn; " + stats.rmsestd + "</td>");
			row.append("<td>" + stats.mae + " &plusmn; " + stats.maestd + "</td>");
		}
	}

	var appendTwoStats = function(table, stats) 
	{
		table.append("<tr><td colspan='5'>" + stats.basket.name + "</td></tr>");
		console.log(stats);
		var row = $("<tr/>");
		row.append("<td>Estimated</td>");
		appendStats(row, stats.estimatedStatistics);
		table.append(row);

		row = $("<tr/>");
		row.append("<td>Actual</td>");
		appendStats(row, stats.statistics);
		table.append(row);
	}

	var ajax = new QSPR.Ajax();
	$("#estimated-statistics").html("<img src='img/roller_transparent.gif'/>");
	ajax.send({
		url: "model/getEstimatedStatistics.do?" + mid,
		success: function (response) {
			var sets = array(response.modelProfileData.setStats);
			var table = qualitativeProperty ?
					$("<table class='tiny'><tr><td></td><td>Accuracy</td><td>Balanced accuracy</td><td>MCC</td></tr></table>") : 
						$("<table class='tiny'><tr><td></td><td>R2</td><td>Q2</td><td>RMSE</td><td>MAE</td></tr></table>");

					$("#estimated-statistics").html(table);
					for (var i = 0; i < sets.length; i++) {
						appendTwoStats(table, sets[i]);
					}
		}
	});
}

function errorsShow(recalculated)
{
	var suffix = recalculated ? "&recalculated=1" : "";
	openTab("Errors of the model", webRoot+"modeldot/errorshow.do?render-mode=popup" + mid + suffix);
}

function replaceOriginalModel()
{
	if (!window.confirm("The original model will be overwritten by the recalculated model. This action cannot be undone. Are you sure you want to continue?"))
		return false;
	location.href = 'model/profile.do?render-mode=popup&storeRecalculatedModel=1' + mid;
}

function discardRecalculatedModel()
{
	if (!window.confirm("Are you sure you want to discard the recalculated model?"))
		return false;
	location.href = 'model/profile.do?render-mode=popup&deleteRecalculatedModel=1' + mid;
}

function excludeDuplicates()
{
	location.href = 'model/profile.do?render-mode=popup&exclude-duplicates=1' + mid;	
}

function restorePoints(recalculated)
{
	location.href = 'model/profile.do?render-mode=popup&restore-basket=' + (recalculated ? "recalculated" : "original") + mid;	
}

function deleteValidationSet()
{
	location.href = 'model/profile.do?render-mode=popup&delete-validation-set='+ valSetId + mid;
}

Action = {exclude: "exclude", include: "include"};
var SelectionHandler = function()
{
	var self = this;
	$(this.selector + "-selection").get(0).plot = this;
	this.ajax = new QSPR.Ajax("modeldot/action.do");

	this.onSelectionChanged = function(value)
	{
		var num = Object.size(this.selectedPoints);
		var sn = Math.round (10* num / (value * this.dataArray[0].data.length), 2)/10;

		$(this.selector + "-selection").html(
				'' + num + ' records selected; S/N = ' + sn + ' ' +
				'<a href="#" onclick="SelectionHandler.getPlot(this).unselectAll(); return false;">[unselect]</a>' + 
				'<a href="#" title="Exclude the selected points from the training dataset" onclick="SelectionHandler.getPlot(this).excludeSelected(Action.exclude); return false;">[exclude]</a>' + 
				'<a href="#" title="Include the selected points back to the training dataset" onclick="SelectionHandler.getPlot(this).excludeSelected(Action.include); return false;">[include]</a>'+
				'<a href="#" title="Browse the selected points in the compound property browser" onclick="SelectionHandler.getPlot(this).browseSelected(); return false;">[browse]</a>'
		);
		$(this.selector + "-selection").setClass("invisible", Object.size(this.selectedPoints) == 0);
	}

	this.excludeSelected = function(action)
	{
		this.ajax.send({
			data: "model_id=" + modelId + "&action=" + action + "&id=" + this.getSelectedIds(),
			success: function(reply)
			{
				window.alert(reply.message.message + "\nThe page will now be reloaded.");
				window.location.reload();
			}
		});
	}

	this.getSelectedIds = function()
	{
		var ids = new Array();
		for (var hash in this.selectedPoints)
		{
			var p = this.selectedPoints[hash];
			ids.push(this.dataArray[p.seriesNum].data[p.pointNum][2]);
		} 

		return ids.join(",");
	}

	this.browseSelected = function()
	{
		var selector = this.selector + "-selection";  
		$(selector).children().last().remove();
		var roller = "<img src=\"img/roller_small.gif\"/>";
		$(selector).append("<i>"+roller+"preparing selection basket...</i>");
		var basketAjax = new QSPR.Ajax("epbrowser/action.do");
		basketAjax.send({
			data: "action=addbasket&clearbasket=1&basket-name=Model profile selected records&id=" + this.getSelectedIds(),
			success: function(reply)
			{
				$(selector).children().last().remove();
				$(selector).append('<a href="#" title="Browse the selected points in the compound property browser" onclick="SelectionHandler.getPlot(this).browseSelected(); return false;">[browse]</a>');
				openTab("Selected records", "epbrowser/show.do?basket-select=Model profile selected records");
			}
		});
	}
}

SelectionHandler.getPlot = function(link)
{
	return $(link).parent().get(0).plot;
}

//A fancy plugin to synchronize highlights on different plots
var PlotSynchronizer = function()
{
	var self = this;
	this.friendlyPlots = new Array();

	$(this.selector).bind("plothover", function(e, pos, item)
			{
		if (item)
		{
			var data = self.dataArray[item.seriesIndex].data[item.dataIndex];
			for (var plot in self.friendlyPlots)
			{
				self.friendlyPlots[plot].clearAutoHighlights("plothover");
				var dublicate = self.friendlyPlots[plot].pointHash['a' + data[2]];
				if (dublicate)
					self.friendlyPlots[plot].flot.highlight(dublicate.seriesNum, dublicate.pointNum, "plothover");
			}
		}

		return false;
			});

	this.createPointHash = function()
	{
		this.pointHash = new Array();
		for (var i = 0; i < this.dataArray.length; i++)
			for (var k = 0; k < this.dataArray[i].data.length; k++)
				this.pointHash['a' + this.dataArray[i].data[k][2]] = {seriesNum: i, pointNum: k};
	}
}

function selectOutliers(plot, adIntervals, adErrors, quantile, value)
{
	plot.unselectAll();

	for (var set = 0; set < plot.dataArray.length; set++)
	{	
		if (plot.dataArray[set].label && !plot.dataArray[set].label.startsWith("validation"))
			for (var i = 0; i < plot.dataArray[set].data.length; i++)
			{
				var interval = 0;
				while (interval < adIntervals.length && plot.dataArray[set].data[i][0] > adIntervals[interval])
					interval++;

				if (plot.dataArray[set].data[i][1] > adErrors[interval] * quantile)
					plot.selectPointLight(set, i);
			}
	}
	plot.onSelectionChanged(value);
}

function showSizes() 
{
	var ajax = new QSPR.Ajax();
	ajax.send({
		url: "model/getSizeSummary.do",
		data: mid,
		dataType: "json",
		success: function(response) {
			$("#size-summary").removeClass("invisible");
			$("#size-summary").html("Configuration: " + response.configuration + "<br/>Model: " + response.model + "<br/>Descriptors: " + response.descriptors + "<br/>Statistics: " + response.statistics);
		}
	});
}

