var MMP = function() {};

MMP.isClassification = false;
MMP.labels = new Object();
			
MMP.decipherLabel = function(label) {
	label = label.toLowerCase();
	if (typeof MMP.labels[label] != 'undefined')
		return MMP.labels[label];
	if (label.indexOf("non") == 0)
		return MMP.labels[label] = -1;
	if (label.indexOf("inhibit") == 0)
		return MMP.labels[label] = 1;
	if (label.indexOf("in") == 0)
		return MMP.labels[label] = -1;
	if (label.indexOf("active") == 0)
		return MMP.labels[label] = 1;
	return MMP.labels[label] = 0;
}

function colorPairs()
{
	$(".pair TR.values").each(function(){
		var vals = [$(this).find("TD").eq(0).html(), $(this).find("TD").eq(1).html()];
		vals = vals.map(MMP.decipherLabel);
		if (vals[0] > vals[1])
			$(this).closest(".pair").addClass("increased");
		else if (vals[0] < vals[1])
			$(this).closest(".pair").addClass("decreased");
	});
}

MMP.drawHistogram  = function(query) {
	new QSPR.Ajax().send({
		url: "matchedpairs/getHistogram.do?" + query,
		success: function(response) {
			console.log(response);
			var dt = [];
			for (var i = 0; i < response.bins.length; i++)
				dt.push([response.bins[i], response.frequencies[i]]);
			var plot = new Plot();
			plot.addData(dt);
			plot.currentData.lines = {show: false};
			plot.currentData.points = {show: false};
			plot.currentData.bars = {show: true, barWidth: response.bins[1] - response.bins[0]};
			
			plot.drawLine(0, plot.minY, 0, plot.maxY);
			$(".histogram").removeClass("invisible");
			plot.render("#histogram");
			$("#xAxis").html($("tr[property].selected SPAN").html());
			
		}
	});
}

