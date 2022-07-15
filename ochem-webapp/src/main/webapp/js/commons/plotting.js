function Plot()
{
	var self = this;
	
	this.dataArray = new Array();
	this.minX = 100000;
	this.maxX = -10000;
	this.minY = 10000;
	this.maxY = -10000;
	this.defaultChartColors = new Array("#C00", "#0C0", "#00C");
	this.series = new Object();
	this.defaultColor = "#CCC";
	this.selectedPoints = new Object();
	this.shiftKey = false;
	this.altKey = false;
	this.normalizeAspectRatio = false;
	this.enlarged = false;
	
	this.adaptSize = true;
	
	
	this.bindEvents = function(plot, eventHolder)
   	{
		eventHolder.mousedown(function(e){
			self.shiftKey = e.shiftKey;
			self.altKey = e.altKey;
		});
   	}
	
	this.options = 	
	{
   		xaxis: {autoscaleMargin: 0.02},
   		yaxis: {autoscaleMargin: 0.02},
   		grid:  {clickable: true, hoverable: true},
   		selection: { mode: "xy" },
   		legend: 
   		{
   			labelFormatter: function(label, series) {
		    	return "<span series='" + series.num + "'>" + (series.highlight ? '<b>' + label + '</b>' : label) + "</span>";
  			}
   		},
   		hooks: { bindEvents: [this.bindEvents] }
   	};
   	
   	
   	
   	this.getData = function(x1, y1, x2, y2)
   	{
   		for (var i = 0; i < this.dataArray; i++)
   			for (var k = this.dataArray.data.length - 1; k >= 0; k--)
   				if (this.dataArray.data[k][0] < x1 || this.dataArray.data[k][0] > x2 || this.dataArray.data[k][1] < y1 || this.dataArray.data[k][1] > y2)
   					this.dataArray.data.splice(k, 1);
   	}
   	
   	this.addPoint = function(seriesId, point)
   	{
   		if (!this.series[seriesId])
   			this.dataArray.push(this.series[seriesId] = {
	   			color: this.defaultChartColors[this.dataArray.length],
	       		data: new Array(),
	       		points: { show: true, fill: true },
	       		clickable: true
		    });
   		this.series[seriesId].data.push(point);
   		this.correctMinMax(point);
   	}
   	
   	this.clear = function()
   	{
   		this.dataArray = new Array();
   	}
   	
   	this.setPoints = function(_points)
   	{
   		this.currentData.points = _points;
   	}
   	
   	this.setLabel = function(label)
   	{
   		this.currentData.title = this.currentData.label = label;
   		return this;
   	}
   	
	this.addData = function(_data, _color, _clickable)
	{
		if (!_color)
			_color = this.defaultChartColors[this.dataArray.length];
		this.dataArray.push(
			this.currentData = 
			{
	   			color: _color,
	       		data: _data,
	       		points: { show: true, fill: true },
	       		clickable: _clickable,
	       		num: this.dataArray.length
		    }
		);
		
		var index = this.dataArray.length-1;
		
		if (this.adaptSize)
			for (var i = 0; i < this.dataArray[index].data.length; i++)
				this.correctMinMax(this.dataArray[index].data[i]);
		
		this.dataRange = new Object();
		this.dataRange.minX = this.minX;
		this.dataRange.maxX = this.maxX;
		this.dataRange.maxY = this.maxY;
		this.dataRange.minY = this.minY;
			
	    return this;
	}
	
	this.correctMinMax = function(point)
	{
		if (point[0] > this.maxX)
			this.maxX = point[0];
			
		if (point[0] < this.minX)
			this.minX = point[0];

		if (point[1] > this.maxY)
			this.maxY = point[1];
			
		if (point[1] < this.minY)
			this.minY = point[1];	
	}
	
	this.onLegendClick = function(series)
	{
		
	}
	
	this.addBaseLine = function()
	{
		return this.drawLine(this.minX, this.minX, this.maxX, this.maxX);
	}
	
	this.fillArea = function(minX, minY, maxX, maxY)
	{
		var as = this.adaptSize;
		this.adaptSize = false;
		this.addData([[minX, maxY, minY]]);
		this.adaptSize = as;
		
		this.currentData.points = {show: false};
		this.currentData.lines = {show: false};
		this.currentData.bars = { 
	        show: true,
	        barWidth : maxX - minX,
	        fill: true,
            fillColor: { colors: [ "rgba(0, 100, 0, 0.06)", "rgba(0, 100, 0, 0.06)" ] }
	    }
		this.currentData.color= "rgba(0, 100, 0, 0.01)";
		this.currentData.clickable= false;
		this.currentData.hoverable= false;
	}
	
	this.drawLine = function(x1, y1, x2, y2)
	{
		var line = [[x1, y1], [x2, y2]];

		this.addData(line, this.defaultColor, false);
		this.setPoints({line: true});
		return this;	
	}
	
	this.verticalLine = function(x)
	{
		this.drawLine(x, this.minY, x, this.maxY);
	}
	
	this.abs = function(x)
	{
		if (x<0)
			return -x;
		else
			return x;
	}
	
	this.toggleFullScreen = function()
	{
		self.enlarged = !self.enlarged;
		if (self.enlarged)
		{
			self.originalSelector = self.selector;
			var bigBigDiv = $("<div/>");
			bigBigDiv.addClass("enlarged-plot");
			bigBigDiv.html("aaa");
			
			var wrapperDiv = $("<div/>");
			wrapperDiv.css({position: "relative"})
			var onclick = "Plot.getPlot('"+self.selector+"').toggleFullScreen(); return false;";
			wrapperDiv.html('<a onclick="'+onclick+'" href="javascript:void()" title="Toggle full screen mode" class="full-screen-link"><img src="img/icons/full-screen-out.gif"/></a>');
			wrapperDiv.prepend(bigBigDiv);
			wrapperDiv.addClass("enlarged-plot-wrapper");
			
			$("BODY").append(wrapperDiv);
			self.render(".enlarged-plot");
		}
		else
		{
			$(".enlarged-plot-wrapper").remove();
			self.render(self.originalSelector);
		}
	}
	
	this.onClick = function(e, pos)
	{
		var xdistance;
		var ydistance;
		
		if (self.altKey)
		{
			self.toggleFullScreen();
			return;
		}
		
		for (var k = self.dataArray.length-1; k>=0; k--)
		{
			if (!(self.dataArray[k].clickable == false))	
			for (var i = self.dataArray[k].data.length-1; i>=0;  i--)
			{
				xdistance = self.abs(self.dataArray[k].data[i][0]-pos.x)
				ydistance = self.abs(self.dataArray[k].data[i][1]-pos.y);
				if ((xdistance < (self.maxX-self.minX)/35) && (ydistance < (self.maxY-self.minY)/25))
				{
					if (self.shiftKey)
						self.selectPoint(k, i);
					else
						self.pointClicked.call(self, k, i, self.dataArray[k].data[i]);
					return;
				}	
			}
		}	
		
		//self.enlarge();
		//$.plot($(self.selector), self.dataArray, self.options);
		self.zoomTo(self.dataRange.minX, self.dataRange.maxX, self.dataRange.minY, self.dataRange.maxY, false);
	}
	
	this.selectPoint = function(seriesNum, pointNum)
	{
		var hash = "" + seriesNum + ":" + pointNum;
		if (this.selectedPoints[hash])
		{
			this.flot.unhighlight(seriesNum, pointNum);
			delete this.selectedPoints[hash];
		}
		else
		{
			this.flot.highlight(seriesNum, pointNum);
			var point = new Object();
			point.seriesNum = seriesNum;
			point.pointNum = pointNum;	
			this.selectedPoints[hash] = point;
		}
		
		this.onSelectionChanged();
	}
	
	this.selectPointLight = function(seriesNum, pointNum)
	{
		var hash = "" + seriesNum + ":" + pointNum;
		this.flot.highlight(seriesNum, pointNum);
		var point = new Object();
		point.seriesNum = seriesNum;
		point.pointNum = pointNum;	
		this.selectedPoints[hash] = point;
	}
	
	this.unselectAll = function()
	{
		this.flot.unhighlight();
		this.selectedPoints = new Object();
		this.onSelectionChanged();
	}
	
	this.selectRange = function(ranges)
	{
		var p;
		for (var k = self.dataArray.length-1; k>=0; k--)
		{
			if (!(self.dataArray[k].clickable == false))	
			for (var i = self.dataArray[k].data.length-1; i>=0;  i--)
			{
				p = self.dataArray[k].data[i];
				if (p[0] > ranges.xaxis.from && p[0] < ranges.xaxis.to && p[1] > ranges.yaxis.from && p[1] < ranges.yaxis.to)
				{
					var hash = "" + k + ":" + i;
					this.flot.highlight(k, i);
					var point = new Object();
					point.seriesNum = k;
					point.pointNum = i;	
					this.selectedPoints[hash] = point;	
				}
			}
		}
		this.onSelectionChanged();
	}
	
	this.onSelectionChanged = function()
	{
			
	}
	
	this.highlightSelected = function()
	{
		for (var point in this.selectedPoints)
			this.flot.highlight(this.selectedPoints[point].seriesNum, this.selectedPoints[point].pointNum);
	}
	
	this.pointClicked = function(arrayIndex, pointIndex)
	{
		//window.alert(arrayIndex+":"+pointIndex);
	}
	
	this.getSeriesByNum = function(num)
	{
		for (var i = 0; i < self.dataArray.length; i++)
			if (this.dataArray[i].num == 1.0*num)
				return this.dataArray[i];
	}
	
	this.render = function(selector)
	{
		if (selector)
			this.selector = selector;	
		selector = this.selector;
		
		var min = this.minX < this.minY ? this.minX : this.minY;
		var max = this.maxX > this.maxY ? this.maxX : this.maxY;
		var w = max - min;
		//min -= w / 10;
		//max += w / 10;
		
		var options = this.options;
		if (this.normalizeAspectRatio)
			options = $.extend(true, {}, this.options, {
                xaxis: { min: min, max: max },
                yaxis: { min: min, max: max }
            });
		else
			options = $.extend(true, {}, this.options, {
                xaxis: { min: this.minX, max: this.maxX },
                yaxis: { min: this.minY, max: this.maxY }
            });
		this.flot = $.plot($(selector), this.dataArray, options);
		$(selector).bind("plotclick", this.onClick);
		$(selector).get(0).plot = this;
		$(selector).bind("plotselected", function (event, ranges) {
	        // clamp the zooming to prevent eternal zoom
	        if (ranges.xaxis.to - ranges.xaxis.from < 0.00001)
	            ranges.xaxis.to = ranges.xaxis.from + 0.00001;
	        if (ranges.yaxis.to - ranges.yaxis.from < 0.00001)
	            ranges.yaxis.to = ranges.yaxis.from + 0.00001;
	        
			if (self.shiftKey)
				self.selectRange(ranges);				
			else
			{
		        // do the zooming
		        self.zoomTo(ranges.xaxis.from, ranges.xaxis.to, ranges.yaxis.from, ranges.yaxis.to);
			}
	      
   		 });
		
		this.highlightSelected();
		
		
		// Tooltips
		if (this.tooltips)
		{
			$(selector).bind("plothover", function (event, pos, item) 
			{
		        $("#x").text(pos.x.toFixed(2));
		        $("#y").text(pos.y.toFixed(2));
		
	            if (item) 
	            {
	                if (self.previousPoint != item.datapoint) 
	                {
	                    self.previousPoint = item.datapoint;
	                    
	                    $("#flot-tooltip").remove();
	                    var x = item.datapoint[0].toFixed(2),
	                        y = item.datapoint[1].toFixed(2);
	                    
	                    self.showTooltip(item.pageX, item.pageY,
	                                item.series.label);
	                }
	            }
	            else 
	            {
	                $("#flot-tooltip").remove();
	                self.previousPoint = null;            
	            }
    		});
		}
		
		// Legend clicks
		$(selector).find("span[series]").click(function(){
			for (var i = 0; i < self.dataArray.length; i++)
				if (self.dataArray[i].num == 1.0*$(this).attr("series"))
					self.onLegendClick(self.dataArray[i]);
		});
	}
	
	this.enlarge = function()
	{
		//$(this.selector).css({width: "700", height: "500"});
		//this.flot = $.plot($(this.selector), this.dataArray, this.options);
	}
	
	this.clearAutoHighlights = function(autoTag)
	{
		if (self.flot.highlights) // Uses my custom flot hack - reveal the private variable "highlights"
			for (var i = 0; i < self.flot.highlights.length; ++i) 
			{
	            var h = self.flot.highlights[i];
	            if (h.auto == autoTag)
	            	self.flot.unhighlight(h.series, h.point);
	        }
	}
	
	this.zoomTo = function(minX, maxX, minY, maxY, border)
	{
		if (border)
		{
			w = maxX - minX;
			h = maxY - minY;
			minX -= w*0.02;
			maxX += w*0.02;
			minY -= h*0.02;
			maxY += h*0.02;
		}
		this.flot = $.plot($(self.selector),
	        			self.dataArray,
	        			//self.getData(ranges.xaxis.from, ranges.yaxis.from, ranges.xaxis.to, ranges.yaxis.to),
	                      $.extend(true, {}, self.options, {
	                          xaxis: { min: minX, max: maxX },
	                          yaxis: { min: minY, max: maxY }
	                      }));
		this.highlightSelected();
	        
	        self.minX = minX;
	        self.maxX = maxX;
	        self.minY = minY;
	        self.maxY = maxY;	
	}
	
	this.showTooltip = function(x, y, contents) {
        $('<div id="flot-tooltip">' + contents + '</div>').css( {
            position: 'absolute',
            display: 'none',
            top: y + 5,
            left: x + 5,
            border: '1px solid #fdd',
            padding: '2px',
            'background-color': '#fee',
            opacity: 0.80
        }).appendTo("body").fadeIn(200);
    }

}

Plot.getPlot = function(selector)
{
	return $(selector).get(0).plot;
}

function PlotSet()
{
	this.plots = new Array();
	
	this.minX = 100000;
	this.maxX = -10000;
	this.minY = 10000;
	this.maxY = -10000;
	
	this.alignRanges = function()
	{
		for (var i = 0; i < this.plots.length; i++)
		{
			if (this.plots[i].dataRange.minX < this.minX)
				this.minX = this.plots[i].dataRange.minX;
			if (this.plots[i].dataRange.maxX > this.maxX)
				this.maxX = this.plots[i].dataRange.maxX;
			if (this.plots[i].dataRange.minY < this.minY)
				this.minY = this.plots[i].dataRange.minY;
			if (this.plots[i].dataRange.maxY > this.maxY)
				this.maxY = this.plots[i].dataRange.maxY;
		}
		
		for (var i = 0; i < this.plots.length; i++)
		{
			this.plots[i].dataRange.minX = this.minX;
			this.plots[i].dataRange.maxX = this.maxX;
			this.plots[i].dataRange.minY = this.minY;
			this.plots[i].dataRange.maxY = this.maxY;
			this.plots[i].zoomTo(this.minX, this.maxX, this.minY, this.maxY);
			//this.plots[i].zoomTo(-1, 1, -1, 1);
		}
		window.alert(this.minY);
	}
}

Object.size = function(obj) {
    var size = 0, key;
    for (key in obj) {
        if (obj.hasOwnProperty(key)) size++;
    }
    return size;
};

var setColors = new Object();
setColors["training"] = "#C00";
setColors["validation"] = "#0C0";
setColors["excluded"] = "#00C";
setColors["modified-on-the-fly"] = "#CC0";
setColors["deleted"] = "#666";

