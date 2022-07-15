function skim(images)
{
	var step = 0;
	var size = images.length;
	if(images.length > 50)
		size = 50;
		
	var stepBy = Math.floor(200 / size);
	var imgkey = "data_key";
	
	var _block1 = '';
	var _block3 = '';
	var _block2 = "<br><img class='block-image' src='depiction.jsp?id="+ images[0] +"' width='200' height='200' usemap='#example"+ imgkey +"' name='example"+ imgkey +"' alt=''><map name='example"+ imgkey +"'>"
	
	for (var i=0; i < size; i++)
	{
		_block1 += '<img src="depiction.jsp?id='+images[i] +'" width="1" height="1" alt=""/ >';		 	 
		
		var head = '<area shape="rect" COORDS="';
		var coord = step + ",0," + eval(step + eval(stepBy - 1)) + ',200"';
		var href = ' href="" onClick="return false" onmouseover="example'+imgkey+'.src=';
		var mouse = "'depiction.jsp?id="+ images[i] +"'";
		var tail = ';" alt="" >';
		_block3 += head + coord + href + mouse + tail;
		step += stepBy;
	}
	var finalString = _block1+_block2+_block3+"</map></br><small>Molecules:- <b>"+images.length+"</b></small>";
	$("#skim").append(finalString);
}
