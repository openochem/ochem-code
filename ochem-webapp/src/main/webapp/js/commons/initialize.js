	// Midnighter

	// Create main popup menu
	
	YAHOO.util.Event.onDOMReady(
	function()
	{
		cmMain = new YAHOO.widget.ContextMenu("maincontextmenu", { trigger: "main-menu" } );
    	cmMain.addItems(
    	[
    		[
    			{text:"Articles list", url:"article.list"},
    			{text: "Introduce new article", url:"article.new"}
    		],
    		[
    			{text:"Authors list", url:"author.list"},
    			{text: "Introduce new author", url:"author.new"}
    		],
    		[
    			{text:"Home",url:"home/show.do"},
    		]
    		]);
    
    	cmMain.render(document.body);
    }
    );
    
    // "Cloth" tablecloth tables
    YAHOO.util.Event.onDOMReady(
    	function()
		{
			$(".tc table tr:even").addClass("even");
			$(".tc table tr:even").addClass("odd");
			$(".tc table tr td").mouseover
			(
				function()
				{
					$(this).parent().find('td').addClass("over");
				}
			);
			$(".tc table tr td").mouseout
			(
				function()
				{
					$(this).parents().find('td').removeClass("over");
				}
			);
		}
	);
	
	// Handle promptable editboxed
    	promptInputsHandlers = function(parent)
		{
			parent.find('input[prompt]').each(function(){promptInputHandler($(this))});
		}
		
		promptInputHandler = function(input)
				{
					if (input.attr('value') == undefined)
					{
						input.addClass("grayed");
						input.attr('value', input.attr('prompt'));					
						//log(input+" "+input.attr('prompt'));
					}
					
					input.focus
					(
						function() 
						{
							if ($(this).attr('value') == $(this).attr('prompt'))
							{
								$(this).removeClass("grayed");
								$(this).attr('value','');
							}
						}
					);
					input.blur
					(
						function() 
						{
							if ($(this).attr('value') == undefined)
							{
								$(this).addClass("grayed");
								$(this).attr('value', $(this).attr('prompt'));
							}
						}
					);
					

				}
		
		//YAHOO.util.Event.onDOMReady(function(){promptInputsHandlers($(document));});
		
		function log(st)
		{
			$('#log').append(st+"<br/>");
		}
		
		// Opened-closed blocks
		YAHOO.util.Event.onDOMReady(
    	function()
		{
			$('.openable h1').click(function(){
				$(this).parent('.openable').toggleClass('opened');
			});
		}
	);