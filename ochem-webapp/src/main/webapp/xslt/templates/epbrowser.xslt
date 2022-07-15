<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	version="1.0">
	<xsl:include href="../helper.xslt" />
	<xsl:include href="inc/substructure-search.xslt" />

	<xsl:template name="content">
		<title>Compounds properties browser</title>
		<style type="text/css">
			TABLE.cfilter {margin-top: 5px;}
			img {vertical-align:
			middle;}
			.padded {margin-left: 20px !important;}
			.yui-skin-sam
			.yui-dialog .ft span.default button {
			color:#000000; !important;
			font-weight:bold; !important;
			}

			#EPFilters {position: relative;}
			#EPFilters H2 {font-weight: bold; margin-top:
			10px;}

			#show-filters
			{
			position: fixed;
			left: 0px;
			top: 50%;
			background-color: #44F;
			color: white;
			}

			.selection
			{
			background-color: #e3eaf6;
			border: 1px solid #222;
			padding: 10px;
			clear: both;
			}

			.shift {margin-left: 15px;}

			.content-message {
			text-align: center;
			font-size: 130%;
			color: #888;
			padding: 50px 0px;
			border: 1px solid gray;
			}

			#sorting-discarded {
			border: 1px solid gray;
			background-color: #EEE;
			color: #900;
			padding: 5px;
			}

			SMALL.gray {
			color: #666;
			}


		</style>

		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/browsers/ep-browser.js?ver=2.3.5"></script>
		<script language="javascript" src="js/commons/long-operation.js"></script>

		<script language="javascript">
			include.plugins('view');
			var compoundBrowser =
			new CompoundBrowser();
			var conditionsFilter = new ConditionsFilter();
			//var basketActionable = new BasketActionable();
			$(document).ready(function() {
			if (getParams["trash"])
			$("#deleteselected").attr("title", "Delete selected records permanently");
			compoundBrowser.options.highlight_row = "no";
			compoundBrowser.initialize();
			compoundBrowser.recordMenu = new YAHOO.widget.Menu("basicmenu");
			compoundBrowser.recordMenu.render();
			compoundBrowser.publishMenu = new YAHOO.widget.Menu("publishmenu");
			compoundBrowser.publishMenu.render();
			conditionsFilter.initialize(true);
			sketcherDialog.filters = compoundBrowser.filters;

			if (getParams["xemistry-smiles"])
			sketcherDialog.setMolecule(getParams["xemistry-smiles"]);

			//$('a[title]').tooltip({showURL: false, delay: 0});
			//basketActionable.initialize();
			});
		</script>
		<table width="100%" height="100%" cellspacing="0">
			<tr>
				<td class="itunes-up" colspan="2">

					<img src="img/icons/compounds.png" id="page-icon" />
					<xsl:call-template name="area-of-interest" />
					<h1>
						<doc term="Compound+property+browser">Compounds properties browser</doc>
					</h1>
					Search for numerical compounds properties linked to scientific
					articles
				</td>
			</tr>
			<tr>
				<td class="browserScope itunes-left" id="EPFilters">
					<div id="EPFiltersDiv">
						<a action="hidefilters" class="close-button" title="Hide filters panel">x</a>
						<div class="selection invisible">
							<nobr>
								<img src="img/icons/cart.png" />
								Your saved selection contains
								<a
									href="epbrowser/show.do?basket-select=selected&amp;render-mode=popup">
									<span id="selection-size"></span>
									records
								</a>
								<a class="delete-link" action="clearselection">[clear]</a>
							</nobr>
						</div>
						<h1>Filters</h1>
						<div class="openable opened">
							<h1>Source</h1>
							<div class="openable-content" style="min-height: 74px;">
								Article/Source
								<a href="javascript:void(0)" action="editarticle" name="article-link">[select]</a>
								<div id="ac-source" style="height: 25px;"></div>
								<table>
									<tr>
										<td>N: (id)</td>
										<td>Table</td>
									</tr>
									<tr>
										<td>
											<input class="small right" name="idarticle" type="text"
												prompt="All" value="" filter="1" />
										</td>
										<td>
											<input class="small right" name="table" type="text"
												prompt="All" value="" filter="1" />
										</td>
									</tr>
								</table>
							</div>
						</div>
						<div class="openable opened">
							<h1>Property</h1>
							<div class="openable-content">
								Activity/Property
								<a action="editproperty" name="property-link">[select]</a>
								<div id="ac-property" style="height: 25px;"></div>
								<input type="checkbox" filter="1" name="hidedummy" />
								Hide records without property
							</div>
						</div>
						<div class="openable">
							<h1>Conditions</h1>
							<div class="openable-content conditionsscope">
								<a action="newfilter" restrict="conditionsfilter">[add]</a>
								<div id="ConditionsFilter"></div>
							</div>
						</div>

						<div class="openable opened">
							<h1>Molecule filters</h1>
							<div class="openable-content">
								<nobr>
									Name / OCHEM ID
									<a help="OCHEM-id-hint" class="infolink">&#160;</a>
									/ Inchi-Key &#160;
									<input type="text" name="name" prompt="" value="" filter="1" />
								</nobr>
								<span class="invisible" id="OCHEM-id-hint">
									Examples of OCHEM IDs:
									<br />
									<br />
									M1234 for a specific molecule
									<br />
									M1234* for a molecule and all its isomers
									<br />
									R1234 for a particular record
									<br />
									R1234* for a record and all the similar records
									<br />
								</span>

								<div style="clear: both;">
									<xsl:call-template name="substructure-search" />
								</div>
								<br class="cb" />
								<br />
								<br />

								Molecular mass
								<a class="infolink"
									title="as unified atomic mass unit (amu / u) or dalton (Da)"></a>
								<br />
								<nobr>
									between
									<input type="text" class="w50" name="mmlo" prompt=""
										value="" filter="1" />
									and
									<input type="text" class="w50" name="mmhi" prompt=""
										value="" filter="1" />
								</nobr>
							</div>

						</div>

						<div class="openable">
							<h1>Advanced molecule filters</h1>
							<div class="openable-content">
								<input type="checkbox" id="structure-basket-cb" />
								<label for="structure-basket-cb">Only the molecules present in a basket</label>
								<a class="infolink"
									title="Only the molecules present in the basket will be filtered. Check the box to select the basket with molecules"></a>
								<div id="structure-basket-div" class="invisible">
									<input type="hidden" name="structure-basket-dofilter"
										filter="1" value="0" />
									<input type="hidden" name="structure-basket" filter="1"
										value="-1" />
									Basket with molecules:
									<a action="selectstructurebasket" name="selectstructurebasket-link">[...]</a>
								</div>
							</div>
						</div>

						<div class="openable opened">
							<h1>Miscellaneous</h1>
							<div id="modification-time" class="invisible">
								<div title="User changes">Modification record id, minutes</div>
								<input type="text" name="interval" prompt="" value=""
									filter="1" />
								<br />
							</div>
							<div class="openable-content">
								<doc term="Working+with+datasets"
									title="Filter the records from a particular set ('basket')">Current set</doc>

								<select id="basket-menu" name="basket-select" filter="1">
									<option value="-1">Loading...</option>
								</select>
								<input type="hidden" name="compare-baskets" value=""
									filter="1" />
								<input type="hidden" name="join-baskets" value="" filter="1" />
								<input type="hidden" name="unit" filter="1" />
								<br />

								<br />
								<b>Data origin and quality:</b>
								<br />
								<div class="shift">
									<xsl:if test="/model/session/user">
										Data introducers:
										<select name="introducer" filter="1">
											<option value="">All users</option>
											<optgroup label="My records">
												<option value="my-records">Records I introduced or modified</option>
												<option value="{/model/session/user/@id}">Records I introduced</option>
												<option value="modified-by-me">Records I modified</option>
												<option value="modified-by-others">My records modified by others</option>
											</optgroup>
											<xsl:if test="/model/session/user/group">
												<optgroup label="{/model/session/user/group/name}">
													<option value="group">
														<xsl:value-of select="/model/session/user/group/name" />
														records
													</option>
													<xsl:for-each select="/model/session/user/group/member">
														<option value="{@id}">
															<xsl:value-of select="@login" />
														</option>
													</xsl:for-each>
												</optgroup>
											</xsl:if>
										</select>
										<br />
									</xsl:if>

									Data visibility:
									<select name="visibility" filter="1" style="margin-left: 40pt;">
										<option value="all">All data</option>
										<option value="public">Freely available</option>
										<option value="private">Private only</option>
									</select>
									<br />

									Data from other users:&#160;&#160;
									<select name="approval-status" filter="1">
										<option value="all">All data</option>
										<option value="only-approved">Only approved data</option>
										<option value="only-awaiting-approval">Only awaiting approval</option>
										<option value="rejected">Rejected by moderator</option>
									</select>
									<br />
									<input type="checkbox" name="experimental" value="true"
										filter="1" class="checkbox" />
									Original records
									<br />
									<input type="checkbox" name="primaryrecords" value="true"
										filter="1" class="checkbox" />
									Primary records
									<br />
								</div>
								<!-- <input type="checkbox" name="validated" value="true" filter="1" 
									class="checkbox"/> Not validated <br/> -->
								<br />
								<b>Discover issues with the data:</b>
								<br />
								<div class="shift">
									<input type="checkbox" name="error" value="true" filter="1"
										class="checkbox" />
									Error records
									<br />
									<input type="checkbox" name="inchierror" value="true"
										filter="1" class="checkbox" />
									Error inchies
									<br />
									<input type="checkbox" name="nameerror" value="true"
										filter="1" class="checkbox" />
									Mismatching names
									<br />
									<input type="checkbox" name="yesstereo" value="true"
										filter="1" class="checkbox padded" disabled="true" />
									Include stereochem.
									<br />
									<input type="checkbox" name="emptymol" value="true"
										filter="1" class="checkbox" />
									Empty molecules

									<div id="duplicate-filter" class="invisible">
										<br />
										<input type="checkbox" name="duplicates" value="true"
											filter="1" class="checkbox" />
										Show only duplicates
										<br />
										<input type="checkbox" name="nostereo" value="true"
											filter="1" class="checkbox padded" disabled="true" />
										No stereochemistry
									</div>
								</div>
								<br />
								<br />
								Sort by:
								<br />
								<nobr>
									<select name="sortby" filter="1">
										<option value="">Creation time</option>
										<option value="time">Modification time</option>
										<option value="value">Numerical value</option>
										<option value="canonicalValue">Converted numerical value</option>
										<option value="molprop">Compound</option>
										<option value="mol.molWeight">Molecular weight</option>
										<option value="property">Property</option>
										<option value="artMolId">Identifier in article</option>
										<option value="artLineNum">Line in article</option>
									</select>
									<input type="checkbox" name="order" value="asc" filter="1" />
									&#160;Ascending order
								</nobr>
								<div id="sorting-discarded" class="invisible"></div>
							</div>
						</div>
						<input type="hidden" name="predicate" value="" filter="1"
							class="invisible" />
						<input type="hidden" name="id" value="" filter="1" class="invisible"
							id="fixed-id" />
						<input type="hidden" name="molid" filter="1"
							value="{filter[@name='molid']/@value}" class="invisible" />
						<input type="hidden" name="value" filter="1"
							value="{filter[@name='value']/@value}" class="invisible" />
						<input type="hidden" name="canonicalValue" filter="1"
							value="{filter[@name='canonicalValue']/@value}" class="invisible" />
						<br />
						<a href="javascript:void()" class="button-link"
							onclick="compoundBrowser.request(true); return false;">refresh</a>
						<a action="reset" class="button-link">reset</a>
						<br />

					</div>
				</td>
				<td class="itunes-right">
					<div class="browserScope" style="position: relative;">
						<a action="showfilters" class="show-button invisible" title="Show filters panel">&lt;</a>
						<div class="upper-command-panel">
							<label class="no-trash">Basket</label>
							<a action="basketaction" type="addbasket"
								title="Add selected records matching current filters to basket"
								class="no-trash">
								<img src="img/icons/add_to_set.gif" />
							</a>
							<a action="basketaction" type="removebasket"
								title="Remove selected records matching current filters from a basket"
								class="no-trash">
								<img src="img/icons/remove_from_set.gif" />
							</a>
							<a action="cleanbasketaction" type="cleanbasket" id="clean-basket"
								class="invisible" title="select duplicates within the basket">
								<img src="img/icons/duplicate.png" />
							</a>
							<label>Records</label>
							<a action="addselect" title="Select all records matching current filters"
								class="invisible">
								<img src="img/icons/select_all.gif" />
							</a>
							<a action="selectpage" title="Select records on currently visible page">
								<img src="img/icons/select_page.gif" />
							</a>
							<a action="removeselect" title="Unselect all records matching current filters">
								<img src="img/icons/unselect_all.gif" />
							</a>
							<a action="edit" title="Click to create new empty record"
								class="no-trash">
								<img src="img/icons/new.gif" />
							</a>
							<a action="publishmenu" class="no-trash" id="publish-menu-link">
								<img src="img/icons/public32.png" />
							</a>
							<a href="batchupload30/show.do?render-mode=full" target="_parent"
								title="Upload a set of properties/molecules from an excel or sdf file"
								class="no-trash">
								<img src="img/icons/batch_upload.gif" />
							</a>
							<!-- <a tab="Batch editing" href="batchedit/show.do?apply-basket=true" 
								title="Batch edit selected molecules" type="batchedit" class="no-trash"><img 
								src="img/icons/edit.gif"/></a> -->
							<a action="batchedit" title="Batch edit selected molecules"
								type="batchedit" class="no-trash">
								<img src="img/icons/edit.gif" />
							</a>
							<a action="deleteselected"
								title="Move selected records matching current filters to the trash"
								id="deleteselected">
								<img src="img/icons/delete.gif" />
							</a>
							<span class="trash invisible">
								<a action="restoreselected">
									<img src="img/icons/restore24.png" />
								</a>
							</span>
							<xsl:if test="session/user">
								<a href="epbrowser/show.do?trash=1&amp;render-mode=full"
									title="View your trash" target="_parent" class="no-trash" id="trash-button">
									<img src="img/icons/trash24.png" />
								</a>
							</xsl:if>
							<label>Tags</label>
							<a action="selecttag"
								title="Add a tag to selected compounds matching current filters">
								<img src="img/icons/tag.png" />
							</a>
							<a action="selecttag"
								title="Remove a tag from selected compounds matching current filters"
								removetag="1">
								<img src="img/icons/tag-orange.png" />
							</a>

							<div style="float: right" class="no-trash">
								<xsl:if
									test="/model/session/@moderator = 'true' or /model/session/user/superuser = 'true'">
									<a action="approve"
										title="Approve the selected records matching current filters">
										<img src="img/icons/approve32.png" />
									</a>
									<a action="disapprove"
										title="Reject the selected records matching current filters. The records will remain public and will have to be reviewed by the user who submitted them.">
										<img src="img/icons/disapprove32.png" />
									</a>
									<a action="unapprove"
										title="Mark the previously approved or rejected records as 'awaiting approval'.">
										<img src="img/icons/unapprove32.png" />
									</a>
									<!-- <input type="checkbox" filter="1" name="awaiting-approval"/>&#160;Awaiting 
										approval&#160;&#160; <span id="unmoderated"><input type="checkbox" filter="1" 
										name="unmoderated"/>&#160;Unmoderated data</span> -->
								</xsl:if>
							</div>
						</div>
						<div class="pager-strip">
							<span>
								<b class="showed">none</b>
								of
								<b class="total">none</b>
							</span>
							<div id="pager" class="pgr">
							</div>
						</div>
						<div id="Browser">
						</div>
						<div>
							<div id="query-status">
							</div>
						</div>
						<div class="pager-strip">
							<span>
								<b class="showed">none</b>
								of
								<b class="total">none</b>
							</span>
							<div id="pager" class="pgr">
							</div>
						</div>

					</div>
				</td>
			</tr>
		</table>

		<div id="basketDialog">
			<div class="hd">Select your molecule set</div>
			<div class="bd">
				<form>
					<label for="n-basket">Molecule set:</label>
					<select name="n-basket" onChange="compoundBrowser.onSelectchange()">
						<option>Loading items...</option>
					</select>
					<div id="basketName" style="display: none;">
						<label for="basket-name">New set name:</label>
						<input type="text" name="basket-name" id="basketNameEdit" />
						<br />
					</div>
					<div id="basketOption">
						<input type="radio" name="mol-option" value="3" checked="true" />
						all selected record
						<br />
						<input type="radio" name="mol-option" value="1" />
						molecule (preserve stereochemistry)
						<br />
						<input type="radio" name="mol-option" value="2" />
						molecule (ignore stereochemistry)
					</div>
				</form>
			</div>
		</div>

		<div id="molViewDialog">
			<div class="hd">Molecule structure</div>
			<div class="bd">
				<a href="javascript:void()" width="300" height="300"
					onclick="compoundBrowser.molViewDialogClose(); return false;">
					<img src="" />
				</a>
			</div>
		</div>

		<div id="cleanBasketDialog">
			<div class="hd">Select record</div>
			<div class="bd">
				<form>
					<div id="recordOption">
						<input type="radio" name="withStreo" value="true" checked="true" />
						molecule (preserve stereochemistry)
						<br />
						<input type="radio" name="withStreo" value="false" />
						molecule (ignore stereochemistry)
						<br />
						<input type="checkbox" name="withValue" value="true" />
						with value
						<br />
					</div>
				</form>
			</div>
		</div>

		<div id="waitingDialog">
			<div class="hd">Please wait</div>
			<div class="bd" style="text-align: center;">
				Please wait until action is completed.
				<br />
				It may take a while.
				<br />
				<img src="img/roller_small.gif" />
			</div>
		</div>

		<div id="restoreDialog">
			<div class="hd">Message</div>
			<div class="bd">
				<form>
					<input type="hidden" name="rec-id" value=""></input>
					<div id="dupRecords">duplicate records</div>
				</form>
			</div>
		</div>

		<div id="basicmenu" class="yuimenu browserScope">
			<div class="bd">
				<ul class="first-of-type">
					<li class="yuimenuitem no-trash">
						<a action="edit" class="yuimenuitemlabel" href="#trackchanges">
							<img src="img/icons/edit.gif" />
							Edit record
						</a>
					</li>
					<li class="yuimenuitem no-trash">
						<a action="clone" class="yuimenuitemlabel" href="#clone">
							Clone record
						</a>
					</li>
					<li class="yuimenuitem no-trash">
						<a action="compounddetails" class="yuimenuitemlabel" href="#clone">
							Open the molecule profile
						</a>
					</li>
					<li class="yuimenuitem">
						<a class="yuimenuitemlabel" href="#changes">
							Changes
						</a>
						<div id="changes" class="yuimenu">
							<div class="bd">
								<ul>
									<li class="yuimenuitem">
										<a action="recordchanges" interval="10" class="yuimenuitemlabel">
											&#x00B1; 10 mins
										</a>
									</li>
									<li class="yuimenuitem">
										<a action="recordchanges" interval="20" class="yuimenuitemlabel">
											&#x00B1; 20 mins
										</a>
									</li>
									<li class="yuimenuitem">
										<a action="recordchanges" interval="30" class="yuimenuitemlabel">
											&#x00B1; 30 mins
										</a>
									</li>
									<li class="yuimenuitem">
										<a action="recordchanges" interval="60" class="yuimenuitemlabel">
											&#x00B1; 60 mins
										</a>
									</li>
								</ul>
							</div>
						</div>
					</li>
				</ul>
				<ul>
					<li class="yuimenuitem no-trash">
						<a action="showsimilarcompounds" class="yuimenuitemlabel" href="#showsimilar">
							Find same compounds (ignore stereochemistry)
						</a>
					</li>
					<li class="yuimenuitem no-trash">
						<a action="showsamecompounds" class="yuimenuitemlabel" href="#showsimilar">
							Find same compounds (preserve stereochemistry)
						</a>
					</li>
					<li class="yuimenuitem no-trash">
						<a action="showdublicates" class="yuimenuitemlabel" href="#showdublicates">
							Find duplicates
						</a>
					</li>
					<li class="yuimenuitem">
						<a action="trackchanges" class="yuimenuitemlabel" href="#trackchanges">
							Track changes
						</a>
					</li>
				</ul>
			</div>
		</div>

		<div id="publishmenu" class="yuimenu browserScope">
			<div class="bd">
				<h6 class="first-of-type">
					Publish the selected records matching current filters.
					<br />
					Please, select the desired option:
				</h6>
				<ul class="first-of-type">
					<!--
					<li class="yuimenuitem no-trash">
						<a action="publish" class="yuimenuitemlabel" href="#publish">
							Publish
						</a>
					</li>
					-->
					<li class="yuimenuitem no-trash">
						<a action="publish_freely_available" class="yuimenuitemlabel"
							href="#publish-freely-avaliable">
							Publish and make freely available for download
						</a>
					</li>

				</ul>
				<ul>
					<li class="yuimenuitem no-trash">
						<a class="yuimenuitemlabel"
							href="https://docs.ochem.eu/display/MAN/Access+levels%3A+private%2C+public+and+free+data"
							target="_blank">
							<img src="img/icons/info-16.png" />
							Read more about public data on OCHEM
						</a>
					</li>
				</ul>
			</div>
		</div>

	</xsl:template>
</xsl:stylesheet>