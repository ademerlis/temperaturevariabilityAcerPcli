// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, true, false, false ];
var arrayMetadata    = [ [ "1", "Pcli-003", "C", "3", "4/20/22", "Untreated" ], [ "2", "Pcli-004", "A", "4", "4/20/22", "Untreated" ], [ "3", "Pcli-006", "B", "6", "4/20/22", "Treated" ], [ "4", "Pcli-008", "A", "8", "4/20/22", "Treated" ], [ "5", "Pcli-011", "A", "11", "3/22/22", "Initial" ], [ "6", "Pcli-015", "C", "15", "4/20/22", "Treated" ], [ "7", "Pcli-018", "C", "18", "4/20/22", "Treated" ], [ "8", "Pcli-024", "B", "24", "4/20/22", "Untreated" ], [ "9", "Pcli-026", "A", "26", "3/22/22", "Initial" ], [ "10", "Pcli-028", "C", "28", "3/22/22", "Initial" ], [ "11", "Pcli-029", "B", "29", "3/22/22", "Initial" ], [ "12", "Pcli-033", "A", "33", "3/22/22", "Initial" ], [ "13", "Pcli-038", "C", "38", "4/20/22", "Untreated" ], [ "14", "Pcli-039", "C", "39", "4/20/22", "Treated" ], [ "15", "Pcli-040", "C", "40", "3/22/22", "Initial" ], [ "16", "Pcli-045", "B", "45", "3/22/22", "Initial" ], [ "17", "Pcli-049", "B", "49", "4/20/22", "Untreated" ], [ "18", "Pcli-050", "A", "50", "4/20/22", "Untreated" ], [ "19", "Pcli-052", "B", "52", "3/22/22", "Initial" ], [ "20", "Pcli-055", "B", "55", "3/22/22", "Initial" ], [ "21", "Pcli-058", "A", "58", "4/20/22", "Untreated" ], [ "22", "Pcli-064", "C", "64", "4/20/22", "Untreated" ], [ "23", "Pcli-066", "B", "66", "4/20/22", "Untreated" ], [ "24", "Pcli-067", "C", "67", "3/22/22", "Initial" ], [ "25", "Pcli-069", "A", "69", "4/20/22", "Untreated" ], [ "26", "Pcli-070", "B", "70", "3/22/22", "Initial" ], [ "27", "Pcli-075", "A", "75", "3/22/22", "Initial" ], [ "28", "Pcli-076", "C", "76", "4/20/22", "Treated" ], [ "29", "Pcli-077", "A", "77", "3/22/22", "Initial" ], [ "30", "Pcli-079", "B", "79", "4/20/22", "Treated" ], [ "31", "Pcli-091", "C", "91", "3/22/22", "Initial" ], [ "32", "Pcli-094", "C", "94", "3/22/22", "Initial" ], [ "33", "Pcli-096", "B", "96", "4/20/22", "Treated" ], [ "34", "Pcli-101", "A", "101", "4/20/22", "Treated" ], [ "35", "Pcli-103", "C", "103", "3/22/22", "Initial" ], [ "36", "Pcli-105", "B", "105", "3/22/22", "Initial" ], [ "37", "Pcli-109", "B", "109", "4/20/22", "Treated" ], [ "38", "Pcli-111", "A", "111", "3/22/22", "Initial" ], [ "39", "Pcli-114", "B", "114", "3/22/22", "Initial" ], [ "40", "Pcli-115", "C", "115", "3/22/22", "Initial" ], [ "41", "Pcli-120", "B", "120", "3/22/22", "Initial" ], [ "42", "Pcli-124", "B", "124", "4/20/22", "Untreated" ], [ "43", "Pcli-126", "A", "126", "3/22/22", "Initial" ], [ "44", "Pcli-128", "A", "128", "3/22/22", "Initial" ], [ "45", "Pcli-130", "A", "130", "4/20/22", "Treated" ], [ "46", "Pcli-132", "A", "132", "4/20/22", "Treated" ], [ "47", "Pcli-135", "C", "135", "4/20/22", "Untreated" ], [ "48", "Pcli-148", "C", "148", "3/22/22", "Initial" ] ];
var svgObjectNames   = [ "pca", "dens" ];

var cssText = ["stroke-width:1; stroke-opacity:0.4",
               "stroke-width:3; stroke-opacity:1" ];

// Global variables - these are set up below by 'reportinit'
var tables;             // array of all the associated ('tooltips') tables on the page
var checkboxes;         // the checkboxes
var ssrules;


function reportinit() 
{
 
    var a, i, status;

    /*--------find checkboxes and set them to start values------*/
    checkboxes = document.getElementsByName("ReportObjectCheckBoxes");
    if(checkboxes.length != highlightInitial.length)
	throw new Error("checkboxes.length=" + checkboxes.length + "  !=  "
                        + " highlightInitial.length="+ highlightInitial.length);
    
    /*--------find associated tables and cache their locations------*/
    tables = new Array(svgObjectNames.length);
    for(i=0; i<tables.length; i++) 
    {
        tables[i] = safeGetElementById("Tab:"+svgObjectNames[i]);
    }

    /*------- style sheet rules ---------*/
    var ss = document.styleSheets[0];
    ssrules = ss.cssRules ? ss.cssRules : ss.rules; 

    /*------- checkboxes[a] is (expected to be) of class HTMLInputElement ---*/
    for(a=0; a<checkboxes.length; a++)
    {
	checkboxes[a].checked = highlightInitial[a];
        status = checkboxes[a].checked; 
        setReportObj(a+1, status, false);
    }

}


function safeGetElementById(id)
{
    res = document.getElementById(id);
    if(res == null)
        throw new Error("Id '"+ id + "' not found.");
    return(res)
}

/*------------------------------------------------------------
   Highlighting of Report Objects 
 ---------------------------------------------------------------*/
function setReportObj(reportObjId, status, doTable)
{
    var i, j, plotObjIds, selector;

    if(doTable) {
	for(i=0; i<svgObjectNames.length; i++) {
	    showTipTable(i, reportObjId);
	} 
    }

    /* This works in Chrome 10, ssrules will be null; we use getElementsByClassName and loop over them */
    if(ssrules == null) {
	elements = document.getElementsByClassName("aqm" + reportObjId); 
	for(i=0; i<elements.length; i++) {
	    elements[i].style.cssText = cssText[0+status];
	}
    } else {
    /* This works in Firefox 4 */
    for(i=0; i<ssrules.length; i++) {
        if (ssrules[i].selectorText == (".aqm" + reportObjId)) {
		ssrules[i].style.cssText = cssText[0+status];
		break;
	    }
	}
    }

}

/*------------------------------------------------------------
   Display of the Metadata Table
  ------------------------------------------------------------*/
function showTipTable(tableIndex, reportObjId)
{
    var rows = tables[tableIndex].rows;
    var a = reportObjId - 1;

    if(rows.length != arrayMetadata[a].length)
	throw new Error("rows.length=" + rows.length+"  !=  arrayMetadata[array].length=" + arrayMetadata[a].length);

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = arrayMetadata[a][i];
}

function hideTipTable(tableIndex)
{
    var rows = tables[tableIndex].rows;

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = "";
}


/*------------------------------------------------------------
  From module 'name' (e.g. 'density'), find numeric index in the 
  'svgObjectNames' array.
  ------------------------------------------------------------*/
function getIndexFromName(name) 
{
    var i;
    for(i=0; i<svgObjectNames.length; i++)
        if(svgObjectNames[i] == name)
	    return i;

    throw new Error("Did not find '" + name + "'.");
}


/*------------------------------------------------------------
  SVG plot object callbacks
  ------------------------------------------------------------*/
function plotObjRespond(what, reportObjId, name)
{

    var a, i, status;

    switch(what) {
    case "show":
	i = getIndexFromName(name);
	showTipTable(i, reportObjId);
	break;
    case "hide":
	i = getIndexFromName(name);
	hideTipTable(i);
	break;
    case "click":
        a = reportObjId - 1;
	status = !checkboxes[a].checked;
	checkboxes[a].checked = status;
	setReportObj(reportObjId, status, true);
	break;
    default:
	throw new Error("Invalid 'what': "+what)
    }
}

/*------------------------------------------------------------
  checkboxes 'onchange' event
------------------------------------------------------------*/
function checkboxEvent(reportObjId)
{
    var a = reportObjId - 1;
    var status = checkboxes[a].checked;
    setReportObj(reportObjId, status, true);
}


/*------------------------------------------------------------
  toggle visibility
------------------------------------------------------------*/
function toggle(id){
  var head = safeGetElementById(id + "-h");
  var body = safeGetElementById(id + "-b");
  var hdtxt = head.innerHTML;
  var dsp;
  switch(body.style.display){
    case 'none':
      dsp = 'block';
      hdtxt = '-' + hdtxt.substr(1);
      break;
    case 'block':
      dsp = 'none';
      hdtxt = '+' + hdtxt.substr(1);
      break;
  }  
  body.style.display = dsp;
  head.innerHTML = hdtxt;
}
