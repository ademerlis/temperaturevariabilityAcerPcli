// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, true, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false ];
var arrayMetadata    = [ [ "1", "Acer-005", "SI-C", "5", "3/22/22", "Pre-treatment", "control", "Day_0", "control_Day_0" ], [ "2", "Acer-009", "BC-8b", "9", "3/22/22", "Pre-treatment", "variable", "Day_0", "variable_Day_0" ], [ "3", "Acer-019", "BC-8b", "19", "3/22/22", "Pre-treatment", "variable", "Day_0", "variable_Day_0" ], [ "4", "Acer-021", "MB-B", "21", "3/22/22", "Pre-treatment", "variable", "Day_0", "variable_Day_0" ], [ "5", "Acer-022", "BC-8b", "22", "3/22/22", "Pre-treatment", "control", "Day_0", "control_Day_0" ], [ "6", "Acer-029", "SI-C", "29", "3/22/22", "Pre-treatment", "variable", "Day_0", "variable_Day_0" ], [ "7", "Acer-030", "BC-8b", "30", "4/20/22", "last day of treatment", "variable", "Day_29", "variable_Day_29" ], [ "8", "Acer-034", "SI-C", "34", "3/22/22", "Pre-treatment", "control", "Day_0", "control_Day_0" ], [ "9", "Acer-038", "MB-B", "38", "4/20/22", "last day of treatment", "variable", "Day_29", "variable_Day_29" ], [ "10", "Acer-041", "MB-B", "41", "4/20/22", "last day of treatment", "variable", "Day_29", "variable_Day_29" ], [ "11", "Acer-042", "MB-B", "42", "3/22/22", "Pre-treatment", "control", "Day_0", "control_Day_0" ], [ "12", "Acer-043", "BC-8b", "43", "3/22/22", "Pre-treatment", "control", "Day_0", "control_Day_0" ], [ "13", "Acer-045", "BC-8b", "45", "3/22/22", "Pre-treatment", "variable", "Day_0", "variable_Day_0" ], [ "14", "Acer-055", "MB-B", "55", "3/22/22", "Pre-treatment", "variable", "Day_0", "variable_Day_0" ], [ "15", "Acer-056", "MB-B", "56", "4/20/22", "last day of treatment", "control", "Day_29", "control_Day_29" ], [ "16", "Acer-057", "MB-B", "57", "4/20/22", "last day of treatment", "control", "Day_29", "control_Day_29" ], [ "17", "Acer-059", "MB-B", "59", "4/20/22", "last day of treatment", "variable", "Day_29", "variable_Day_29" ], [ "18", "Acer-063", "BC-8b", "63", "3/22/22", "Pre-treatment", "control", "Day_0", "control_Day_0" ], [ "19", "Acer-064", "SI-C", "64", "3/22/22", "Pre-treatment", "variable", "Day_0", "variable_Day_0" ], [ "20", "Acer-068", "BC-8b", "68", "4/20/22", "last day of treatment", "variable", "Day_29", "variable_Day_29" ], [ "21", "Acer-069", "SI-C", "69", "4/20/22", "last day of treatment", "variable", "Day_29", "variable_Day_29" ], [ "22", "Acer-072", "BC-8b", "72", "4/20/22", "last day of treatment", "control", "Day_29", "control_Day_29" ], [ "23", "Acer-074", "MB-B", "74", "4/20/22", "last day of treatment", "variable", "Day_29", "variable_Day_29" ], [ "24", "Acer-078", "SI-C", "78", "3/22/22", "Pre-treatment", "control", "Day_0", "control_Day_0" ], [ "25", "Acer-079", "SI-C", "79", "3/22/22", "Pre-treatment", "variable", "Day_0", "variable_Day_0" ], [ "26", "Acer-082", "SI-C", "82", "4/20/22", "last day of treatment", "control", "Day_29", "control_Day_29" ], [ "27", "Acer-083", "MB-B", "83", "3/22/22", "Pre-treatment", "variable", "Day_0", "variable_Day_0" ], [ "28", "Acer-089", "BC-8b", "89", "4/20/22", "last day of treatment", "control", "Day_29", "control_Day_29" ], [ "29", "Acer-091", "BC-8b", "91", "4/20/22", "last day of treatment", "variable", "Day_29", "variable_Day_29" ], [ "30", "Acer-094", "MB-B", "94", "3/22/22", "Pre-treatment", "variable", "Day_0", "variable_Day_0" ], [ "31", "Acer-095", "SI-C", "95", "3/22/22", "Pre-treatment", "control", "Day_0", "control_Day_0" ], [ "32", "Acer-096", "SI-C", "96", "4/20/22", "last day of treatment", "variable", "Day_29", "variable_Day_29" ], [ "33", "Acer-097", "MB-B", "97", "4/20/22", "last day of treatment", "control", "Day_29", "control_Day_29" ], [ "34", "Acer-102", "MB-B", "102", "3/22/22", "Pre-treatment", "control", "Day_0", "control_Day_0" ], [ "35", "Acer-106", "SI-C", "106", "4/20/22", "last day of treatment", "control", "Day_29", "control_Day_29" ], [ "36", "Acer-108", "BC-8b", "108", "4/20/22", "last day of treatment", "control", "Day_29", "control_Day_29" ], [ "37", "Acer-110", "BC-8b", "110", "3/22/22", "Pre-treatment", "variable", "Day_0", "variable_Day_0" ], [ "38", "Acer-112", "MB-B", "112", "3/22/22", "Pre-treatment", "control", "Day_0", "control_Day_0" ], [ "39", "Acer-113", "SI-C", "113", "4/20/22", "last day of treatment", "variable", "Day_29", "variable_Day_29" ], [ "40", "Acer-114", "SI-C", "114", "4/20/22", "last day of treatment", "control", "Day_29", "control_Day_29" ], [ "41", "Acer-123", "SI-C", "123", "4/20/22", "last day of treatment", "control", "Day_29", "control_Day_29" ], [ "42", "Acer-124", "BC-8b", "124", "3/22/22", "Pre-treatment", "control", "Day_0", "control_Day_0" ], [ "43", "Acer-137", "BC-8b", "137", "4/20/22", "last day of treatment", "variable", "Day_29", "variable_Day_29" ], [ "44", "Acer-139", "SI-C", "139", "3/22/22", "Pre-treatment", "variable", "Day_0", "variable_Day_0" ], [ "45", "Acer-141", "SI-C", "141", "4/20/22", "last day of treatment", "variable", "Day_29", "variable_Day_29" ], [ "46", "Acer-143", "MB-B", "143", "4/20/22", "last day of treatment", "control", "Day_29", "control_Day_29" ], [ "47", "Acer-146", "BC-8b", "146", "4/20/22", "last day of treatment", "control", "Day_29", "control_Day_29" ], [ "48", "Acer-150", "MB-B", "150", "3/22/22", "Pre-treatment", "control", "Day_0", "control_Day_0" ] ];
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
