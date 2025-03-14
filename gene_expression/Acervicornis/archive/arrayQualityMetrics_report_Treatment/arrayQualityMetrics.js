// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, true, false, false, false, false, false, true, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false ];
var arrayMetadata    = [ [ "1", "Acer.005", "SI_C", "5", "3/22/22", "Initial" ], [ "2", "Acer.009", "BC_8b", "9", "3/22/22", "Initial" ], [ "3", "Acer.019", "BC_8b", "19", "3/22/22", "Initial" ], [ "4", "Acer.021", "MB_B", "21", "3/22/22", "Initial" ], [ "5", "Acer.022", "BC_8b", "22", "3/22/22", "Initial" ], [ "6", "Acer.029", "SI_C", "29", "3/22/22", "Initial" ], [ "7", "Acer.030", "BC_8b", "30", "4/20/22", "Treated" ], [ "8", "Acer.034", "SI_C", "34", "3/22/22", "Initial" ], [ "9", "Acer.038", "MB_B", "38", "4/20/22", "Treated" ], [ "10", "Acer.041", "MB_B", "41", "4/20/22", "Treated" ], [ "11", "Acer.042", "MB_B", "42", "3/22/22", "Initial" ], [ "12", "Acer.043", "BC_8b", "43", "3/22/22", "Initial" ], [ "13", "Acer.045", "BC_8b", "45", "3/22/22", "Initial" ], [ "14", "Acer.055", "MB_B", "55", "3/22/22", "Initial" ], [ "15", "Acer.056", "MB_B", "56", "4/20/22", "Untreated" ], [ "16", "Acer.057", "MB_B", "57", "4/20/22", "Untreated" ], [ "17", "Acer.059", "MB_B", "59", "4/20/22", "Treated" ], [ "18", "Acer.063", "BC_8b", "63", "3/22/22", "Initial" ], [ "19", "Acer.064", "SI_C", "64", "3/22/22", "Initial" ], [ "20", "Acer.068", "BC_8b", "68", "4/20/22", "Treated" ], [ "21", "Acer.069", "SI_C", "69", "4/20/22", "Treated" ], [ "22", "Acer.072", "BC_8b", "72", "4/20/22", "Untreated" ], [ "23", "Acer.074", "MB_B", "74", "4/20/22", "Treated" ], [ "24", "Acer.078", "SI_C", "78", "3/22/22", "Initial" ], [ "25", "Acer.079", "SI_C", "79", "3/22/22", "Initial" ], [ "26", "Acer.082", "SI_C", "82", "4/20/22", "Untreated" ], [ "27", "Acer.083", "MB_B", "83", "3/22/22", "Initial" ], [ "28", "Acer.089", "BC_8b", "89", "4/20/22", "Untreated" ], [ "29", "Acer.091", "BC_8b", "91", "4/20/22", "Treated" ], [ "30", "Acer.094", "MB_B", "94", "3/22/22", "Initial" ], [ "31", "Acer.095", "SI_C", "95", "3/22/22", "Initial" ], [ "32", "Acer.096", "SI_C", "96", "4/20/22", "Treated" ], [ "33", "Acer.097", "MB_B", "97", "4/20/22", "Untreated" ], [ "34", "Acer.102", "MB_B", "102", "3/22/22", "Initial" ], [ "35", "Acer.106", "SI_C", "106", "4/20/22", "Untreated" ], [ "36", "Acer.108", "BC_8b", "108", "4/20/22", "Untreated" ], [ "37", "Acer.110", "BC_8b", "110", "3/22/22", "Initial" ], [ "38", "Acer.112", "MB_B", "112", "3/22/22", "Initial" ], [ "39", "Acer.113", "SI_C", "113", "4/20/22", "Treated" ], [ "40", "Acer.114", "SI_C", "114", "4/20/22", "Untreated" ], [ "41", "Acer.123", "SI_C", "123", "4/20/22", "Untreated" ], [ "42", "Acer.124", "BC_8b", "124", "3/22/22", "Initial" ], [ "43", "Acer.137", "BC_8b", "137", "4/20/22", "Treated" ], [ "44", "Acer.139", "SI_C", "139", "3/22/22", "Initial" ], [ "45", "Acer.141", "SI_C", "141", "4/20/22", "Treated" ], [ "46", "Acer.143", "MB_B", "143", "4/20/22", "Untreated" ], [ "47", "Acer.146", "BC_8b", "146", "4/20/22", "Untreated" ], [ "48", "Acer.150", "MB_B", "150", "3/22/22", "Initial" ] ];
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
