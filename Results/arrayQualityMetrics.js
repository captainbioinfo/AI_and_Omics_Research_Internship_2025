// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false ];
var arrayMetadata    = [ [ "1", "GSM4731674", "T1-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 1", "Tumor" ], [ "2", "GSM4731675", "T2-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 2", "Tumor" ], [ "3", "GSM4731676", "T4-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 4", "Tumor" ], [ "4", "GSM4731677", "T5-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 5", "Tumor" ], [ "5", "GSM4731678", "T6-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 6", "Tumor" ], [ "6", "GSM4731679", "T7-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 7", "Tumor" ], [ "7", "GSM4731680", "T8-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 8", "Tumor" ], [ "8", "GSM4731681", "T9-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 9", "Tumor" ], [ "9", "GSM4731682", "T10-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 10", "Tumor" ], [ "10", "GSM4731683", "T16-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 16", "Tumor" ], [ "11", "GSM4731684", "T17-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 17", "Tumor" ], [ "12", "GSM4731685", "T18-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 18", "Tumor" ], [ "13", "GSM4731686", "T19-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 19", "Tumor" ], [ "14", "GSM4731687", "T20-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 20", "Tumor" ], [ "15", "GSM4731688", "T21-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 21", "Tumor" ], [ "16", "GSM4731689", "T22-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 22", "Tumor" ], [ "17", "GSM4731690", "T23-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 23", "Tumor" ], [ "18", "GSM4731691", "T24-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 24", "Tumor" ], [ "19", "GSM4731692", "T26-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 26", "Tumor" ], [ "20", "GSM4731693", "T27-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 27", "Tumor" ], [ "21", "GSM4731694", "T28-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 28", "Tumor" ], [ "22", "GSM4731696", "T30-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 30", "Tumor" ], [ "23", "GSM4731697", "T31-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 31", "Tumor" ], [ "24", "GSM4731698", "T32-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 32", "Tumor" ], [ "25", "GSM4731699", "T33-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 33", "Tumor" ], [ "26", "GSM4731700", "T34-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 34", "Tumor" ], [ "27", "GSM4731701", "T35-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 35", "Tumor" ], [ "28", "GSM4731702", "T36-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 36", "Tumor" ], [ "29", "GSM4731703", "T37-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 37", "Tumor" ], [ "30", "GSM4731704", "T38-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 38", "Tumor" ], [ "31", "GSM4731705", "T39-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 39", "Tumor" ], [ "32", "GSM4731706", "T40-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 40", "Tumor" ], [ "33", "GSM4731707", "T41-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 41", "Tumor" ], [ "34", "GSM4731708", "T42-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 42", "Tumor" ], [ "35", "GSM4731709", "T43-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 43", "Tumor" ], [ "36", "GSM4731710", "T44-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 44", "Tumor" ], [ "37", "GSM4731711", "T45-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 45", "Tumor" ], [ "38", "GSM4731712", "T46-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 46", "Tumor" ], [ "39", "GSM4731713", "T47-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 47", "Tumor" ], [ "40", "GSM4731714", "T48-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 48", "Tumor" ], [ "41", "GSM4731715", "T49-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 49", "Tumor" ], [ "42", "GSM4731716", "T50-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 50", "Tumor" ], [ "43", "GSM4731717", "T51-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 51", "Tumor" ], [ "44", "GSM4731718", "T52-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 52", "Tumor" ], [ "45", "GSM4731719", "T53-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 53", "Tumor" ], [ "46", "GSM4731720", "T54-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 54", "Tumor" ], [ "47", "GSM4731721", "T55-RNA", "Tumor from colorectal cancer patient", "tissue: Tumor", "Colorectal cancer patient 55", "Tumor" ], [ "48", "GSM4731746", "N1-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 1", "Native tissue" ], [ "49", "GSM4731747", "N2-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 2", "Native tissue" ], [ "50", "GSM4731748", "N4-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 4", "Native tissue" ], [ "51", "GSM4731749", "N5-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 5", "Native tissue" ], [ "52", "GSM4731750", "N6-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 6", "Native tissue" ], [ "53", "GSM4731751", "N7-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 7", "Native tissue" ], [ "54", "GSM4731752", "N8-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 8", "Native tissue" ], [ "55", "GSM4731753", "N9-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 9", "Native tissue" ], [ "56", "GSM4731754", "N10-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 10", "Native tissue" ], [ "57", "GSM4731755", "N16-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 16", "Native tissue" ], [ "58", "GSM4731756", "N17-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 17", "Native tissue" ], [ "59", "GSM4731757", "N18-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 18", "Native tissue" ], [ "60", "GSM4731758", "N19-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 19", "Native tissue" ], [ "61", "GSM4731759", "N20-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 20", "Native tissue" ], [ "62", "GSM4731760", "N21-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 21", "Native tissue" ], [ "63", "GSM4731761", "N22-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 22", "Native tissue" ], [ "64", "GSM4731762", "N23-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 23", "Native tissue" ], [ "65", "GSM4731763", "N24-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 24", "Native tissue" ], [ "66", "GSM4731764", "N26-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 26", "Native tissue" ], [ "67", "GSM4731765", "N27-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 27", "Native tissue" ], [ "68", "GSM4731766", "N28-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 28", "Native tissue" ], [ "69", "GSM4731767", "N29-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 29", "Native tissue" ], [ "70", "GSM4731768", "N30-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 30", "Native tissue" ], [ "71", "GSM4731769", "N31-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 31", "Native tissue" ], [ "72", "GSM4731770", "N32-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 32", "Native tissue" ], [ "73", "GSM4731771", "N33-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 33", "Native tissue" ], [ "74", "GSM4731772", "N34-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 34", "Native tissue" ], [ "75", "GSM4731773", "N35-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 35", "Native tissue" ], [ "76", "GSM4731774", "N36-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 36", "Native tissue" ], [ "77", "GSM4731775", "N37-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 37", "Native tissue" ], [ "78", "GSM4731776", "N38-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 38", "Native tissue" ], [ "79", "GSM4731777", "N39-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 39", "Native tissue" ], [ "80", "GSM4731778", "N40-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 40", "Native tissue" ], [ "81", "GSM4731779", "N41-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 41", "Native tissue" ], [ "82", "GSM4731780", "N42-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 42", "Native tissue" ], [ "83", "GSM4731781", "N43-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 43", "Native tissue" ], [ "84", "GSM4731782", "N44-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 44", "Native tissue" ], [ "85", "GSM4731783", "N45-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 45", "Native tissue" ], [ "86", "GSM4731784", "N46-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 46", "Native tissue" ], [ "87", "GSM4731785", "N47-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 47", "Native tissue" ], [ "88", "GSM4731786", "N48-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 48", "Native tissue" ], [ "89", "GSM4731787", "N49-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 49", "Native tissue" ], [ "90", "GSM4731788", "N50-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 50", "Native tissue" ], [ "91", "GSM4731789", "N51-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 51", "Native tissue" ], [ "92", "GSM4731790", "N52-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 52", "Native tissue" ], [ "93", "GSM4731791", "N53-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 53", "Native tissue" ], [ "94", "GSM4731792", "N54-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 54", "Native tissue" ], [ "95", "GSM4731793", "N55-RNA", "Native tissue from colorectal cancer patient", "tissue: Native tissue", "Colorectal cancer patient 55", "Native tissue" ] ];
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
