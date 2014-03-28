$( function() {
  
  var tt = function(x) {
    return "<tt>" + x + "</tt>";
  }
  
  var helpMessages = {
    
    "gpr_files-help": {
      title: "Upload a Set of GenePix Results Files",
      message: "Select a set of GenePix results (.gpr) files you wish to " +
      "analyze."
    },
    
    "mapping_file-help": {
      title: "Upload a Mapping File",
      message: "This should be a tabular, comma-separated " +
        "file with at least three columns, labeled " +
        "<ul>" +
        "<li><tt>filename</tt></li>" +
        "<li><tt>ptid</tt></li>" +
        "<li><tt>visit</tt></li>" +
        "</ul>" +
        "For each <tt>ptid</tt> (Patient ID), the visit column should have at least one " +
        "<tt>pre</tt> and one <tt>post</tt> sample. " +
        "<br/><br/>" +
        "Any additional columns will be kept in the resulting " +
        "peptide set and can be used later on for grouping."
    },
    
    "summarizePeptides_summary-help": {
      title: "Summary Method",
      message: "The method used for merging replicates."
    },
    
    "summarizePeptides_position-help": {
      title: "Additional Information",
      message: "A peptide collection -- currently, the collections available " +
        "in <tt>PEP.db</tt>. Please consult the documentation available for the " +
        "<tt>PEP.db</tt> package for more information.",
    },
    
    "slidingMean_width-help": {
      title: "Sliding Mean",
      message: "The width of the sliding window used when smoothing intensities."
    },
    
    "slidingMean_split_by_space-help": {
      title: "Smooth over Space",
      message: "If checked, the peptides will be smoothed by space. " +
        "When checked, peptides are smoothed within groups defined by the space of the " +
        "<tt>RangedData</tt> object occupying the <tt>FeatureRange</tt> slot of the " +
        "constructed <tt>peptideSet</tt>." +
        "<br/><br/>" +
        "This is useful if, for example, spaces indicate peptides drawn from separate " +
        "proteins that should not be smoothed together."
    },
    
    "makePeptideSet_check_row_order-help": {
      title: "Reduce to a Common Set of Peptides?",
      message: "If checked, the final set of probes is taken to be those with " +
        "ids found in all arrays that were read."
    }
      
    
  };
  
  $(".help-icon").each( function(index, value) {
    var id = $(this).attr("id");
    console.log(id);
    var msg = helpMessages[id];
    var help_message = "<strong>" + msg.title + "</strong><br/><br/>" + msg.message;
    console.log(msg.title);
    console.log(msg.message);
    new Opentip("#" + id, help_message, { fixed: true });
  })
  
  $("#R_input").keypress( function(e) {
    if (e.which == 13) {
      $("#R_send").click();
    }
  });
  
});