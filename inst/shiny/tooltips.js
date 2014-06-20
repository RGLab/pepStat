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
      title: "Position Database",
      message: "A peptide collection -- currently, the collections available " +
        "in <tt>pepDat</tt>. Please consult the documentation available for the " +
        "<tt>pepDat</tt> package for more information.",
    },
    
    "slidingMean_width-help": {
      title: "Sliding Mean",
      message: "The width of the sliding window used when smoothing intensities."
    },
    
    "slidingMean_split_by_clade-help": {
      title: "Split by Clade",
      message: "If checked, peptides are smoothed within clades. Otherwise, a peptide at a given position will borrow information from the neighboring peptides as well as ones from other clades around this position."
    },
    
    "makePeptideSet_check_row_order-help": {
      title: "Reduce to a Common Set of Peptides?",
      message: "If checked, the final set of probes is taken to be those with " +
        "ids found in all arrays that were read."
    },
    
    "makeCalls_method-help": {
      title: "Method used to make Positivity Calls",
      message: "<h6>Absolute</h6>" +
        "<p>Intensities exceeding the threshold are labelled as positive.</p>" +
        "<h6>False Discovery Rate</h6>" +
        "<p>A left-tail method is used to generate a threshold, controlling the " +
        "false discovery rate at level <tt>cutoff</tt>.</p>"
    },
    
    "makeCalls_cutoff-help": {
      title: "Positivity Cutoff",
      message: "<h6>Absolute</h6>" +
        "<p>Intensities exceeding the threshold are labelled as positive.</p>" +
        "<h6>False Discovery Rate</h6>" +
        "<p>A left-tail method is used to generate a threshold, controlling the " +
        "false discovery rate at level <tt>cutoff</tt>.</p>"
    },
    
    "makeCalls_group-help": {
      title: "Group by...",
      message: "Select this option to group slides according to a certain " +
        "variable in your data."
    }
      
    
  };
  
  $(".help-icon").each( function(index, value) {
    var id = $(this).attr("id");
    var msg = helpMessages[id];
    var help_message = "<strong>" + msg.title + "</strong><br/><br/>" + msg.message;
    new Opentip("#" + id, help_message, { fixed: true, background: "#FFFFBA", borderColor: "#EAEAEA" });
  })

});
