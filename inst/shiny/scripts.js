$( function() {
  
  $("#R_input").keypress( function(e) {
    if (e.which == 13) {
      $("#R_send").click();
    }
  });
  
});