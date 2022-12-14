  
  
  ############################################################################
  ### RP PACKAGE - GENERIC METHODS 
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     01.07.2019
  # First version:    01.07.2019
  # --------------------------------------------------------------------------
  
  # METHODS:
  # show
  
  # --------------------------------------------------------------------------
  setMethod(f = "show",
            signature = "VBODY",
            definition = function(x)
            {
              slot_names <- slotNames(x)
              for (i in seq(along = slot_names) ) {
                tmp <- slot(x, slot_names[i])
                cat( "@", slot_names[i])
                show(tmp)
              }
            })
  
 
  
  
  
  
  
  
  
  
  