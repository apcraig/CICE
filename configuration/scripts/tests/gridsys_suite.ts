# Test         Grid    PEs        Sets    BFB-compare
smoke          gx3     8x2        diag1,run5day
restart        gx3     4x2        debug,diag1
smoke          gbox80  1x1        box2001
smoke          gbox80  1x1        boxslotcyl
smoke	       gbox80  1x1        boxsymn
smoke	       gbox80  1x1        boxsyme
smoke	       gbox80  1x1        boxsymne
smoke	       gbox80  1x1        boxsymn,kmtislands
smoke	       gbox80  1x1        boxsyme,kmtislands
smoke	       gbox80  1x1        boxsymne,kmtislands
smoke	       gbox80  1x1        boxislands_n
smoke	       gbox80  1x1        boxislands_e
smoke	       gbox80  1x1        boxislands_ne
smoke	       gbox80  1x1        boxadv


smoke          gx3     8x2        diag1,run5day,gridcd
restart        gx3     4x2        debug,diag1,gridcd
smoke          gbox80  1x1        box2001,gridcd
smoke          gbox80  1x1        boxslotcyl,gridcd
smoke	       gbox80  1x1        boxsymn,gridcd
smoke	       gbox80  1x1        boxsyme,gridcd
smoke	       gbox80  1x1        boxsymne,gridcd
smoke	       gbox80  1x1        boxsymn,kmtislands,gridcd
smoke	       gbox80  1x1        boxsyme,kmtislands,gridcd
smoke	       gbox80  1x1        boxsymne,kmtislands,gridcd
smoke	       gbox80  1x1        boxislands_n,gridcd
smoke	       gbox80  1x1        boxislands_e,gridcd
smoke	       gbox80  1x1        boxislands_ne,gridcd
smoke	       gbox80  1x1        boxadv,gridcd
