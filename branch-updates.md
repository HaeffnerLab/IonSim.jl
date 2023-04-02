# Ions.jl

+ Reorganized the code in ions.jl.
    - This was done for (my subjective opinion of) readability (with go-ahead from Neil). I tried to stick to the format (from top of file to bottom) 
        1) object fields
        2) object definitions 
        3) helper functions for object definitions 
        4) object modifiers
        5) user-exposed functions that take new objects as arguments
        6) base overrides
    But we should probably split this file up at some point.
+ `full_level_structure` renamed `level_structure`
+ `full_transitions` renamed `transitions`
+ Removed `default_sublevels`
+ Removed `sublevel_aliases`