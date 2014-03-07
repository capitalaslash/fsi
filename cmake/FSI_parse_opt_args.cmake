# macro to parse optional arguments to a function/macro.
# args is a list of arguments to be checked, that will be set to the specified
# value if found in the optional arguments
macro(FSI_parse_opt_args args)
  #message(STATUS "args = ${args}\nARGN = ${ARGN}")

  list(APPEND opt_arg ${args})
  set(found_arg_name)
  foreach(arg ${ARGN})
    if(found_arg_name)
      #message(STATUS "FOUND ${found_arg_name} ${arg}")
      set(${found_arg_name} ${arg})
      set(found_arg_name)
    else()
      #message(STATUS "SEARCH ${opt_arg} ${arg}")
      list(FIND opt_arg "${arg}" is_arg)
      if(is_arg GREATER -1)
        set(found_arg_name ${arg})
      endif()
    endif()
  endforeach()

endmacro()

