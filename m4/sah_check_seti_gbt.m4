# SETI_BOINC is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2, or (at your option) any later
# version.

AC_PREREQ([2.54])

AC_DEFUN([SAH_CHECK_SETI_GBT],[
  AC_ARG_VAR([SETI_GBT],[SETI_GBT directory])
  if test -z "$HEAD"
  then
    AC_PATH_PROG(HEAD,head)
  fi
  if test -z "$FIND"
  then
    AC_PATH_PROG(FIND,find)
  fi
  thisdir=`pwd`
  AC_MSG_CHECKING([for SETI_GBT])
  seti_gbt_search_path=[`echo $SETI_GBT seti_gbt* ../seti_gbt* ../lib/seti_gbt* $HOME/lib/seti_gbt* /usr/local/seti_gbt* /usr/local`]
  for seti_gbt_dir in $seti_gbt_search_path
  do
    if test -d $seti_gbt_dir/lib
    then
      if test -f $seti_gbt_dir/lib/libsetigbt.so
      then
        cd $seti_gbt_dir
        SETI_GBT=`pwd`
	cd $thisdir
	break
      else
        if $FIND $seti_gbt_dir -name "lib/libsetigbt.so" 2>/dev/null >/dev/null
	then
	  SETI_GBT=`$FIND $seti_gbt_dir -name "lib/libsetigbt.so" -print | $HEAD -1 | sed 's/\/lib\/libsetigbt.so//'`         
          cd $SETI_GBT
          SETI_GBT=`pwd`
	  cd $thisdir
	  break
	fi
      fi
    fi
  done
  if test -n "$SETI_GBT" 
  then
    AC_MSG_RESULT($SETI_GBT)
  else
    AC_MSG_RESULT(not found)
    no_seti_gbt=yes
  fi
])
