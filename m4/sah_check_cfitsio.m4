# SETI_BOINC is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2, or (at your option) any later
# version.

AC_PREREQ([2.54])

AC_DEFUN([SAH_CHECK_CFITSIO],[
  AC_ARG_VAR([CFITSIO],[CFITSIO directory])
  if test -z "$HEAD"
  then
    AC_PATH_PROG(HEAD,head)
  fi
  if test -z "$FIND"
  then
    AC_PATH_PROG(FIND,find)
  fi
  thisdir=`pwd`
  AC_MSG_CHECKING([for CFITSIO])
  cfitsio_search_path=[`echo $CFITSIO cfitsio* ../cfitsio* ../lib/cfitsio* $HOME/lib/cfitsio* /usr/local/cfitsio* /usr/local`]
  for cfitsio_dir in $cfitsio_search_path
  do
    if test -d $cfitsio_dir
    then
      if test -f $cfitsio_dir/libcfitsio.a
      then
        cd $cfitsio_dir
        CFITSIO=`pwd`
	cd $thisdir
	break
      else
        if $FIND $cfitsio_dir -name "lib/libcfitsio.a" 2>/dev/null >/dev/null
	then
	  CFITSIO=`$FIND $cfitsio_dir -name "lib/libcfitsio.a" -print | $HEAD -1 | sed 's/\/lib\/libcfitsio.a//'`         
          cd $CFITSIO
          CFITSIO=`pwd`
	  cd $thisdir
	  break
	fi
      fi
    fi
  done
  if test -n "$CFITSIO" 
  then
    AC_MSG_RESULT($CFITSIO)
  else
    AC_MSG_RESULT(not found)
    no_cfitsio=yes
  fi
])
	

