#! /bin/csh

# Run this with the desired version number as the only operand.  Ex:
# collect2_line 2.27

/usr/local/gcc/bin/../lib/gcc-lib/sparc-sun-solaris2.7/3.0.4/collect2 -V -Y P,/usr/ccs/lib:/usr/lib -Qy -o setiathome_${1}_sparc-sun-solaris2.7 /usr/local/gcc/bin/../lib/gcc-lib/sparc-sun-solaris2.7/3.0.4/crt1.o /usr/local/gcc/bin/../lib/gcc-lib/sparc-sun-solaris2.7/3.0.4/crti.o /usr/ccs/lib/values-Xa.o /usr/local/gcc/bin/../lib/gcc-lib/sparc-sun-solaris2.7/3.0.4/crtbegin.o -L. -L../../boinc/lib -L/usr/local/gcc/bin/../lib/gcc-lib/sparc-sun-solaris2.7/3.0.4 -L/usr/local/gcc/bin/../lib/gcc-lib -L/usr/local/gcc-3.0.4/lib/gcc-lib/sparc-sun-solaris2.7/3.0.4 -L/usr/ccs/bin -L/usr/ccs/lib -L/usr/local/gcc/bin/../lib/gcc-lib/sparc-sun-solaris2.7/3.0.4/../../.. -L/usr/local/gcc-3.0.4/lib/gcc-lib/sparc-sun-solaris2.7/3.0.4/../../.. main.o analyzeFuncs.o analyzeReport.o analyzePoT.o pulsefind.o gaussfit.o lcgamm.o malloc_a.o seti.o seti_header.o timecvt.o s_util.o version.o worker.o chirpfft.o spike.o progress.o ../db/schema_master_client.o ../db/sqlrow_client.o ../db/sqlblob.o ../db/xml_util.o -looura -B static -lstdc++ -B dynamic -lz -lsocket -lnsl -lelf -ldl -laio -lboinc -lm -lgcc_s -lgcc -lc -lgcc_s -lgcc /usr/local/gcc/bin/../lib/gcc-lib/sparc-sun-solaris2.7/3.0.4/crtend.o /usr/local/gcc/bin/../lib/gcc-lib/sparc-sun-solaris2.7/3.0.4/crtn.o
