mkdir -p Libraries
cp *.sh Libraries
cp *.inc Libraries

cd Libraries

sh INSTALL_Metis.sh
sh INSTALL_SCOTCH.sh
sh INSTALL_OpenBlas.sh
sh INSTALL_MUMPS.sh
