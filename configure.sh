#This unzips the tar file
tar -xzf mpich-3.3.1.tar.gz

#This enters the folder created on unzipping the tar file
cd  mpich-3.3.1

#This will configure the prerequisites for cluster formation
./configure --disable-fortran
make
sudo make install

#If configuration is successful this command will give the version
mpiexec --version
