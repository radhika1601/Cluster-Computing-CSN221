#Add the client computer's IP Addresses to /etc/hosts
cat /etc/hosts;

#This will create a new user named mpiuser 
sudo adduser mpiuser;
#This will let you set password for mpiuser
sudo passwd mpiuser;



#Setting up SSH

#install ssh
sudo dnf install -y openssh-server;
#start ssh service on device
sudo systemctl start sshd.service;
sudo systemctl enable sshd.service;
#login to the created user 'mpiuser'
su - mpiuser

#Generate SSH keys and copy them to the machines' list of authorized_keys
ssh-keygen -t rsa

ssh mpiuser@master mkdir -p .ssh
cat .ssh/id_rsa.pub | ssh mpiuser@master 'cat >> .ssh/authorized_keys'
ssh mpiuser@master 'chmod 700 .ssh'
chmod 640 .ssh/authorized_keys
#This will create a folder in $HOME for mpiuser named cloud that is shared in cluster
mkdir cloud
exit

#Install NFS service
sudo dnf -y install nfs-utils
#mount the shared directory
sudo mount -t nfs master:/home/mpiuser/cloud /home/mpiuser/cloud
#to check the mounted directories
df -h
#Create a systems table in /etc/fstab
echo "Add master:/home/mpiuser/cloud /home/mpiuser/cloud nfs  to /etc/fstab"
sudo vi /etc/fstab

