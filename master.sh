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

ssh mpiuser@client mkdir -p .ssh
cat .ssh/id_rsa.pub | ssh mpiuser@client 'cat >> .ssh/authorized_keys'
ssh mpiuser@client 'chmod 700 .ssh'
chmod 640 .ssh/authorized_keys
#This will create a folder in $HOME for mpiuser named cloud that is shared in cluster
mkdir cloud
exit

#Install NFS service 
sudo dnf -y install nfs-utils
#add the n=given folder to /etc/exports
echo "ADD /home/mpiuser/cloud *(rw,sync,no_root_squash,no_subtree_check) to /etc/exports"
#This will open the /etc/exports folder in vi editor
sudo vi /etc/exports
#This should be run everytime ther is a change in /etc/exports
exportfs -a
#this will start the nfs service
sudo systemctl start rpcbind nfs-server
sudo systemctl enable rpcbind nfs-server


echo "ssh and nfs have been setUp\n add all cluster IPs to /etc/hosts"
