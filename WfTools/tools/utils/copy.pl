use strict;

my $servers=shift;

defined $servers or $servers=1;

system("unalias -a *");  #to avoid "prompt command"
system("chmod +x *.sh"); #to make executable
system("./stopAllServers.sh");

system("rm server1/output/*");
system("rm -r server1/runs/*");

for(my $i=2;$i<=$servers;$i++){
	system("rm -r server$i");
	system("cp -r server1 server$i");
}

system("./restart.sh");

