use strict;

my $conf = "/etc/ochem/";
my $file = "$conf/tasks";
my $MINIMAL_MEMORY = 1024;

my $test = 0;

my (%SERVERS,$JAVA,$METASERVER,$MEMORY,$COUNT,$GPU,$SLEEP);

$SLEEP = 5;

updateTasks($file);
my $tasks = initTasks();

for my $SERVER (keys %SERVERS){
    my @tasks = keys %{$tasks->{$SERVER}};
    print "using TASK=$SERVER JAVA=$JAVA->{$SERVER} MEMORY=$MEMORY->{$SERVER} METASERVER=$METASERVER->{$SERVER} and SERVERS:$SERVERS{$SERVER} TASKS: @tasks \n";
}

for(my $count=0; 1; $count++){
    sleep($SLEEP);
    updateTasks($file);
    $tasks =updateMemoryTasks();
    my @tasks = keys %{$tasks};
    for my $SERVER (keys %SERVERS){
        my $server1 = "server$SERVER"."0";
        my $servers = delServers($SERVER);
        -f "$server1/tests.txt" or next; # do not start untill tests are finished
        $test && print "tests are finished\n";
        my $new = checkTasks($SERVER);
        defined $new or next;
        canRUN($SERVER) or next;
        $test && print "can run on $SERVER\n";
        my $start = startServer($SERVER,$servers);
        $start or next;
        $test && print "started server $SERVER for $new\n";
    }
    if($count == 0){
        system("chmod +x *.sh; ./restart.sh 2> /dev/null &");
    }
    
}

sub updateTasks(){
    my $ file = shift;
    open(I,$file) or die("File specifying server configuration is not available $file\n");
    while(<I>){
        my @a = split;
        scalar @a == 2 or next;
        $SERVERS{$a[0]} = $a[1];
    }close(I);
}

sub canRUN(){
    my $SERVER = shift;
    defined $GPU->{$SERVER} or return 1; # can always run if GPU is not used
    
    open(I,"nvidia-smi --query-gpu=memory.free --format=csv -i $GPU->{$SERVER} |");
    
    my $memory;
    
    while(<I>){
        /memory/ && next;
        my @a = split;
        $memory = $a[0];
        $test && print "GPU $GPU->{$SERVER} memory free: $memory\n";
        last;
    }
    close(I);
    return $memory > 2000  ? 1:0 ; # at least 2MB of free GPU memory
}


sub delServers(){
    my $SERVER = shift;
    open(I,"ps -wax | grep java |");
    my @servers;
    while(<I>){
        /server$SERVER(\d+)/ or next;
        $servers[$1]=1;
    }
    close(I);
    for(my $i = 1; $i<=$SERVERS{$SERVER};$i++){ # at least 10 servers above the specified number
        $servers[$i] && next; # running
        system("mv server$SERVER$i/runs/* server$SERVER"."0/runs/ >/dev/null 2>/dev/null"); # to have a possibility to debug failed tasks
        system("rm -r server$SERVER$i 2> /dev/null");
        # killing remaining python processes, if any
        open(I,"ps -wax |");
        while(<I>){
            /server$SERVER$i\s+/ or /server$SERVER$i\// or next;
            my @a = split;
            $a[0] = $a[0]>0?$a[0]:$a[1];
            system("kill -9 $a[0]");
        }
        close(I);
    }
    return \@servers;
}

sub startServer(){
    my ($SERVER,$servers) = @_;
    my $server1 = "server$SERVER"."0";
    for(my $i = 1;$i <= $SERVERS{$SERVER}; $i++){
        my $datestring = localtime();
        $test && print "$datestring starting $i out of $SERVERS{$SERVER} busy: $servers->[$i]\n";
        $servers->[$i] && next;
        ##sleep($SLEEP->{$SERVER});
        my $server = "server$SERVER$i";
        my $mem = $MEMORY->{$SERVER};
        $mem = $mem > $MINIMAL_MEMORY ? $mem : $MINIMAL_MEMORY;
        print "$datestring starting $server with $SERVER and gpu: $GPU->{$SERVER}\n";
        system("cp -r $server1 $server");
        system("rm $server/output/*");
        system("rm -rf $server/runs/*");
        system("perl -pi -e 's|<processId>.+<\/processId>|<maximumIdleTime>150<\/maximumIdleTime>|' $server/version.xml 2> /dev/null");
        #system("perl -pi -e 's|current-version=\".+\">|current-version=\"local\">|' $server/version.xml 2> /dev/null"); # Does not make sence since such servers are not updated after the release of new code
        system("perl -pi -e 's|<memoryLimit>.+<\/memoryLimit>|<memoryLimit>$mem<\/memoryLimit>|' $server/version.xml 2> /dev/null");
        $mem = defined $GPU->{$SERVER}?$mem*5:$mem; # give 5 times more memory for GPU tasks
        system("cd $server; $JAVA->{$SERVER}/bin/java -Xmx$mem"."m -XX:MinHeapFreeRatio=25 -XX:MaxHeapFreeRatio=40 -cp lib/*:lib/activation-1.1.jar qspr.metaserver.serv.MultiServer /etc/cs_servers/$server >o 2>>o &");
        return 1;
    }
    return 0;
}

sub initTasks(){
    my $tasks;
    for my $SERVER (keys %SERVERS){
        my $server1 ="server$SERVER"."0";
        system("cp -r server1 $server1"); # last one will contribute linked files
        my $file = "$conf/$SERVER.xml";
        open(I,$file) or die("cannot open $file\n");
        system("cp $file $server1/version.xml");
        while(<I>){
            if(/<java-home>(\S+)<\/java-home>/){
                $JAVA->{$SERVER} = $1;
            }
            if(/<gpuCard>(\S+)<\/gpuCard>/){
                $GPU->{$SERVER} = $1;
            }
            if(/<metaserverURL>(\S+)<\/metaserverURL>/){
                $METASERVER->{$SERVER} = $1;
            }
            if(/<memoryLimit>(\S+)<\/memoryLimit>/){
                $MEMORY->{$SERVER} = $1;
            }
            if(/<application>(\S+)<\/application>/){
                $tasks->{$SERVER}->{lc($1)}=1;
                lc($1) eq "corinalocal" && do {$tasks->{$SERVER}->{"corina"}=1;};
            }
            #if(/<sleepTime>(\S+)<\/sleepTime>/){
            #    $tasks->{$SLEEP}->{lc($1)}=1;
            #}
        }
        close(I);
        system("perl -pi -e 's|<memoryLimit>.+<\/memoryLimit>|<memoryLimit>512<\/memoryLimit>|' $server1/version.xml 2> /dev/null");
        
    }
    return $tasks;
}

sub updateMemoryTasks(){
    my $tasks;
    for my $SERVER (keys %SERVERS){
        my $file = "$conf/$SERVER.xml";
        open(I,$file) or die("cannot open $file\n");
        while(<I>){
            if(/<java-home>(\S+)<\/java-home>/){
                $JAVA->{$SERVER} = $1;
            }
            if(/<gpuCard>(\S+)<\/gpuCard>/){
                $GPU->{$SERVER} = $1;
            }
            if(/<metaserverURL>(\S+)<\/metaserverURL>/){
                $METASERVER->{$SERVER} = $1;
            }
            if(/<memoryLimit>(\S+)<\/memoryLimit>/){
                $MEMORY->{$SERVER} = $1;
            }
            if(/<application>(\S+)<\/application>/){
                $tasks->{$SERVER}->{lc($1)}=1;
                lc($1) eq "corinalocal" && do {$tasks->{$SERVER}->{"corina"}=1;};
            }
        }
        close(I);
    }
    return $tasks;
}


sub checkTasks{
    my $SERVER = shift;
    my $datestring = localtime();
    $test && print("$datestring trying $SERVER\n");
    open(I,"wget --no-check-certificate --proxy=off --timeout=$SLEEP -q -O - $METASERVER->{$SERVER}?action=queue |");
    $test && print("$datestring got $SERVER\n");
    my $found;
    $COUNT->{$SERVER}=0;
    while(<I>){
        chomp;
        length or next;
        $test && print "$METASERVER->{$SERVER}?action=queue : $_\n";
        my @a=split(/,/,lc($_));
        if($tasks->{$SERVER}->{$a[0]} && $a[1]*1024 <= $MEMORY->{$SERVER}){
            $found = uc($a[0]);
            $COUNT->{$SERVER} = $a[2]; # number of tasks
        }
        defined $found && do{ $test && print "found: $found\n"; last;};
    }close(I);
    return $found;
}
