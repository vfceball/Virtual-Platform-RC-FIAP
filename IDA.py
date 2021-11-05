#set pid [getPID]
#set np [getNP]
from datetime import datetime
from ReadRecord import ReadRecord
from getSaT import getSaT
from load_PEERNGA_record import load_PEERNGA_record

start_time = datetime.now()

firstInt = float(self.ui.firstInt.text())  # This is the first intensity to run the elastic run (e.g. 0.05g)
incrStep = float(self.ui.incrStep.text()) # This is the increment used during the hunting phase (e.g. 0.10g)
maxRuns = float(self.ui.maxRuns.text()) # This is the maximum number of runs to use (e.g. 20)
dCap = float(self.ui.dCap.text()) # Drift capacity
xi = float(self.ui.xi.text()) # Damping


firstInt = 0.05  # This is the first intensity to run the elastic run (e.g. 0.05g)
incrStep = 0.10  # This is the increment used during the hunting phase (e.g. 0.10g)
maxRuns = 1  # This is the maximum number of runs to use (e.g. 20)
IMtype = 1  # Intensity measure with which to conduct the IDA.
dCap = 6.  # Drift capacity in %
xi = 0.05  # Damping
outsdir = "IDA" # Where to print the outputs
mdlfile = "mod_nl.tcl"
pflag = 1
IDA_HTF_error_log = []
# Open an error file that will log the IDA_HTF errors
print('^^^^^^^^ STARTING IDA HTF ^^^^^^^^')
IDA_HTF_error_log.append('^^^^^^^^ STARTING IDA HTF ^^^^^^^^')

# Get the ground motion set information
nmsfile = self.ui.nmsfile.text()

with open('Earthquakes/' + nmsfile + '.txt', "r") as archivo:
    eqnms_list = [linea.rstrip() for linea in archivo]

nrecs = len(eqnms_list)
for ind in range(nrecs):
    # EQname = 'Earthquakes/' + eqnms_list[ind] + '.at2'
    # acc, dt, npts, eqname = load_PEERNGA_record(EQname)

    EQname = 'Earthquakes/' + eqnms_list[ind]
    inFile = EQname + '.at2'
    outFile = EQname + '.g4'
    dt, npts = ReadRecord(inFile, outFile)

    IMtype = self.ui.comboBoxIM.currentText()
    if IMtype == 'PGA':
        Sa, Sv, Sd, pga, amax = getSaT(outFile, dt, 0.0, xi, npts)  # Get the PGA of the record in the X direction
        IMgeomean = pga

    elif IMtype == 'Sa(T1)':
        Sa, Sv, Sd, pga, amax = getSaT(outFile, dt, T1m, xi, npts)  # Get the PGA of the record in the X direction
        IMgeomean = Sa
        print('T1m =', T1m, 'Sa =', Sa, 'amax', amax)

# set eqnms_listX     [read [open "Sismos/nmsfileX5.txt" "r"]];
# set dts_list        [read [open "Sismos/dtsfile5.txt" "r"]];
# set durs_list       [read [open "Sismos/dursfile5.txt" "r"]];
# set nrecs           [llength $dts_list];
# set g 9.81;
#
# for {set i 1} {$i <= $nrecs} {incr i 1} {
#   #if {[expr ($i-1) % $np] == $pid} {
#     set IM_log [open $outsdir/IM_${i}.txt "w"];
#
#     # Load the info about the record
#     set EQnameX [lindex $eqnms_listX $i-1];     # Get the name of the record1
#     set dt  [lindex $dts_list $i-1];    # Current dt
#     set dur [lindex $durs_list $i-1];   # Current duration
#     puts "dt=$dt dur=$dur"
#
#     # Establish the IM
#     if {$IMtype==1} {
#         # IM is PGA
#         # Now get the spectral ordinates
#         getSaT Sismos/$EQnameX.g3 [lindex $dts_list $i-1] 0.0 $xi; # Get the PGA of the record in the X direction
#         set IMgeomean $pga
#     } elseif {$IMtype==2} {
#         # IM is Sa at a specified period
#         # Need to get the conditioning period
#         # Now get the spectral ordinates
#         getSaT Sismos/$EQnameX.g3 [lindex $dts_list $i-1] $T1m $xi; # Get the Sa(T1,5%) of the record in the X direction
#         set IMgeomean $Sa
#     } elseif {$IMtype==3} {
#         # IM is Sa,gm at a specified periods
#         # Period list [T2m, min[(T2m+T1m)/2, 1.5T2m], T1m, 1.5T1m]
#         set p_list {};
#         lappend p_list $T2m;
#         if {[expr 0.5*($T2m+$T1m)]<[expr 1.5*$T2m]} {
#             lappend p_list [expr 0.5*($T2m+$T1m)];
#         } else {
#             lappend p_list [expr 1.5*$T2m];
#         };
#         lappend p_list $T1m;
#         lappend p_list [expr 1.5*$T1m];
#
#         # Get Spectral values at each
#         set Sa_listX {};
#         for {set pn 1} {$pn<=[llength $p_list]} {incr pn 1} {
#         getSaT Sismos/$EQnameX.g3 [lindex $dts_list $i-1] [lindex $p_list $pn-1] $xi; # Get the Sa(T1,5%) of the record in the X direction
#             lappend Sa_listX $Sa
#         };
#         # Get the geometric mean of these
#         set IMgeomean [expr pow([lindex $Sa_listX 0]*[lindex $Sa_listX 1]*[lindex $Sa_listX 2]*[lindex $Sa_listX 3],1/4.0)];
#     } elseif {$IMtype==4} {
#         # IM is Sa,gm at a specified periods
#         # Period list [T2m, min[(T2m+T1m)/2, 1.5T2m], T1m, 1.5T1m, 2T1m]
#         set p_list {};
#         lappend p_list $T2m;
#         if {[expr 0.5*($T2m+$T1m)]<[expr 1.5*$T2m]} {
#             lappend p_list [expr 0.5*($T2m+$T1m)];
#         } else {
#             lappend p_list [expr 1.5*$T2m];
#         };
#         lappend p_list $T1m;
#         lappend p_list [expr 1.5*$T1m];
#         lappend p_list [expr 2.0*$T1m];
#
#         # Get Spectral values at each
#         set Sa_listX {};
#         for {set pn 1} {$pn<=[llength $p_list]} {incr pn 1} {
#         getSaT Sismos/$EQnameX.g3 [lindex $dts_list $i-1] [lindex $p_list $pn-1] $xi; # Get the Sa(T1,5%) of the record in the X direction
#             lappend Sa_listX $Sa
#         };
#
#         # Get the geometric mean of these
#         set IMgeomean [expr pow([lindex $Sa_listX 0]*[lindex $Sa_listX 1]*[lindex $Sa_listX 2]*[lindex $Sa_listX 3]*[lindex $Sa_listX 4],1/5.0)];
#     } elseif {$IMtype==5} {
#         # IM is Sa,gm at a specified periods
#         # Period list [T2m, min[(T2m+T1m)/2, 1.5T2m], T1m, 1.5T1m, 2T1m]
#         set p_list {};
#         lappend p_list $T2m;
#         if {[expr 0.5*($T2m+$T1m)]<[expr 1.5*$T2m]} {
#             lappend p_list [expr 0.5*($T2m+$T1m)];
#         } else {
#             lappend p_list [expr 1.5*$T2m];
#         };
#         lappend p_list $T1m;
#         lappend p_list [expr 1.5*$T1m];
#         lappend p_list [expr 2.0*$T1m];
#
#         # Get Spectral values at each
#         set Sa_listX {};
#         for {set pn 1} {$pn<=[llength $p_list]} {incr pn 1} {
#         getSaT Sismos/$EQnameX.g3 [lindex $dts_list $i-1] [lindex $p_list $pn-1] $xi; # Get the Sa(T1,5%) of the record in the X direction
#             lappend Sa_listX $Sa
#         };
#
#         # Get the geometric mean of these
#         set SagmX [expr pow([lindex $Sa_listX 0]*[lindex $Sa_listX 1]*[lindex $Sa_listX 2]*[lindex $Sa_listX 3]*[lindex $Sa_listX 4],1/5.0)];
#
#         # Get the time between 5 and 75% of Arias Intensity
#         Arias Sismos/$EQnameX.g3 [lindex $dts_list $i-1] 0.05 0.75;
#         set t12X $t12;
#
#         # Get the IMs in the two directions
#         set IMgeomean [expr $SagmX*pow($t12X,0.2)];
#     };
#     puts "IMgeomean =$IMgeomean "
#     # Set up the initial indices for HTF
#     set j   1;
#     set IM  {};     # Initialise the list of IM used for printing
#     set IMlist  {};     # This is just a list that will be used in filling
#     set hFlag   1;      # Hunting flag (1 for when we're hunting)
#     set tFlag   0;      # Tracing flag (0 at first)
#     set fFlag   0;      # Filling flag (0 at first)
#     puts "j=$j maxRuns=$maxRuns"
#
#     set Tabla_IDA [open "$outsdir/Tabla_IDA_Record${i}.txt" "w"]
#
#     while {$j<=$ maxRuns} {
#         # As long as the hunting flag is 1, meaning we havent reached a collapse
#         if {$hFlag==1} {
#             # Determine the intensity to run at during the hunting (Andiamo a cacciare!)
#             if {$j==1} {
#                 lappend IM $firstInt;
#             } else {
#                 lappend IM [expr [lindex $IM $j-2]+($j-1)*$incrStep]; # Ramp it up!
#             };
#             # Determine the scale factor that needs to be applied to the record
#             set sfX     [expr [lindex $IM $j-1]/$IMgeomean*$g];
#             set     run "Record${i}_Run${j}";           # This is a tag that outputs will be labelled with
#             set     log [open $outsdir/log_IDA_${run}.txt "w"];
#
#             # The hunting intensity has been determined, now we can analyse
#             source  $mdlfile
#             source Exci_pattern.tcl
#             if {$pflag>0} {puts [format "Record:$i  Run:$j IM:%.3f" [lindex $IM $j-1] ]}
#             runNRHA3D $dt $dur $dCap $tNode $bNode $log $pflag $inip $inic $ncols
#             close $log
#             set Sa_IDA [expr $maxVB/$Wtotal]
#             if {$cIndex==0} {puts $Tabla_IDA "[lindex $IM $j-1] $mdrft $Sa_IDA"}
#             incr j 1;
#
#             # Check the hunted run for collapse
#             if {$cIndex>0} {
#                 set     hFlag   0;  # Stop hunting
#                 set     tFlag   1;  # Start tracing
#                 incr    j       -1; # Reduce by 1 because j was increased at end of hunting and we want to redo that point
#                 set     jhunt   $j; # The value of j we hunted to
#                 if {$jhunt==2} {puts $error_log "WARNING: ${run} - Collapsed achieved on first increment, reduce increment..."};
#             } else {
#                 puts $IM_log [format "%.3f" [lindex $IM $j-2]]; # j-2 because we've already increased j, but need to know if collapsed
#             };
#             wipe;
#         }; # Close hunting
#
#         # When the first collapse is reached, we start tracing between last convergence and the first collapse
#         if {$tFlag==1} {
#             # The first phase is to trace the last DeltaIM to get within the resolution
#             if {$j==$jhunt} {
#                 set firstC  [lindex $IM $j-1];              # This is the IM of the hunting collapse
#                 set IM  [lreplace $IM $j-1 $j-1];       # Remove that value of IM from the array (it's already been appended)
#             };
#             set diff    [expr $firstC-[lindex $IM $j-2]];   # Determine the difference between the hunting's noncollapse and collapse IM
#             set inctr   [expr 0.20*$diff]; # Take 0.2 of the difference
#             if {$inctr<0.05} {set inctr 0.025}; # Place a lower threshold on the increment so it doesnt start tracing too fine
#             set IMtr    [expr [lindex $IM $j-2]+$inctr];    # Calculate new tracing IM, which is previous noncollapse plus increment
#             lappend IM $IMtr
#             set sfX     [expr [lindex $IM $j-1]/$IMgeomean*$g];
#             puts $IM_log  [format "%.3f" $IMtr];
#             set     run "Record${i}_Run${j}";           # This is a tag that outputs will be labelled with
#             set     log [open $outsdir/log_IDA_${run}.txt "w"];
#
#             # The trace intensity has been determined, now we can analyse
#             source  $mdlfile
#             source Exci_pattern.tcl
#             if {$pflag>0} {puts [format "Record:$i  Run:$j IM:%.3f" $IMtr ]}
#             runNRHA3D $dt $dur $dCap $tNode $bNode $log $pflag $inip $inic $ncols
#
#             close $log
#             set Sa_IDA [expr $maxVB/$Wtotal]
#             if {$cIndex==0} {puts $Tabla_IDA "$IMtr $mdrft $Sa_IDA"}
#             if {$cIndex>0} {
#                 # Not sure if this is the best way, to just trace back up to collapse again
#                 set tFlag 0; # Stop tracing
#                 set fFlag 1; # Start filling
#                 set jtrace  $j; # The value of j we traced to
#                 set IMlist $IM; # Get the list of IMs
#                 if {$j==$jhunt} {
#                     # This means the first trace collapsed, should reduce the increment
#                     puts $error_log "WARNING: ${run} - First trace for collapse resulted in collapse..."
#                 };
#             };
#             incr j 1;
#             wipe;
#         }; # Close the tracing
#
#         # When the required resolution is reached, we start filling
#         if {$fFlag==1} {
#             # Reorder the list so we can account for filled runs
#             set IMlist [lsort -real $IMlist];
#
#             # Determine the biggest gap in IM for the hunted runs
#             set gap 0.0;
#             # We go to the end of the list minus 1 because, if not we would be filling between a noncollapsing and a collapsing run,
#             # for which we are not sure if that filling run would be a non collapse - In short, does away with collapsing fills
#             for {set ii 1} {$ii<[expr [llength $IMlist]-1]} {incr ii 1} {
#                 set temp [expr [lindex $IMlist $ii]-[lindex $IMlist $ii-1]];    # Find the running gap of hunted runs
#                 if {$temp>$gap} {
#                     set gap $temp
#                     set IMfil [expr [lindex $IMlist $ii-1]+$gap/2];     # Determine new filling IM
#                 };  # Update to maximum gap
#             };
#
#             lappend IM $IMfil
#             lappend IMlist $IMfil
#             set sfX     [expr [lindex $IM $j-1]/$IMgeomean*$g];
#             puts $IM_log  [format "%.3f" $IMfil];
#             set     run "Record${i}_Run${j}";           # This is a tag that outputs will be labelled with
#             set     log [open $outsdir/log_IDA_${run}.txt "w"];
#
#             # The trace intensity has been determined, now we can analyse
#             source  $mdlfile
#             source Exci_pattern.tcl
#             if {$pflag>0} {puts [format "Record:$i  Run:$j IM:%.3f" $IMfil ]}
#             runNRHA3D $dt $dur $dCap $tNode $bNode $log $pflag $inip $inic $ncols
#
#             close $log
#             set Sa_IDA [expr $maxVB/$Wtotal]
#             if {$cIndex==0} {puts $Tabla_IDA "$IMfil $mdrft $Sa_IDA"}
#             incr j 1;
#             wipe;
#         }; # Close the filling
#
#         # Wrap it up and finish
#         if {$j==$maxRuns && $hFlag==1} {
#             puts $error_log "WARNING: ${run} - Collapse not achieved, increase increment or number of runs..."
#         };
#         if {$j==$maxRuns && $fFlag==0} {
#             puts $error_log "WARNING: ${run} - No filling, algorithm still tracing for collapse (reduce increment & increase runs)..."
#         };
#         wipe;
#     }; # Close the maxRuns while loop
#     close $IM_log
#     close $Tabla_IDA
#  #};
# };
# puts "^^^^^^^^ FINISHED IDA HTF ^^^^^^^^"
# puts $error_log "^^^^^^^^ FINISHED IDA HTF ^^^^^^^^"
# close $error_log


time_elapsed = datetime.now() - start_time
print('Time elapsed (hh:mm:ss.ms) {}'.format(time_elapsed))
