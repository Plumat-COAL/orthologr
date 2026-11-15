#!/usr/bin/perl
#
# pal2nal_batch.pl  -  PAL2NAL for many gene pairs in one Perl process (orthologr)
#
#   perl pal2nal_batch.pl <pal2nal.pl> <jobs> <status>
#
#   <jobs>   : one line per gene pair, fields separated by tabs: the output file, then the arguments of pal2nal.pl
#              (e.g. pep.aln, nuc.fasta, -output, fasta)
#   <status> : written at the end, one line per job: the wait status of the job ($?: 0 when PAL2NAL exited with
#              status 0, otherwise exit status * 256 + signal), or -1 when no process could be started for it
#
# pal2nal.pl is read and compiled once. Every job runs in a child process of its own (fork), which starts from the
# compiled but never run script: no variable of one gene pair is seen by another (PAL2NAL keeps its state in
# global variables), and exit / die end the job as they end 'perl pal2nal.pl <arguments> > <output file>'. The
# standard output of the job goes to its output file; its standard error is discarded. The caller checks the
# status and the output of every job and computes a job that did not succeed again with 'perl pal2nal.pl'.
# (No 'use strict' here: pal2nal.pl is compiled in the scope of this file and uses global variables.)

$batch_script = $ARGV[0];
$batch_jobfile = $ARGV[1];
$batch_statusfile = $ARGV[2];
if (!defined($batch_statusfile)) {
    print STDERR "usage: perl pal2nal_batch.pl <pal2nal.pl> <jobs> <status>\n";
    exit 2;
}

open(BATCH_SRC, "<", $batch_script) || die "Can't open $batch_script: $!\n";
{
    local $/;
    $batch_src = <BATCH_SRC>;
}
close(BATCH_SRC);

# the script as the body of a subroutine; '#line' keeps its file name and line numbers in messages
$batch_main = eval("sub {\n#line 1 \"$batch_script\"\n" . $batch_src . "\n;}");
if (!defined($batch_main)) {
    print STDERR "pal2nal.pl could not be compiled: $@\n";
    exit 3;
}
undef($batch_src);

open(BATCH_JOBS, "<", $batch_jobfile) || die "Can't open $batch_jobfile: $!\n";
@batch_jobs = <BATCH_JOBS>;
close(BATCH_JOBS);

@batch_status = ();
foreach $batch_job (@batch_jobs) {
    chomp($batch_job);
    @batch_args = split(/\t/, $batch_job, -1);
    $batch_out = shift(@batch_args);
    $batch_pid = fork();
    if (!defined($batch_pid)) {
        push(@batch_status, -1);
        next;
    }
    if ($batch_pid == 0) {
        # the job: the arguments and the standard output of 'perl pal2nal.pl <arguments> > <output file>'
        open(STDOUT, ">", $batch_out) || exit(126);
        open(STDERR, ">", "/dev/null");
        @ARGV = @batch_args;
        $0 = $batch_script;
        undef(@batch_jobs);
        undef(@batch_status);
        $? = 0;
        $! = 0;
        $@ = "";
        &$batch_main();
        exit(0);
    }
    waitpid($batch_pid, 0);
    push(@batch_status, $?);
}

open(BATCH_STATUS, ">", $batch_statusfile) || die "Can't open $batch_statusfile: $!\n";
print BATCH_STATUS join("\n", @batch_status), "\n";
close(BATCH_STATUS);
exit 0;
