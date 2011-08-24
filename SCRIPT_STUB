#!/usr/bin/perl

use strict;
use warnings;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------

my $usage = "

Synopsis:

  script_stub script_name

Description:

  This script will create a new perl script stub file named script_name.

";

my $file = shift;

die $usage unless $file;

upgrade($file) if -e $file;

open (my $OUT, '>', $file) or die "Can't open OUT: $!\n";

print $OUT <<"END";
#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my \$usage = "

Synopsis:

$file <How the hell do you use this thing>

Description:

No really, how the hell do you use this thing!

";


my (\$help);
my \$opt_success = GetOptions('help'    => \\\$help,
			      );

die \$usage if \$help || ! \$opt_success;

my \$file = shift;
die \$usage unless \$file;
open (my \$IN, '<', \$file) or die "Can't open \$file for reading\\n\$!\\n";

while (<\$IN>) {

}

#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------

sub {

}

END
    1;
finish();

sub upgrade {
	my $file = shift;
	open (my $IN, '<', $file) or die "Can't open $file for reading: $!\n";

	my $main_text = "
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my \$usage = \"

Synopsis:

$file <How the hell do you use this thing>

Description:

No really, how the hell do you use this thing!

\";
";

	my $subs_text = "
#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------

sub {

}
";

	my @lines = (<$IN>);
	my ($main, $subs, $count);
	for my $line (@lines) {
		$count++;
		$main = $count if $line =~ /^use+/;
		$subs ||= $count if $line =~ /^sub\s+/;
		exec ('emacs', $file) or die "Couldn't run 'emacs $file':\n$!\n"
		  if $line =~ /^\#----------------------------------- MAIN -------/;
	}
	$subs ||= $count;
	close $IN;
	exec ('emacs', $file) or die "Couldn't run 'emacs $file':\n$!\n" unless defined $main;
	my $backup = $file . '.bak';
	`mv $file $backup`;

	open (my $OUT, '>', $file) or die "Can't open $file for writing: $!\n";
	$count = 0;
	for my $line (@lines) {
		$count++;
		print $OUT $main_text if $count == $main + 1;
		print $OUT $line;
		print $OUT $subs_text if $count == $subs;
	}
	close $OUT;
	finish();
}

sub finish {

	`chmod +x $file`;
	exec ('emacs', $file) or die "Couldn't run 'emacs $file':\n$!\n";
}
