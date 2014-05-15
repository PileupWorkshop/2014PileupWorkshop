#!/usr/bin/perl -w
#
# take a file from STDIN and place side-by-side entries corresponding
# to different indices.
#
# usage:
#
#  mergeidx.pl "indexstring1" "indexstring2" [...] < file
#
# OR
#
#  mergeidx.pl "indexstring1" "indexstring2" [...] -f file1 -f file2 [...]
#
# the indexstrings must appear in the same order in the command line
# and the file. In the resulting file the inner loop will be over indices
# the outer loop over files.
#
# Other options:
#
#  -nblank N
#   usually 2 blank lines are required to separate different groups
#   This allows you to modify that requirement.    
#
#  -td 
#   works with topdrawer format
#----------------------------------------------------------------------
#
# Think about how one might process this when there is the option of 
# multiple files. Option struct might be 
# 
#      mergeidx.pl index1 index2 -f file1 -f file2 index3
#

@lines = ();

# sort out the files from the indices among the command-line
# arguments
@indices=();
@filehandles=();
$nfile=0;
$requirednblank=2;
$tdformat=0;
while ($arg = shift @ARGV) {
  if ($arg eq "-f") {
    $file = shift @ARGV;
    if ($file =~ /\.gz$/) {
      open($filehandles[$nfile], "gunzip -c $file|") || die "Could not open $file\n";
    } else {
      open($filehandles[$nfile], "<$file") || die "Could not open $file\n";
    }
    $nfile++;
  } elsif ($arg eq "-nblank") {
    $requirednblank = shift @ARGV;
  } elsif ($arg eq "-td") {
    $tdformat = 1
  } else {
    push @indices, $arg;
  }
}
# if no files were specified, then use STDIN
if ($#filehandles <0) {push @filehandles, STDIN;}

$nrun=-1;

foreach $handle (@filehandles) {
  @localIndices = @indices;

  $active=0;
  $nblank=0;
  
  $indexstring = shift @localIndices;
  
  while ($line = <$handle>) {
    if ($indexstring =~ /^[0-9]+$/) {$indexstring = "index $indexstring";}
    if (!$active && $line =~ /$indexstring/) {
      $nrun++;
      $active = 1;
      $iline = 0;
    }
  
    if ($active) {
      chomp($line);
      if ($nrun == 0) {
        if ($line =~ /^\s*[a-z(]/i) {$lines[$iline] = "# ".$line;}
        else                        {$lines[$iline] = $line;}
      } else {
        $lines[$iline] .= " ".$line;
      }
      $iline++;
      # find end of index session
      $deactivate = 0;
      if ($line =~ /^\s*$/) {
        $nblank++;
        if ($nblank == $requirednblank) {
          $deactivate = 1;
        }
      } elsif ($line =~ /^\s*PLOT\s*$/i) {
        $deactivate = 1;
      } elsif ($nblank < $requirednblank) {
        $nblank=0;
      }

      if ($deactivate) {
        $active = 0;
        $nblank = 0;
        if (!($indexstring = shift @localIndices)) {last;}
      }

    }
  }      # loop over indices
}        # loop over filehandles

# now output
$out = join("\n",@lines);
print $out."\n";
