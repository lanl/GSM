#!/usr/bin/perl
#
# fxref [-xref|-tree] [-file] [-nofunc] [-ignore list]  [-root name] [files...]
#
#   fbb, 6/24/01
#
########################################################################
$FPP   = "$ENV{HOME}/bin/fpp";
$DEBUG = 1;

$MAXDEPTH = 25;          ### limit for depth of call-tree
$SUBPROGS = "all";       ### flag, scan for both subr's & fcn's (all)
$LEN      = 3;           ### length for routine name indents

$DIR  = `pwd`;     chop( $DIR  );
$DATE = localtime;
$PRT  = "stdout";

$TREE = "yes";
$XREF = "yes";
$PROG = "";
$DEP  = "";
undef @NO_DESCEND;

while( $ARGV[0] =~ /^-/ ) {
  $_ = shift;
  if( /^-xref/ ) { $TREE="no";   $XREF="yes";  next; }
  if( /^-tree/ ) { $TREE="yes";  $XREF="no";   next; }
  if( /^-root/ ) { $PROG=shift;  next; }
  if( /^-nofunc/ ) { $SUBPROGS="no-functions";  next; }
  if( /^-file/ ) { $PRT = "no";  next; }
  if( /^-ignore/ ) { $d=shift; push(@NO_DESCEND,$d); next;}
  if( /^-D/ ) { $DEP = "$DEP ".$_;      next; }
  if( /^-U/ ) { $DEP = "$DEP ".$_;      next; }
  if( /^--/ ) { next; }
  die( "USAGE:  fxref [-xref|-tree] [-file] [-ignore list]  [-root name] [--] [files...]\n" );
}
@flist = <@ARGV>;
if( $flist[0] eq '' ) { @flist = <*.F90 *.f90 *.F *.f *.F95 *.f95>; }
$IN = "";
#-----------------------------------------------------------------------
###
### scan files once, to find 
###     name, type, first-line, last-line, calls-to-subs,
###        include files, use files
### for each prog, func, subr, entr, module
###
print STDERR ".....scanning files\n";
while( $FILE = shift(@flist) ) {
    $INSIDE = -1;
    $INTERFACE_SKIP = 0;
    if( $FILE =~ /\.F90$/ ) {
      open( F, "$FPP $DEP -- $FILE |") || die("fxref: can't open file $FILE\n");
    }
    else {
      open( F, "<$FILE" ) || die("fxref: can't open file $FILE\n");
    }
    while( <F> ) {
      chop;
      #
      ### remove comments (simpleminded fix)
      s/!.*//;
      ### remove strings
      if( ! /^\s*#\s*include/ ) { s/"[^"]*"//g;  s/'[^']*'//g; }
      ### skip some types of statement
      if( /^[cC*]/ ) { next; }
      if( /^\s*$/  ) { next; }
      if( /^\s*end\s*(if|do|select)/i ) { next; }
      if( /module\s+procedure/i ) { next; }

      if(    /^\s*interface/i      ) { $INTERFACE_SKIP = 1; next; }
      elsif( /^\s*end\s*interface/ ) { $INTERFACE_SKIP = 0; next; }
      if( $INTERFACE_SKIP ) { next; }
        
      #
      if( /^\s*(SUBROUTINE)\s+(\w+)/i   ||
          /^\s*(FUNCTION)\s+(\w+)/i     ||
          /^\s*(PROGRAM)\s+(\w+)/i      ||
          /^\s*(MODULE)\s+(\w+)/i       ||
          /^\s*(ENTRY)\s+(\w+)/i          ) {
        ($n=$2)             =~ tr/A-Z/a-z/;
        ($t=substr($1,0,4)) =~ tr/A-Z/a-z/;
        push(@name,$n);
        push(@type,$t);
        push(@line,$.);
        push(@last, 0);
        push(@file,$FILE);
        $seqn{$n} = $#name;
        $cto = "";
        if( $t eq 'entr' ) {
          if( $INSIDE != -1 ) { $cto = "$name[$INSIDE],"; }
        }
        else {
          if( $INSIDE != -1 ) { $last[$INSIDE] = $.-1; }
          $INSIDE = $#name;
        }
        push(@callto, $cto);
        push(@callby, "");
      }
      elsif( /^\s*#\s*include\s+"(.*)"/i  ||
             /^\s*#\s*include\s+<(.*)>/i  ||
	     /^\s*include\s+'(.*)'/i      ||
	     /^\s*use\s+(\w+)/i    ) {
        if( $INSIDE != -1 ) { $include[$INSIDE] .= "$1,"; }
      }
      elsif( /\Wcall\s+(\w+)/i ) {
        if( $INSIDE != -1 ) {
	  ($c=$1) =~ tr/A-Z/a-z/;
          if( index(",$callto[$INSIDE]", ",$c,") < 0 ) {
	      $callto[$INSIDE] .= "$c,";
          }
        }
      }
      elsif( /^\s*(end|contains)\s*$/i ) {
        if( $INSIDE != -1 ) {
          $last[$INSIDE] = $.;
          $INSIDE = -1;
        }
      }
    }
    close(F);
}
print STDERR "\n";
if( $#name < 0 ) { die("***** no subprogs -- is this Fortran?"); }
#-----------------------------------------------------------------------
###
### scan files again, to add functions to calls-to list
###
if( $SUBPROGS eq "all" ) {
    undef @funcs;
    for( $i=0; $i<@name; $i++ ) {
	if( $type[$i] eq "func" ) {
	    push( @funcs, $name[$i]);
	}
	elsif( $type[$i] eq "entr" ) {
	    chop( $n=$callto[$i] );
	    $j = $seqn{$n};
	    if( $type[$j] eq "func" ) {
		push( @funcs, $n);
	    }
	}
    }
    if( @funcs>0 ) {
	print STDERR ".....re-scanning files for function calls\n";
	for( $i=0; $i<@name; $i++ ) {
	    $filei = $file[$i];
	    $linei = $line[$i];
	    $lasti = $last[$i];
	    $ctoi  = $callto[$i];
	    print STDERR ".";
            if( $filei =~ /\.F90$/ ) {
              open( F, "$FPP $DEP -- $filei |") || die("can't open file $filei\n");
            }
            else {
              open( F, "<$filei" ) || die("can't open file $filei\n");
            }
            $INTERFACE_SKIP = 0;
	    while( <F> ) {
		if( $. >= $lasti ) { last; }
		if( $. <= $linei ) { next; }
		### remove comments (simpleminded fix)
		s/!.*//;
		### remove strings
		if( ! /^\s*#\s*include/ ) { s/"[^"]*"//g; s/'[^']*'//g;}
	        ### skip some types of statement
		if( /^[cC*]/        ) { next; }
		if( /^\s*#/         ) { next; }
		if( /^\s*external/i ) { next; }
		if( /^\s*entry/i    ) { next; }
		if( /^\s*$/         ) { next; }
		if( /^\s*end\s*(if|do|select|interface)/i ) { next; }

                if(    /^\s*interface/i      ) { $INTERFACE_SKIP = 1; next; }
                elsif( /^\s*end\s*interface/ ) { $INTERFACE_SKIP = 0; next; }
                if( $INTERFACE_SKIP ) { next; }
		foreach $namej (@funcs) {
		    if( /\W$namej\s*\([^=]*/i ) {
			if( index(",$ctoi", ",$namej,") <0 ) {
			    $ctoi .= "$namej,";
			}
		    }
		}
	    }
	    $callto[$i] = $ctoi;
	    close(F);
	}
	print STDERR "\n";
    }
}
else {
    print STDERR "...SKIPPING re-scan for function calls\n";
}
#-----------------------------------------------------------------------
###
### set up the call-by lists
###
for( $i=0; $i<@name; $i++ ) {
    chop( $ctoi = $callto[$i] );
    undef @to;
    @to = split(/,/, $ctoi);
    foreach $n (@to) {
        $j = $seqn{$n};
        if( $j ne "" ) {
            if( index(",$callby[$j]", "$name[$i],") < 0 ) {
                $callby[$j] .= "$name[$i],";
            }
        }
    }
}
#-----------------------------------------------------------------------
###
### strip trailing comma from callto, callby, & include lists
###
for( $i=0; $i<@name; $i++ ) {
	chop( $callto[$i] );
	chop( $callby[$i] );
	chop( $include[$i] );
}
#-----------------------------------------------------------------------
###
### find first program name
###
    if( $PROG eq "" ) {
	$PROG = "?";
	for( $i=0; $i<@name; $i++ ) {
	    if( $type[$i] eq "prog" ) { $PROG=$name[$i]; last; }
	}
    }
#-----------------------------------------------------------------------
###
### print info, with call-to, call-by, & include files
###
sub print_xref_info {
  my( $label, $comma_lst ) = @_;
  my( @lst, $L, $j, $n, $t );
  printf P "%s%40s%-20s", $IN, " ", $label;
  @lst = reverse( sort( split(/,/, $comma_lst) ) );
  $L = length($IN)+60;
  for( $j=0, $n=$#lst; $j<=$n; $j++ ) {
    $t = pop(@lst);
    if( $L+length($t)>120 || $L+11>120 ) {
	printf P "\n%s%60s", $IN, " ";
	$L = length($IN)+60;
    }
    printf P "%-10s ", $t;
    $L = $L + (length($t)>10?length($t):10)+1;
  }
  printf P "\n";
}

sub print_xref {
  my( $n ) = $_[0];
  my( $i );

  $i = $seqn{$n};
  printf P "%s%-4s %-16s %s(%d,%d)\n",
    $IN, $type[$i],$name[$i],$file[$i],$line[$i],$last[$i];

  &print_xref_info( "Uses modules:  ", $include[$i] );
  &print_xref_info( "Makes calls to:", $callto[$i]  );
  &print_xref_info( "Is called from:", $callby[$i]  );
}

if( $XREF eq "yes" ) {
    my( $n );
    if( $PRT eq "stdout" ) {
        open( P, ">&STDOUT");
    }
    else {
	print STDERR "\n.....xref output to file:  ${PROG}_xref.txt\n\n";
	open( P, ">${PROG}_xref.txt") || die("Can't open ${PROG}_xref.txt");
    }
    foreach $n (sort(keys(%seqn))) {
      &print_xref( $n );
    }
    close(P);
}
#-----------------------------------------------------------------------
###
### print calling tree
###
sub print_tree {
    my( $root, $indent, $depth, $chain ) = @_;
    my( $i, $n, @to, $ldr, $str );


    $depth = $depth + 1;
    $n     = $seqn{$root};
    if( $n ) {
      $str = sprintf("%-s", $indent.$root);
      printf P "%-s\n", $str .  " " x (80-length($str)) . $file[$n];
      #printf P "%-s\n", $indent.$root."\t\t\t".$file[$n];
    }
    else {
      printf P "%-s\n", $indent.$root;
    }

    $indent .= "." . " " x $LEN;
        
    foreach $i (@NO_DESCEND) {
        if( $root =~ /^$i/i ) {
	      printf P $indent."---etc---\n";
            return;
        }
    }
    if( $root =~ /^\.\.\./ ) { 
#	printf P "\n";
	return;
    }
    elsif( index(",$chain", ",$root") >= 0 ) {
    	printf P $indent."RECURSIVE\n";
    }
    elsif( $depth > $MAXDEPTH ) {
	$ldr = "\n**********";
	printf P "$ldr call-tree-depth > $MAXDEPTH";
	printf P "$ldr Recursive ??? program bug ??? fxref bug ???\n\n";
	$callto[$n] = "...??";
    }
    elsif( $n ne ""  &&  $callto[$n] ne "" ) {
	#@to = sort( split(/,/, $callto[$n])  );
	@to = split( /,/,  $callto[$n] );
	for( $i=0; $i<@to; $i++ ) {
	    #if( $i > 0 ) { printf P "$indent"; }
	    &print_tree( $to[$i], $indent, $depth, $chain.",$root" );
	}
	$callto[$n] = "...etc...";
    }
   #else {
   #    printf P "\n";
   #}
}

if( $TREE eq "yes" ) {
    my( $i, $n, $t, $c );
    if( $PRT eq "stdout" ) {
	open( P, ">&STDOUT");
    }
    else {
	print STDERR "\n.....call-tree output to file:  ${PROG}_tree.txt\n\n";
	open( P, ">${PROG}_tree.txt") || die("Can't open ${PROG}_tree.txt");
    }
	
    for( $i=0; $i<@name; $i++ ) {
	$n = $name[$i];
	$t = $type[$i];
	$c = $callby[$i];
	if( $n eq "$PROG" ) {
	    &print_tree( $n, "$IN", 0, "" );
	}
    }
    close(P);
}
#-----------------------------------------------------------------------
exit;
