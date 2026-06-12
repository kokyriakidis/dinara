#!/usr/bin/perl
# Gate 1 diagnostic: measure collapse/tangle in a Shasta2AnchorGraph GFA.
#
# Reads a GFA written by Shasta2AnchorGraph::writeGfa:
#   S <vid> * LN:i:1 wn:i:<window> ws:Z:<fw|rc>
#   L <src> + <dst> + 0M RC:i:<cov> [sp:i: sn:i: sr:i:]
#
# Reports baseline metrics for the de-novo detangling work:
#   - vertices/edges, distinct normalized windows
#   - vertices-per-window distribution (collapse proxy)
#   - intra- vs inter-window edges
#   - per-window predecessor/successor counts (tangle windows)
#   - strand-strand edges (fw<->rc of same window: hairpin/collapse signature)
#   - hairpin windows (connect to their own normalized window)
#   - span-product (sp*sn) distribution + log2 histogram (drives T threshold)
#   - shared-read (sr) distribution (drives N threshold)
#
# Usage: perl Gate1AnchorGraphCollapseStats.pl Shasta2AnchorGraph.gfa [edges.csv]
#   If edges.csv is given, a per-edge dump of sp,sn,product,sr is written there
#   so T can be re-derived later without re-running.

use strict;
use warnings;

my $file = $ARGV[0] or die "Usage: $0 <anchorgraph.gfa>\n";
open(my $fh, '<', $file) or die "Cannot open $file: $!\n";

my %vWindow;   # vid -> normalized window
my %vStrand;   # vid -> fw|rc
my %winVerts;  # window -> count of vertices (any strand)
my $sLines = 0;
my $lLines = 0;

# Edge bookkeeping at normalized-window granularity.
my %succ;      # win -> { succWin => count }  (different window only)
my %pred;      # win -> { predWin => count }
my $intra = 0; # edges within same (window,strand)
my $interSameWinDiffStrand = 0; # fw<->rc of same normalized window
my $interDiffWin = 0;

my @edges; # [srcWin, srcStrand, dstWin, dstStrand]

# Inter-window support evidence (drives N and T thresholds).
my @spanProducts;  # sp*sn per inter-window edge that carries sp/sn
my @sharedReads;   # sr per inter-window edge that carries sr
my $csvOut;        # optional per-edge CSV dump
my $csvName = $ARGV[1];
if (defined $csvName) {
    open($csvOut, '>', $csvName) or die "Cannot open $csvName for writing: $!\n";
    print $csvOut "srcWin,srcStrand,dstWin,dstStrand,sp,sn,product,sr\n";
}

while (my $line = <$fh>) {
    if ($line =~ /^S\t(\S+)\t/) {
        $sLines++;
        my $vid = $1;
        my ($wn) = $line =~ /\bwn:i:(\d+)/;
        my ($ws) = $line =~ /\bws:Z:(fw|rc)/;
        if (defined $wn) {
            $vWindow{$vid} = $wn;
            $vStrand{$vid} = defined($ws) ? $ws : '?';
            $winVerts{$wn}++;
        }
    } elsif ($line =~ /^L\t(\S+)\t\S\t(\S+)\t/) {
        $lLines++;
        my ($src, $dst) = ($1, $2);
        next unless exists $vWindow{$src} && exists $vWindow{$dst};
        my $sw = $vWindow{$src}; my $ss = $vStrand{$src};
        my $dw = $vWindow{$dst}; my $ds = $vStrand{$dst};
        push @edges, [$sw, $ss, $dw, $ds];

        # Support tags (present only on inter-window edges).
        my ($sp) = $line =~ /\bsp:i:(\d+)/;
        my ($sn) = $line =~ /\bsn:i:(\d+)/;
        my ($sr) = $line =~ /\bsr:i:(\d+)/;

        if ($sw == $dw) {
            if ($ss eq $ds) { $intra++; }
            else            { $interSameWinDiffStrand++; }
        } else {
            $interDiffWin++;
            $succ{$sw}{$dw}++;
            $pred{$dw}{$sw}++;

            if (defined $sp && defined $sn) {
                my $prod = $sp * $sn;
                push @spanProducts, $prod;
                push @sharedReads, $sr if defined $sr;
                if ($csvOut) {
                    printf $csvOut "%s,%s,%s,%s,%d,%d,%d,%s\n",
                        $sw, $ss, $dw, $ds, $sp, $sn, $prod,
                        (defined $sr ? $sr : '');
                }
            }
        }
    }
}
close($fh);
close($csvOut) if $csvOut;

# Distinct windows.
my %allWins;
$allWins{$_} = 1 for keys %winVerts;
$allWins{$_} = 1 for keys %succ;
$allWins{$_} = 1 for keys %pred;
my $nWins = scalar keys %allWins;

# Tangle windows: >1 distinct predecessor OR >1 distinct successor.
my $tangleWins = 0;
my $multiPred = 0;
my $multiSucc = 0;
my $multiBoth = 0;
for my $w (keys %allWins) {
    my $np = exists $pred{$w} ? scalar keys %{$pred{$w}} : 0;
    my $ns = exists $succ{$w} ? scalar keys %{$succ{$w}} : 0;
    $multiPred++ if $np > 1;
    $multiSucc++ if $ns > 1;
    $multiBoth++ if $np > 1 && $ns > 1;
    $tangleWins++ if $np > 1 || $ns > 1;
}

# Hairpin windows: window connects to its own normalized window
# (via the fw<->rc edges we counted, or any succ/pred == self).
my %hairpinWin;
for my $e (@edges) {
    my ($sw, $ss, $dw, $ds) = @$e;
    if ($sw == $dw && $ss ne $ds) {
        $hairpinWin{$sw} = 1;
    }
}
my $nHairpin = scalar keys %hairpinWin;

# Vertices-per-window distribution (collapse proxy).
my @counts = sort { $a <=> $b } values %winVerts;
my $nW = scalar @counts;
my ($vmin, $vmax, $vsum) = (0, 0, 0);
if ($nW) {
    $vmin = $counts[0];
    $vmax = $counts[-1];
    $vsum += $_ for @counts;
}
my $vmean = $nW ? $vsum / $nW : 0;
my $median = $nW ? $counts[int($nW/2)] : 0;
my $p99 = $nW ? $counts[int($nW*0.99)] : 0;

# Collapse proxy: windows whose vertex count is an outlier (> 3x median).
my $collapseThresh = $median > 0 ? 3 * $median : ($vmean > 0 ? 3 * $vmean : 1e9);
my $collapsedWins = 0;
$collapsedWins++ for grep { $_ > $collapseThresh } @counts;

printf "=== Gate 1: Anchor-graph collapse/tangle stats ===\n";
printf "File: %s\n\n", $file;

printf "S lines (vertices):            %d\n", $sLines;
printf "L lines (edges):               %d\n", $lLines;
printf "Distinct normalized windows:   %d\n\n", $nWins;

printf "-- Edge breakdown --\n";
printf "Intra-window (same win+strand):        %d\n", $intra;
printf "Strand-strand (same win, fw<->rc):     %d   <- hairpin/collapse signature\n", $interSameWinDiffStrand;
printf "Inter-window (different window):        %d\n\n", $interDiffWin;

printf "-- Window vertex counts (collapse proxy) --\n";
printf "min / median / mean / p99 / max:  %d / %d / %.1f / %d / %d\n", $vmin, $median, $vmean, $p99, $vmax;
printf "Windows > 3x median verts:        %d   <- candidate collapsed windows\n\n", $collapsedWins;

printf "-- Tangle windows --\n";
printf ">1 predecessor:                   %d\n", $multiPred;
printf ">1 successor:                     %d\n", $multiSucc;
printf ">1 pred AND >1 succ:              %d\n", $multiBoth;
printf "Tangle (>1 pred OR >1 succ):      %d  (%.1f%% of windows)\n",
    $tangleWins, ($nWins ? 100.0*$tangleWins/$nWins : 0);
printf "Hairpin windows (self fw<->rc):   %d   <- strand-strand contacts\n\n", $nHairpin;

# ---- Support-evidence distributions (drive N and T thresholds) ----

sub pct {
    my ($aref, $q) = @_;
    my $n = scalar @$aref;
    return 0 unless $n;
    my $i = int($q * ($n - 1) + 0.5);
    $i = 0 if $i < 0; $i = $n - 1 if $i >= $n;
    return $aref->[$i];
}

sub log2bucket {
    my ($v) = @_;
    return '0' if $v <= 0;
    my $b = 0;
    my $lo = 1;
    while ($lo * 2 <= $v) { $lo *= 2; $b++; }
    my $hi = $lo * 2;
    return sprintf("%d-%d", $lo, $hi - 1);  # e.g. 1-1, 2-3, 4-7, 8-15
}

printf "-- Span-product (sp*sn) distribution on inter-window edges --\n";
printf "    (T = span-product threshold is read from here)\n";
my $nSP = scalar @spanProducts;
if ($nSP) {
    my @sp = sort { $a <=> $b } @spanProducts;
    my $sum = 0; $sum += $_ for @sp;
    printf "edges with sp/sn:                 %d\n", $nSP;
    printf "min / p10 / p25 / median:         %d / %d / %d / %d\n",
        $sp[0], pct(\@sp, 0.10), pct(\@sp, 0.25), pct(\@sp, 0.50);
    printf "p75 / p90 / p95 / p99 / max:      %d / %d / %d / %d / %d\n",
        pct(\@sp, 0.75), pct(\@sp, 0.90), pct(\@sp, 0.95), pct(\@sp, 0.99), $sp[-1];
    printf "mean:                             %.1f\n", $sum / $nSP;

    # Log-scale histogram (power-of-2 buckets), ordered by bucket low edge.
    my %hist;
    $hist{ log2bucket($_) }++ for @sp;
    my %lowOf;
    for my $b (keys %hist) {
        my ($lo) = $b =~ /^(\d+)/;
        $lowOf{$b} = $lo;
    }
    printf "  span-product histogram (log2 buckets):\n";
    for my $b (sort { $lowOf{$a} <=> $lowOf{$b} } keys %hist) {
        my $c = $hist{$b};
        my $bar = '#' x int(60 * $c / $nSP + 0.5);
        printf "    %-12s %8d  %s\n", $b, $c, $bar;
    }
} else {
    printf "  (no inter-window edges carried sp/sn tags)\n";
}
printf "\n";

printf "-- Shared-read count (sr) distribution on inter-window edges --\n";
printf "    (N = minimum corroborating-read count is read from here)\n";
my $nSR = scalar @sharedReads;
if ($nSR) {
    my @sr = sort { $a <=> $b } @sharedReads;
    my %srHist;
    for my $v (@sr) {
        my $k = $v >= 4 ? '4+' : "$v";
        $srHist{$k}++;
    }
    for my $k ('1', '2', '3', '4+') {
        next unless exists $srHist{$k};
        my $c = $srHist{$k};
        my $bar = '#' x int(60 * $c / $nSR + 0.5);
        printf "    sr=%-3s %8d  (%.1f%%)  %s\n",
            $k, $c, 100.0 * $c / $nSR, $bar;
    }
    printf "  sr=1 edges are single-read joins (chimera-prone, count<2).\n";
} else {
    printf "  (no inter-window edges carried sr tags)\n";
}

if (defined $csvName) {
    printf "\nPer-edge CSV written to: %s\n", $csvName;
}
