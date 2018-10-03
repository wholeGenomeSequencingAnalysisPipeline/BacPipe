use strict;
use Bio::TreeIO;

my $in = new Bio::TreeIO(-file => $ARGV[0],
-format => 'newick');
my $out = new Bio::TreeIO(-file => '>'.$ARGV[0].'.txt',
-format => 'tabtree');
while( my $tree = $in->next_tree ) {
 $out->write_tree($tree);
}
my $in = new Bio::TreeIO(-file => $ARGV[0],
-format => 'newick');
my $out = new Bio::TreeIO(-file => '>'.$ARGV[0].'.svg',
-format => 'svggraph');
while( my $tree = $in->next_tree ) {
	$out->write_tree($tree);
}

#system "qlmanage -t -s 1000 -o . mytree.svg";
#system "sips -s format gif mytree.svg.png --out mytree.gif";
