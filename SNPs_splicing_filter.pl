#! /usr/bin/perl 
#ZZ,YZ,YS,OZ,MS,DL,DB,CB
						
###################################################################
%hash;  
%hashbed;
%hashchr;
$chuan="_";
open (T0,"@ARGV[1]");
        while($vcf=<T0>){
        chomp $vcf;
	@VCF =split(/\s+/ ,$vcf);
	if($VCF[0]=~ m/^#/){
	next;
	}else{
	$VCF[0]=~ s/^Chr//;
	$pos=$VCF[0].$chuan.$VCF[1];
	$hash{$pos}="";	
	}
	}

open (T2,"@ARGV[2]");
        while($bed=<T2>){
        chomp $bed;
        @BED =split(/\s+/ ,$bed);
        $hashbed{@BED[1]}=@BED[2];
	$hashchr{@BED[1]}=@BED[0];
        }
%hashline;
open (T1,"@ARGV[0]");
	while($line=<T1>){
	chomp $line;
	@A =split (/\s+|[;\-:]/ ,$line);
	for (5 .. $#A){
	$ZZ=$A[1]+$A[$_];
	$Epos=$A[0].$chuan.$ZZ;
	if(exists($hash{$Epos})){
	next;
	}
	#open (R,">>filter_SNPs.txt");
	$infoline=$A[0]."\t".$ZZ."\t".$A[3]."\t".$A[4]."\t".$A[0].":".$A[1]."-".$A[2]."\n";
	push (@info, $infoline);
	
for $k (keys(%hashchr)){
	if($hashchr{$k} eq $A[0] && $ZZ>=$k && $ZZ<=$hashbed{$k}){
	push (@info2,$infoline);
	}
	}
}
}
# foreach (@info){print "$_";}
%hash2;
foreach $aa (@info2){
$hash2{$aa}=1;
}
for(my $i=0;$i<@info;$i++)
{
delete $info[$i] if($hash2{$info[$i]}==1);
}
print "@info\n";
#foreach (@uniq){print "$_";}
# print @info,"\n";
#for $row in (@info){
#print "$row\n";}
#open (T3,"filter_SNPs.txt");
#        while($info=<T3>){
#        chomp $info;
#        @B =split (/\s+|[;\-:]/ ,$info);
#open (T4,"@ARGV[2]");
#        while($bed=<T4>){
#   	  chomp $bed;
#        @BED =split(/\s+/ ,$bed);
#        if($BED[0] eq $B[0]){ 
#	 $Q1=$BED[1]-$B[1];
#	$Q2=$B[1]-$BED[2];
#	if($Q1>0 || $Q2>0){
#	print "$info\n";
#	}
#	}
#	}
#}
#    	#print "$A[0]\t$ZZ\t$A[3]\t$A[4]\t";
        #print "$A[0]:$A[1]\-$A[2]";
        #print "\n";}
 	#print "$A[0]\t$ZZ\t$A[3]\t$A[4]\t";
        #print "$A[0]:$A[1]\-$A[2]";
        #print "\n";
	
