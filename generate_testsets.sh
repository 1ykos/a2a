#!/bin/bash
#awk 'NF==4{print $1,$2,$3,$4;next}{if ($1<10000) print $0}' premerge > first10000
awk '
NF==4{
  for (i=0;i<=7.5;++i) {
    file=ARGV[1]"_"i;
    print $1,$2,$3,$4 > file;
  }
  e=sqrt(3)*(rand()+rand()+rand()+rand()-2); # standard normal variable
  e*=sqrt(0.5); # so difference between two is 1
  next
}
{
  for (i=0;i<=7.5;++i) {
    file=ARGV[1]"_"i;
    if (int($1)%2) {
      print $1,$2,$3,$4,$5*(1-e*0.5**i),$6*(1-e*0.5**i)**2,$7,$8,$9,$10,$11,$12 > file;
    } else {
      print $1,$2,$3,$4,$5*(1+e*0.5**i),$6*(1+e*0.5**i)**2,$7,$8,$9,$10,$11,$12 > file;
    }
  }
}' first10000
