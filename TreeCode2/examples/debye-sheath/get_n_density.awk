BEGIN{c = 0;}
{
	if($1 ~ /^=*$/) {c++}
	else{
		if(c == n){print $0}
	}
	if(c > n){exit}
}
