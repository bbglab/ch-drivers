{
	pos=$1":"$2"-"$2;
	mch=0;
	cmd="tabix " tabixfile " " pos;
    while ( ( cmd |getline t ) > 0 ) {
    	split(t,r,"\t");
    	if ($3 == r[3] && $4 == r[4]) {
    		mch=1;
    		break;
		}
	};
	close(cmd);
	if (mch == 0) {
		print $0
	}
}