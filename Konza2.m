LT=[9814006.920 297814.320 9813260.53 295620.94 9816258.91 298521.27 9813443.795 297348.473];
L=LT'

A=[9814006.605 297813.7923 1 0
    297813.7923 -9814006.605 0 1
    9813260.212	295620.3984 1 0
295620.3984 -9813260.212	0 1
9816258.587	298520.7603 1 0
298520.7603 -9816258.587	 0 1
9813443.508	297347.9102 1 0
297347.9102 -9813443.508 0 1];

N=A'*A;
Qxx=inv(N);
d=A'*L;
X=Qxx*d


