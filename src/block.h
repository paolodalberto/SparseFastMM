



#define BM_ = 8 
#define BN_ = 8 

#define EMPTY_BLOCK = 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0  \
    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0		       \
    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0		       \
    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0		       \
    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0		       \
    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0		       \
    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0		       \
    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0		       \


static add_b(COOB C,  COOB A , COOB B){
  
  c = C.data; a = A.data; b = B.data;

  for (int i=0; i<BM_; i++) { 
    c[i][0] = add(a[i][0],b[i][0]);
    c[i][1] = add(a[i][1],b[i][1]);
    c[i][2] = add(a[i][2],b[i][2]);
    c[i][3] = add(a[i][3],b[i][3]);
    c[i][4] = add(a[i][4],b[i][4]);
    c[i][5] = add(a[i][5],b[i][5]);
    c[i][6] = add(a[i][6],b[i][6]);
    c[i][7] = add(a[i][7],b[i][7]);

  

}
static mul_b(COOB C,  COOB A , COOB B){
  
  Mat *c = C.data; *a = A.data; *b = B.data;
  Mat a0,a1,a2,a3,a4,a5,a6,a7;
  Mat t0,t1,t2,t3,t4,t5,t6,t7; 
  for (int i=0; i<BM_; i++){
    a0=a[i][0];a1=a[i][1];a2=a[i][2];a3=a[i][3];a4=a[i][4];a5=a[i][5];a6=a[i][6];a7=a[i][7]; 
    for (int j=0; j<BN_; j++) {
      t0=mul(a0,b[0][j]);
      t1=mul(a1,b[1][j]);
      t2=mul(a2,b[2][j]);
      t3=mul(a3,b[3][j]);
      t4=mul(a4,b[4][j]);
      t5=mul(a5,b[5][j]);
      t6=mul(a6,b[6][j]);
      t7=mul(a7,b[7][j]); 
      c[i][j] = t0+ t1 +t2 +t3+t4+t5+t6+t7;
    }

}
