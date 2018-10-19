#include "poly_grid.h"
#include "chebyshev2.h"
#include "real_def.h"
#if DDDD
#include "dd_utils.h"
#include "lapack_dd.h"
#include "str_dd.h"
#elif MPFR
#include "lapack_mpfr.h"
#include "str_mpfr .h"
#else
#include "lapack.h"
#include "str_double.h"
#endif

int main(int argc, char *argv[])
{
    if (argc < 5)
        throw gen_err("\nusage:  ./test_jacobi n, xb, aa, bb, t1, t2, t3 with -1<=t1<=1, t2=1+t1, t3=1-t1\n");
    int n; 
    Real xb = 0.0;
    Real aa = 0.0;
    Real bb = 0.0;
    Real t1 = 0.0, t2 = 0.0, t3 = 0.0;
    Real temp = 0.0;
    sscanf(argv[1], "%d", &n);
    str_to_real(argv[2], xb);
    str_to_real(argv[3], aa);
    str_to_real(argv[4], bb);
    // str_to_real(argv[5], t1);
    //str_to_real(argv[6], t2);
    //str_to_real(argv[7], t3);
    poly_grid<Real> *pg1 = new jacobi_grid<Real>(n,aa,bb,-xb,xb);
    poly_grid<Real> *pg2 = new jacobi_grid<Real>(n,aa,bb,0.0,xb); //0,1,...,n-1
    poly_grid<Real> *pg3 = new jacobi_grid<Real>(n,aa,bb,xb,0.0);
    mp_mat<Real> J(n,n,0.0);
    mp_mat<Real> Phi(n,n);
    valarray<Real> tb2(0.0,2*n), tb3(0.0,2*n);
    // calculating tb2 and tb3 using method 2 (method1 is not as accurate as 2. Also there will be underflow and overflow problem when calculating tb's. Makes it even harder to use)
    if (0) {
      valarray<Real> d0(0.0,n+1), d0_3(0.0,n+1);
      valarray<Real> d1(0.0,n+1), d1_3(0.0,n+1);
      d0[0] = d0_3[0] = 1.0;
      d0[1] = -pg2->a[0];
      d0_3[1] = -pg3->a[0];
      d1[0] = d1_3[0] = 0.0;
      d1[1] = d1_3[1] = 1.0;
      for (int i=2; i<=n; i++) {
	d0[i] = -pg2->a[i-1]*d0[i-1]-pg2->b[i-1]*d0[i-2];
	d1[i] = d0[i-1]-pg2->a[i-1]*d1[i-1]-pg2->b[i-1]*d1[i-2];
	d0_3[i] = -pg3->a[i-1]*d0_3[i-1]-pg3->b[i-1]*d0_3[i-2];
	d1_3[i] = d0_3[i-1]-pg3->a[i-1]*d1_3[i-1]-pg3->b[i-1]*d1_3[i-2];
      }
      for (int i=0; i<n; i++) {
	tb2[2*i+1] = -d0[i+1]/d0[i];
	tb3[2*i+1] = -d0_3[i+1]/d0_3[i];
      }
      for (int i=1; i<n; i++) {
	tb2[2*i] = (d0[i]-d1[i+1]-tb2[2*i+1]*d1[i])/(d1[i]+tb2[2*i-1]*d1[i-1]);
	tb3[2*i] = (d0_3[i]-d1_3[i+1]-tb3[2*i+1]*d1_3[i])/(d1_3[i]+tb3[2*i-1]*d1_3[i-1]);
      }
    }
    // use recurrence relation to find tb
    if (0) {
      for (int i=0; i<n; i++) {
	tb2[2*i+1] = pg2->a[i]-tb2[2*i];
	//cout << i << " " << pg2->a[i] << endl;
	printf("%d is %s\n", 2*i+1, str(tb2[2*i+1],0) );
	//tb3[2*i+1] = pg3->a[i]-tb3[2*i];
	// cout << 2*i+1 << " " << tb2[2*i+1] << endl;
	if (i<n-1) {
	  tb2[2*i+2] = pg2->b[i+1]/tb2[2*i+1];
	  //cout << i << " " << pg2->b[i+1] << endl;
	  printf("%d is %s\n", 2*i+2, str(tb2[2*i+2],0) );
	  //tb3[2*i+2] = pg3->b[i+1]/tb3[2*i+1];
	  // cout << 2*i+2 << " " << tb2[2*i+2] << endl;
	}
      }
    }
    // use formula of Gegenbauer polynomial directly
    Real tmp;
    if (1) {
      tmp = 2*bb+1;
      for (int i=0; i <n; i++) {
	tb2[2*i+1] = (2*i+1+tmp)*(2*i+1+2*aa+tmp)/(4*(2*i+1+aa+bb)*(2*i+2+aa+bb));
	if (i<n-1) {
	  tb2[2*i+2] = (2*i+2)*(2*i+2+2*aa)/(4*(2*i+2+aa+bb)*(2*i+3+aa+bb));
	}
      }
    }
      
    Real err = 0.0;
    Real x = -0.99;
    // check if continued fraction converges to decide if the recurrence relation posesses a minimal solution. == finding continued fraction A_n/B_n where A's and B's are solutions with initial value 1,0 and 0,1
    if (0) {
      for(int j=1; j < 4000; j++) {
	tmp = 0.0;
	for(int i=j; i>1; i--) {
	  tmp = sqrt(tb2[i-1]/tb2[i])/(-x/sqrt(tb2[i])-tmp);
	}
	cout << str(tmp,0) << endl;
      }
    }

    if (0) {
      for(int j=1; j < 2000; j++) {
	tmp = 0.0;
	for(int i=j; i>1; i--) {
	  tmp = sqrt(pg1->b[i-1]/pg1->b[i])/((pg1->a[i-1]-x)/sqrt(pg1->b[i])-tmp);
	}
	cout << str(tmp,0) << endl;
      }
    }

    // find eigenvalues of 2*n matrix with the first row storing 1,0 .., the second 0, 1, ...
    if (0) {
      tmp = 0.01;
      int info = 0;
      int N = 800;
      mp_mat<Real> mat(2,N, 0.0);
      valarray<Real> eigen(0.0, 2);
      mat(0,0) = 1.0;
      mat(1,0) = 0.0;
      mat(0,1) = 0.0;
      mat(1,1) = 1.0;
      for(int i=2; i<N; i++) {
	mat(0,i) = (tmp*mat(0,i-1)-sqrt(tb2[i-1])*mat(0,i-2))/sqrt(tb2[i]);
	mat(1,i) = (tmp*mat(1,i-1)-sqrt(tb2[i-1])*mat(1,i-2))/sqrt(tb2[i]);
	dgesvd('N', 'N', 2, i+1, mat.p, i+1, &eigen[0], NULL, 1, NULL, 1, &info);
	// cout << str(eigen[0]/eigen[1], 0) << endl;
	cout << i << " " <<  eigen[0] << " " << eigen[1] << endl;
      }
      //mat.dump("haha",17);
      //dgesvd('N', 'N', 2, N, mat.p, N, &eigen[0], NULL, 1, NULL, 1, &info);
      //cout << N << " " <<  eigen[0] << " " << eigen[1] << endl;
    }

	
    // 3 xgrids
    valarray<Real> x1(0.0,n), x2(0.0,2*n); //x3(0.0,2*n);
    valarray<Real> E1(0.0,n), E2(0.0,2*n); //E3(0.0,2*n);
    mp_mat<Real> QT(2*n,2*n, 0.0);
    x1[0] = pg1->a[0];
    for (int i=1; i<n; i++) {
      x1[i] = pg1->a[i];
      E1[i-1] = sqrt(pg1->b[i]);
    }
    for (int i=1; i<2*n; i++) {
      E2[i-1] = sqrt(tb2[i]);
      //E3[i-1] = sqrt(tb3[i]);
    }
    int info = 0;
    // struct timespec start, finish;
    // double elapsed = 0.0;
    // clock_gettime(CLOCK_MONOTONIC, &start);
    dstedc('N', n, &x1[0], &E1[0], QT.p, n, &info);
    dstedc('N', 2*n, &x2[0], &E2[0], QT.p, 2*n, &info);
    //dstedc('N', 2*n, &x3[0], &E3[0], QT.p, 2*n, &info);
    // clock_gettime(CLOCK_MONOTONIC, &finish);
    // elapsed += (finish.tv_sec-start.tv_sec);
    // elapsed += (finish.tv_nsec-start.tv_nsec)/1.0e9;
#if DDDD
    ofstream fp("xgrid_dd");
#else
    ofstream fp("xgrid");
#endif
    for (int i=n; i<2*n; i++) {
      fp << str(x1[i-n],0) << endl;
      fp << str(2*x2[i]*x2[i],0) << endl;
      //fp << str(-2*x3[3*n-1-i]*x3[3*n-1-i],0) << endl;
    }
 
    valarray<Real> p1(0.0,n);
    valarray<Real> tp2(0.0,2*n-1), tp3(0.0,2*n-1);
    // evaluate a point near -1, 1 by 3 methods and compare with qaudruple precision values to measure err.
    // method 1: original grid
    // method 2: -1->0, 1->1
    // method 3: -1 -> 1, 1 -> 0

    // phi recurrence
    if (0) {
      valarray<Real> val1(0.0,n);  // stores the value at the left end to normalize
      Real temp2 = sqrt(t2*0.5); // t2 = t1+1
      Real temp3 = sqrt(t3*0.5); // t3 = 1-t1
      // Real val = sqrt(PI/2)*min(Real(1.0),pow(t2,bb/2+0.25))*min(Real(1.0),pow(t3,aa/2+0.25));
      tp2[0] = 1.0/pg2->cc[0];
      tp2[1] = temp2*tp2[0]/sqrt(tb2[1]);
      
      // val2[0] = tp2[0];
      // val2[1] = 0.0;
      
      tp3[0] = 1.0/pg3->cc[0];
      tp3[1] = temp3*tp3[0]/sqrt(tb3[1]);
      
      // val3[0] = tp3[0];
      // val3[1] = tp3[0]/sqrt(tb3[1]);
      
      p1[0] = 1.0/pg1->cc[0];
      p1[1] = (t1-pg1->a[0])/pg1->cc[1];
      val1[0] = 1.0/pg1->cc[0];
      val1[1] = (-1.0-pg1->a[0])/pg1->cc[1];
      for (int i=1; i<2*n-2; i++) {
	tp2[i+1] = (temp2*tp2[i]-sqrt(tb2[i])*tp2[i-1])/sqrt(tb2[i+1]);
	// val2[i+1] = -sqrt(tb2[i])*val2[i-1]/sqrt(tb2[i+1]);
	tp3[i+1] = (temp3*tp3[i]-sqrt(tb3[i])*tp3[i-1])/sqrt(tb3[i+1]);
	// val3[i+1] = (val3[i]-sqrt(tb3[i])*val3[i-1])/sqrt(tb3[i+1]);
      }
      for (int i=1; i<n-1; i++) {
      	p1[i+1] = ((t1-pg1->a[i])*p1[i]-sqrt(pg1->b[i])*p1[i-1])/sqrt(pg1->b[i+1]);
	val1[i+1] = ((-1.0-pg1->a[i])*val1[i]-sqrt(pg1->b[i])*val1[i-1])/sqrt(pg1->b[i+1]);
      }
      temp = 1.0/pow(2,(aa+bb+1)/2);
#if DDDD
      ofstream fp("eval_dd");
#else
      ofstream fp("eval");
#endif
      Real sign = 1.0;
      for (int i=0; i<n; i++) {
	fp << str(p1[i]/val1[i],0) << endl;
	fp << str(temp*tp2[2*i]/val1[i],0) << endl;
        fp << str(sign*temp*tp3[2*i]/val1[i],0) << endl;
	sign *= -1.0;
	// fp << str(temp*tp2[2*i],0) << endl;
	// fp << str(sign*temp*tp3[2*i],0) << endl;
	// sign *= -1.0;
      }
      //#endif
    }
    
    // better way to see the evaluation error: evaluate the 3 methods on a Chebyshev grid on [0,x_1], where x_1 is the first zero of phi_n(x).  A 20-point chebyshev grid is probably enough.  You could normalize phi on that interval by its value at 0. And for the error, you could compute something like sqrt[ 1/20 \sum_{i=1}^20 err_i^2 ] to prevent the other methods from getting lucky with roundoff errors.
    if (1) {
      Real x, t, temp;
      int m = 20;
      valarray<Real> val(0.0,n);  // stores the value at the left end to normalize
      temp = 1.0/pow(2,(aa+bb+1)/2);
#if DDDD
      ofstream fp("eval_dd");
      string line;
      ifstream in("xx17");
      // ifstream inn("val");
      //getline( inn, line);
      //istringstream chuan(line);
      //chuan >> temp;
      //for (int i=0; i<n; i++) {
      //getline( inn, line);
      //istringstream chuan(line);
      //chuan >> val[i];
      //}
      FILE *fp2 = fopen("val","w");
      for (int i=0; i<m; i++) {
	getline( in, line);
	istringstream yan(line);
	yan >> t >> x;
	tp2[0] = 1.0/pg2->cc[0];
	tp2[1] = t*tp2[0]/sqrt(tb2[1]);
	p1[0] = 1.0/pg1->cc[0];
	p1[1] = (2*x-1-pg1->a[0])/pg1->cc[1];
	if (i==0) {
	  val[0] = 1.0/pg1->cc[0];
	  val[1] = (-1.0-pg1->a[0])/pg1->cc[1];
	  fprintf(fp2, "%s\n%s\n", str(val[0],0), str(val[1],0));
	}
	for (int j=1; j<2*n-2; j++) {
	  tp2[j+1] = (t*tp2[j]-sqrt(tb2[j])*tp2[j-1])/sqrt(tb2[j+1]);
	}
	for (int j=1; j<n-1; j++) {
	  p1[j+1] = ((2*x-1-pg1->a[j])*p1[j]-sqrt(pg1->b[j])*p1[j-1])/sqrt(pg1->b[j+1]);
	  if(i==0) {
	    val[j+1] = ((-1.0-pg1->a[j])*val[j]-sqrt(pg1->b[j])*val[j-1])/sqrt(pg1->b[j+1]);
	    fprintf(fp2, "%s\n", str(val[j+1],0));
	  }
	}
      	//getline( inn, line);
	for (int k=0; k<n; k++) {
	  fp << str(p1[k],0) << endl;
	  fp << str(temp*tp2[2*k],0) << endl;
	}
      }
      fclose(fp2);
#else
      chebyshev<Real> z(m);
      ofstream fp("eval");
      FILE *fp1 = fopen("xx17","w");
      //FILE *fp2 = fopen("val","w");
      //fprintf(fp2, "%s\n", str(temp,0));
      for(int j=0; j<m; j++) {
	t = z.xx[j+1]*x2[n];
	x = t*t;
	fprintf(fp1, "%s %s\n", str(t,0), str(x,0));
	
	tp2[0] = 1.0/pg2->cc[0];
	tp2[1] = t*tp2[0]/sqrt(tb2[1]);
	p1[0] = 1.0/pg1->cc[0];
	p1[1] = (2*x-1-pg1->a[0])/pg1->cc[1];
	//if (j==0) {
	//val[0] = 1.0/pg1->cc[0];
	//val[1] = (-1.0-pg1->a[0])/pg1->cc[1];
	//fprintf(fp2, "%s\n%s\n", str(val[0],0), str(val[1],0));
	//}
	for (int i=1; i<2*n-2; i++) {
	  tp2[i+1] = (t*tp2[i]-sqrt(tb2[i])*tp2[i-1])/sqrt(tb2[i+1]);
	}
	for (int i=1; i<n-1; i++) {
	  p1[i+1] = ((2*x-1-pg1->a[i])*p1[i]-sqrt(pg1->b[i])*p1[i-1])/sqrt(pg1->b[i+1]);
	  //if(j==0) {
	  //val[i+1] = ((-1.0-pg1->a[i])*val[i]-sqrt(pg1->b[i])*val[i-1])/sqrt(pg1->b[i+1]);
	  //fprintf(fp2, "%s\n", str(val[i+1],0));
	  //}
	}
	for (int i=0; i<n; i++) {
	  fp << str(p1[i],0) << endl;
	  fp << str(temp*tp2[2*i],0) << endl;
	}
	fp << endl;
      }
      //fclose(fp2);
      fclose(fp1);
#endif
    }
    
    // determine weights at xgrids and compare x/w with golub-welsch and measure the time
    // clock_gettime(CLOCK_MONOTONIC, &finish);
    // elapsed += (finish.tv_sec-start.tv_sec);
    // lapsed += (finish.tv_nsec-start.tv_nsec)/1.0e9;
    // cout << "time for new method: " << elapsed << endl;
    // fpw.close();
    // fpx.close();
    
// #if DDDD
//       ofstream fp("weight_dd");
// #else
//       ofstream fp("weight");
// #endif
//       for (int i=0; i<n; i++) {
// 	fp << str(w1[i],0) << endl;
// 	fp << str(w2[i],0) << endl;
// 	fp << str(w3[i],0) << endl;
//       }
//     }
    // elapsed = 0;
    // clock_gettime(CLOCK_MONOTONIC, &start);
    pg1->compute_xw(0, Phi.p);
    // clock_gettime(CLOCK_MONOTONIC, &finish);
    // elapsed = (finish.tv_sec-start.tv_sec);
    // elapsed += (finish.tv_nsec-start.tv_nsec)/1.0e9;
    // cout << "time for gw: " << elapsed << endl;
#if DDDD
    ofstream fpx("x_dd");
    ofstream fpw("w_dd");
#else
    ofstream fpw("w");
    ofstream fpx("x");
#endif
    valarray<Real> w1(0.0,n), w2(0.0,n), w3(0.0,n), w(0.0,n);
    temp = pow(2,aa+bb+1);
    Real t = 0.0;
    for (int j=0; j<n; j++) {
      t = x2[j+n];
      fpx << str(2*t*t,0) << endl;
      tp2[0] = 1.0/pg2->cc[0];
      tp2[1] = t*tp2[0]/sqrt(tb2[1]);
      for (int i=1; i<2*n-2; i++) {
	tp2[i+1] = (t*tp2[i]-sqrt(tb2[i])*tp2[i-1])/sqrt(tb2[i+1]);
      }
      for (int i=0; i<n; i++) {
	w[j] += tp2[2*i]*tp2[2*i];
      }
      w[j] = temp/w[j];
      fpx << str(pg1->xgrid[j],0) << endl;
      fpw << str(sqrt(w[j]),0) << endl;
      fpw << str(pg1->wgrid[j],0) << endl;
      // printf("%23s\n", str(w[j],0) );
    }
    
    // part 1 method 1 using linear terms to find tb
    if (1) {
      for (int i=0; i<n; i++){
	J(i,i) = pg2->a[i];
	if (i > 0){
	  J(i-1,i) = J(i,i-1) = sqrt(pg2->b[i]);
	}
      }
      mp_mat<Real> tJ(2*n,2*n,0.0);
      for (int i=1; i<2*n; i++) {
	tJ(i-1,i) = tJ(i,i-1) = sqrt(tb2[i]);
      }
      mp_mat<Real> tmp = tJ*tJ;
      err = 0.0;
      for (int i=0; i<n-1; i++) {
	err += pow(tmp(2*i,2*i+2)-J(i,i+1),2);
	// printf("%32s %32s\n", str(tmp(2*i,2*i+2),0), str(J(i,i+1),0));
      }
      printf("err between even submatrix of tJ^2 and J = %32s\n", str(sqrt(err),0));
    }
    
    
    // part 1 method 2 using two equation relations to find tb directly
    if (0) {
      valarray<Real> tb2_2(0.0,2*n);
      for (int i=0; i<n; i++) {
	tb2_2[2*i+1] = pg2->a[i]-tb2_2[2*i];
	if (i<n-1) {
	  tb2_2[2*i+2] = pg2->b[i+1]/tb2_2[2*i+1];
	}
      }
      err = 0.0;
      for (int i=0; i<2*n; i++) {
	temp = tb2_2[i]-tb2[i];
	err += temp*temp;
      }
      printf("err between tb's calculated by 2 methods = %32s\n", str(sqrt(err),0));
    }
}
