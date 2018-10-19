#include "poly_grid.h"
#include "real_def.h"
#include "chebyshev2.h"
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

Real exp_(Real &x, Real cm)
{
  int count = 1;
  while (exp(x/count)<cm) {
    count += 1;
  }
  x = exp(x/count);
  return count;
}


int main(int argc, char *argv[])
{
    if (argc < 3)
        throw gen_err("\nusage:  ./test_jacobi n, aa");
    int n;
    Real aa;
    Real t;
    sscanf(argv[1], "%d", &n);
    str_to_real(argv[2], aa);
    //str_to_real(argv[3], t);
    poly_grid<Real> *pg = new general_laguerre_grid<Real>(n,aa); // 0,1,...,n-1
    // calculate hermite to check our method works correctly at the beginning
    poly_grid<Real> *pg1 = new general_hermite_grid<Real>(2*n,aa+0.5); // 0,1,2,...,2n-1
    Real cM = pg->cM;
    Real cm = pg->cm;
    
    mp_mat<Real> J(n,n,0.0); 
    for (int i=0; i<n; i++){
      J(i,i) = pg->a[i];
      if (i > 0){
    	J(i-1,i) = J(i,i-1) = sqrt(pg->b[i]);
      }
    }
    Real err = 0.0;
    valarray<Real> tb(0.0,2*n);
    // part 1 method 1: d0 and d1 overflows when n is too large
    if (0) {
      valarray<Real> d0(0.0,n+1);
      valarray<int> dk0(0,n+1);
      valarray<Real> d1(0.0,n+1);
      valarray<int> dk1(0,n+1);
      Real dd0,dd1,dd1_;
      d0[0] = 1.0;
      d0[1] = -pg->a[0];
      d1[0] = 0.0;
      d1[1] = 1.0;
      for (int i=2; i<=n; i++) {
	dk0[i] = dk1[i] = dk0[i-1];
	dd0 = d0[i-2];
	mult_by_pow2(dd0,dk0[i-2]-dk0[i-1]);
	d0[i] = -pg->a[i-1]*d0[i-1]-pg->b[i-1]*dd0;
	dd1 = d1[i-1];
	mult_by_pow2(dd1,dk1[i-1]-dk0[i-1]);
	dd1_= d1[i-2];
	mult_by_pow2(dd1_,dk1[i-2]-dk0[i-1]);
        // d0[i] = -pg->a[i-1]*d0[i-1]-pg->b[i-1]*d0[i-2];
	d0[i] = -pg->a[i-1]*d0[i-1]-pg->b[i-1]*dd0;
        // d1[i] = d0[i-1]-pg->a[i-1]*d1[i-1]-pg->b[i-1]*d1[i-2];
	d1[i] = d0[i-1]-pg->a[i-1]*dd1-pg->b[i-1]*dd1_;
	if (abs(d0[i])>cM) { d0[i] *= cm; dk0[i] += 320;}
	if (abs(d1[i])>cM) { d1[i] *= cm; dk1[i] += 320;}
	if (abs(d0[i])<cm) { d0[i] *= cM; dk0[i] -= 320;}
	if (abs(d1[i])<cm) { d1[i] *= cM; dk1[i] -= 320;}
      }
      // for (int i=0; i<n+1; i++)
      // 	printf("%d, d0 = %23s, d1=%23s\n",i,str(d0[i],0),str(d1[i],0));
      for (int i=0; i<n; i++) {
	tb[2*i+1] = -d0[i+1]/d0[i];
	mult_by_pow2(tb[2*i+1], dk0[i+1]-dk0[i]);
      }
      for (int i=1; i<n; i++) {
	dd0 = d0[i];
	mult_by_pow2(dd0,dk0[i]-dk1[i]);
	dd1 = d1[i+1];
	mult_by_pow2(dd1,dk1[i+1]-dk1[i]);
	dd1_ = d1[i-1];
	mult_by_pow2(dd1_,dk1[i-1]-dk1[i]);
	tb[2*i] = (dd0 - dd1 - tb[2*i+1]*d1[i])/(d1[i]+tb[2*i-1]*dd1_);
	// tb[2*i] = (d0[i]-d1[i+1]-tb[2*i+1]*d1[i])/(d1[i]+tb[2*i-1]*d1[i-1]);
      }
      for (int i=1; i<2*n; i++) {
	err += (pg1->b[i]-tb[i])*(pg1->b[i]-tb[i]);
      }
      printf("method 1 err = %32s\n", str(sqrt(err),0));
    }
    // part 1 method 2 and it always gives accurate tb's
    if (1) {
      for  (int i=0; i<n; i++) {
	if (i>0) {
	  tb[2*i] = pg->b[i]/tb[2*i-1];
	  // cout << 2*i << " " << tb[2*i] << endl;
	}
	tb[2*i+1] = pg->a[i]-tb[2*i];
	// cout << 2*i+1 << " " << tb[2*i+1] << endl;
      }
      err = 0.0;
      for (int i=0; i<2*n; i++) {
	err += (pg1->b[i]-tb[i])*(pg1->b[i]-tb[i]);
      }
      printf("method 2 err = %32s\n", str(sqrt(err),0));
    }

    // finding A_n/B_n where A's and B's are solutions with initial value 1,0 and 0,1
    if (0) {
      Real x = 0.01;
      Real tmp1;
      for(int j=40; j < 20000; j++) {
	tmp1 = 0.0;
	for(int i=j; i>1; i--) {
	  tmp1 = -sqrt(tb[i-1]/tb[i])/(x/sqrt(tb[i])+tmp1);
	}
	cout << j << " " << str(tmp1,0) << endl;
      }
    }
    
    // evaluating tphi_2n(t) and compare with phi_n(t*t)
    valarray<Real> p(0.0,n);
    valarray<Real> tp(0.0,2*n-1);
    Real temp = 0.0;
    if (0) {
      tp[0] = 1.0/pg->cc[0];
      tp[1] = t/(pg->cc[0]*sqrt(tb[1]));
      p[0] = 1.0/pg->cc[0];
      p[1] = (t*t-pg->a[0])/pg->cc[1];
      for (int i=1; i<2*n-2; i++) {
	tp[i+1] = (t*tp[i]-sqrt(tb[i])*tp[i-1])/sqrt(tb[i+1]);
      }
      for (int i=1; i<n-1; i++) {
	p[i+1] = ((t*t-pg->a[i])*p[i]-sqrt(pg->b[i])*p[i-1])/sqrt(pg->b[i+1]);
      }
#if DDDD
      ofstream fp("eval_dd");
#else
      ofstream fp("eval");
#endif
      for (int i=0; i<n; i++) {
	fp << str(temp*p[i],0) << endl;
	fp << str(temp*tp[2*i],0) << endl;
      }
    }


    // Golub-Welsch to calculate xgrid and wgrid
    mp_mat<Real> Phi(n,n);
    pg->compute_xw(0, Phi.p);
    
    // 2 xgrids 
    valarray<Real> x2(0.0,2*n);
    valarray<Real> E1(0.0,n), E2(0.0,2*n);
    mp_mat<Real> QT(2*n,2*n, 0.0);
    for (int i=1; i<2*n; i++) {
      E2[i-1] = sqrt(tb[i]);
    }
    int info = 0;
    // eigenvalues in ascending orders are stored in X
    dstedc('N', 2*n, &x2[0], &E2[0], QT.p, 2*n, &info);
#if DDDD
    ofstream fp("xgrid_dd");
    for (int i=0; i<n; i++) {
      fp << str(pg->xgrid[i],0) << endl;
      fp << str(pg->xgrid[i],0) << endl;
    }
#else
    ofstream fp("xgrid");
    for (int i=n; i<2*n; i++) {
      fp << str(pg->xgrid[i-n],0) << endl;
      fp << str(x2[i]*x2[i],0) << endl;
    }
#endif


// evaluate using Chebyshev polynomial
    if (1) {
      Real x, t;
      int m = 20;
#if DDDD
      ofstream fp("eval_dd");
      string line;
      ifstream in("chuan");
      for (int i=0; i<m; i++) {
	getline(in, line);
	istringstream is(line);
	is >> t >> x;
	temp = pow(t,aa)*exp(-t*t/2);
        tp[0] = 1.0/pg->cc[0];
	//tp[1] = t/(pg->cc[0]*sqrt(tb[1]));
	tp[1] = t/(pg->cc[0]*sqrt(pg1->b[1]));
	p[0] = 1.0/pg->cc[0];
	p[1] = (x-pg->a[0])/pg->cc[1];
	for (int i=1; i<2*n-2; i++) {
	  //tp[i+1] = (t*tp[i]-sqrt(tb[i])*tp[i-1])/sqrt(tb[i+1]);
	  tp[i+1] = (t*tp[i]-sqrt(pg1->b[i])*tp[i-1])/sqrt(pg1->b[i+1]);
	}
	for (int i=1; i<n-1; i++) {
	  p[i+1] = ((x-pg->a[i])*p[i]-sqrt(pg->b[i])*p[i-1])/sqrt(pg->b[i+1]);
	}
	for (int i=0; i<n; i++) {
	  fp << str(temp*p[i],0) << endl;
	  fp << str(temp*tp[2*i],0) << endl;
	}
      }
#else
      chebyshev<Real> z(m);
      ofstream fp("eval");
      FILE *fp1 = fopen("chuan","w");
      for (int j=0; j<m; j++) {
	t = z.xx[j+1]*x2[n];
	temp = pow(t,aa)*exp(-t*t/2);
	x = t*t;
	fprintf(fp1, "%s %s\n", str(t,0), str(x,0));
	tp[0] = 1.0/pg->cc[0];
	//tp[1] = t/(pg->cc[0]*sqrt(tb[1]));
	tp[1] = t/(pg->cc[0]*sqrt(pg1->b[1]));
	p[0] = 1.0/pg->cc[0];
	p[1] = (x-pg->a[0])/pg->cc[1];
	for (int i=1; i<2*n-2; i++) {
	  //tp[i+1] = (t*tp[i]-sqrt(tb[i])*tp[i-1])/sqrt(tb[i+1]);
	  tp[i+1] = (t*tp[i]-sqrt(pg1->b[i])*tp[i-1])/sqrt(pg1->b[i+1]);
	}
	for (int i=1; i<n-1; i++) {
	  p[i+1] = ((x-pg->a[i])*p[i]-sqrt(pg->b[i])*p[i-1])/sqrt(pg->b[i+1]);
	}
	for (int i=0; i<n; i++) {
	  fp << str(temp*p[i],0) << endl;
	  fp << str(temp*tp[2*i],0) << endl;
	}
      }
#endif
    }


    // determine weights at xgrid
    valarray<Real> w1(0.0,n), w2(0.0,n), w(0.0,n);
    Real more = 0;
    Real p0 = 0.0;
    Real pp = 1.0;
    Real tmp = 0.0;
    int kp = 0;
    int ksum = 0;
#if DDDD
    ofstream fp1("wgrid_dd");
    for (int i=0; i<n; i++) {
      fp1 << str(pg->wgrid[i],0) << endl;
      fp1 << str(pg->wgrid[i],0) << endl;
    }
#else
    ofstream fp1("wgrid");
    for (int j=0; j<n; j++) {
      p0 = 0.0;
      pp = 1.0/pg->cc[0];
      tmp = 0.0;
      kp = 0;
      ksum = 0;
      cout << pg->xgrid[j] << endl;
      w[j] += pp*pp;
      t = x2[j+n];
      p0 = pp;
      pp = t*p0/sqrt(tb[1]);
      for (int i=2; i<2*n-1; i++) {
	Real p00 = p0;
	p0 = pp;
	int temp = 0;
	pp = (t*p0 - sqrt(tb[i-1])*p00)/sqrt(tb[i]);
	if (i%2==0) {
	  Real absw = abs(w[j]);
	  if (absw > cM) { w[j] *= cm; ksum += 320; temp = 1;}
	  if (absw < cm) { w[j] *= cM; ksum -= 320; temp = -1;}
	  // if (i%32==0) { // chosing 16, some weights(e.g. k=191) will be nan since kp=320,ksum=0
	  //   Real absp = abs(pp) + abs(p0);
	  //   if (absp>cM) { pp *= cm; p0 *= cm; kp += 320;}
	  //   if (absp<cm) { pp *= cM; p0 *= cM; kp -= 320;}
	  //   absp = abs(w[j]);
	  //   if (absp > cM*cM) { w[j] *= cm*cm; ksum += 2*320;}
	  //   if (absp < cm*cm) { w[j] *= cM*cM; ksum -= 2*320;}
	  // }
	  if (temp == 1) { mult_by_pow2(pp,-160); mult_by_pow2(p0,-160);}
	  if (temp == -1) { mult_by_pow2(pp,160); mult_by_pow2(p0,160);}
	    w[j] += pp*pp;
	}
      }
      // p[0] = 1.0/pg->cc[0];
      // p[1] = (x1[j]-pg->a[0])/pg->cc[1];
      // for (int i=1; i<n-1; i++) {
      //   p[i+1] = ((t-pg->a[i])*p[i]-sqrt(pg->b[i])*p[i-1])/sqrt(pg->b[i+1]);
      // }
      // for (int i=0; i<n; i++) {
      //   w[j] += p[i]*p[i];
      // }
      int k1 = 0;
      tmp = -t*t;
      more = exp_(tmp,cm);
      // printf("w[j]0 = %23s\n", str(w[j],0) );
      w[j] = 1.0/(w[j]*pow(t,2*aa)); // if n is pretty big, xgrid value t would be pretty big, wo w[j] might underflow. Need to use cm, cM again. But I didn't implement here.
      for (int i=0; i<more; i++) {
	w[j] /= tmp;
	if (w[j]>cM) { w[j] *= cm; k1 += 320;}
      }
      mult_by_pow2(w[j],k1-ksum);
      fp1 << str(pg->wgrid[j],0) << endl;
      fp1 << str(sqrt(w[j]),0) << endl;
      // printf("%23s\n",str(w[j],0) );
    }
#endif
}
