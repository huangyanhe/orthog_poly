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
    poly_grid<Real> *pg_bigger = new general_laguerre_grid<Real>(n+1,aa);
    // calculate hermite to check our method works correctly at the beginning
    poly_grid<Real> *pg1 = new general_hermite_grid<Real>(2*n+1,aa+0.5); // 0,1,2,...,2n-1
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
    
#if DDDD
    ofstream fp2("tb");
    Real tmp1;
    //string line1;
    //ifstream in1("b");
    cout << "tb generated for double precision usage" << endl;
    for  (int i=0; i<n; i++) {
      if (i>0) {
	//getline(in1, line1);
	//istringstream yan1(line1);
	//yan1 >> tmp1;
	//tb[2*i] = tmp1/tb[2*i-1];
	tb[2*i] = pg->b[i]/tb[2*i-1];
	fp2 << str(tb[2*i],0) << endl;
      }
      //getline(in1, line1);
      //istringstream yan1(line1);
      // yan1 >> tmp1;
      // tb[2*i+1] = tmp1-tb[2*i];
      tb[2*i+1] = pg->a[i]-tb[2*i];
      fp2 << str(tb[2*i+1],0) << endl;
    }
    err = 0.0;
    for (int i=0; i<2*n; i++) {
      err += (pg1->b[i]-tb[i])*(pg1->b[i]-tb[i]);
    }
    printf("method 2 err = %32s\n", str(sqrt(err),0));
#else
    if (0) {
      ofstream fp6("b");
      for (int i=0; i<n; i++) {
	fp6 << str(pg->a[i],0) << endl;
	if (i<n-1) {
	  fp6 << str(pg->b[i+1],0) << endl;
	}
      }
    }
    if (0) {
      string line;
      ifstream in("tb");
      for (int i=1; i<2*n; i++) {
	getline(in, line);
	istringstream yan(line);
	yan >> tb[i];
	//cout << str(tb[i],0) << endl;
      }
    }
    // only for testing the improved method
    if (0) {
      for  (int i=0; i<n; i++) {
	if (i>0) {
	  tb[2*i] = pg->b[i]/tb[2*i-1];
	}
	tb[2*i+1] = pg->a[i]-tb[2*i];
      }
    }
    // using hermite for tb
    if (1) {
      for (int i=1; i<2*n; i++) {
	tb[i] = pg1->b[i];
      }
    }
      
 
#endif
    

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
    valarray<Real> p(0.0,n+1);
    valarray<Real> tp(0.0,2*n+1); //the last 2 entries are for newton's method
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
    valarray<Real> x2(0.0,2*n), alpha(0.0, n+1), der(0.0, n), xx(0.0, n);
    valarray<Real> E1(0.0,n), E2(0.0,2*n);
    mp_mat<Real> QT(2*n,2*n, 0.0);
    for (int i=1; i<2*n; i++) {
      E2[i-1] = sqrt(tb[i]);
    }
    int info = 0;
    // eigenvalues in ascending orders are stored in X
    dstedc('N', 2*n, &x2[0], &E2[0], QT.p, 2*n, &info);
    for (int i=n; i<100+n; i++) {
      xx[i-n] = x2[i]*x2[i];
    }
    for (int i=100; i<n; i++) {
      xx[i] = pg->xgrid[i];
    }
    alpha[n] = 1.0;
    pg_bigger->eval_der(&alpha[0], n, &xx[0], &der[0], 1.0, 0, n+1);
    
#if DDDD
    ofstream fp("xgrid_dd");
    for (int i=0; i<n; i++) {
      fp << str(pg->xgrid[i],0) << endl;
      fp << str(pg->xgrid[i],0) << endl;
    }
#else
    ofstream fp("xgrid");
    for (int i=n; i<n+100; i++) {
      fp << str(pg->xgrid[i-n],0) << endl;
      t = x2[i];
      tp[0] = 1.0/pg->cc[0];
      tp[1] = t*tp[0]/sqrt(pg1->b[1]);
      for (int j=1; j<2*n; j++) {
	tp[j+1] = (t*tp[j]-sqrt(pg1->b[j])*tp[j-1])/sqrt(pg1->b[j+1]);
      }
      //cout << xx[i-n] << " " << tp[2*n] << " " << der[i-n] << endl;
      fp << str(xx[i-n]-tp[2*n]/der[i-n],0) << endl;
    }
    for(int i=100; i<n; i++) {
      t = pg->xgrid[i];
      p[0] = 1.0/pg->cc[0];
      p[1] = (t-pg->a[0])/pg->cc[1];
      for (int i=1; i<n; i++) {
	p[i+1] = ((t-pg->a[i])*p[i]-sqrt(pg->b[i])*p[i-1])/sqrt(pg->b[i+1]);
      }
      fp << str(pg->xgrid[i], 0) << endl;
      fp << str(xx[i]-p[n]/der[i], 0) << endl;
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
	tp[1] = t/(pg->cc[0]*sqrt(tb[1]));
	//tp[1] = t/(pg->cc[0]*sqrt(pg1->b[1]));
	p[0] = 1.0/pg->cc[0];
	p[1] = (x-pg->a[0])/pg->cc[1];
	for (int i=1; i<2*n-2; i++) {
	  tp[i+1] = (t*tp[i]-sqrt(tb[i])*tp[i-1])/sqrt(tb[i+1]);
	  //tp[i+1] = (t*tp[i]-sqrt(pg1->b[i])*tp[i-1])/sqrt(pg1->b[i+1]);
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
      ifstream fp1("chuan");
      if (!fp1.fail()) {
	cout << "eval1 created" << endl;
	ofstream myFile("eval1");
	string line;
	for(int j=0; j<m; j++) {
	  getline(fp1, line);
	  istringstream yan(line);
	  yan >> t >> x;
	  temp = pow(t,aa)*exp(-x/2);
	  tp[0] = 1.0/pg->cc[0];
	  tp[1] = t/(pg->cc[0]*sqrt(tb[1]));
	  p[0] = 1.0/pg->cc[0];
	  p[1] = (x-pg->a[0])/pg->cc[1];
	  for (int i=1; i<2*n-2; i++) {
	    tp[i+1] = (t*tp[i]-sqrt(tb[i])*tp[i-1])/sqrt(tb[i+1]);
	  }
	  for (int i=1; i<n-1; i++) {
	    p[i+1] = ((x-pg->a[i])*p[i]-sqrt(pg->b[i])*p[i-1])/sqrt(pg->b[i+1]);
	  }
	  for (int i=0; i<n; i++) {
	    myFile << str(temp*p[i],0) << endl;
	    myFile << str(temp*tp[2*i],0) << endl;
	  }
	}
      }
      else {
	cout << "eval created" << endl;
	ofstream fp("eval");
	FILE *fp1 = fopen("chuan","w");
	for (int j=0; j<m; j++) {
	  //t = z.xx[j+1]*x2[n];
	  //temp = pow(t,aa)*exp(-t*t/2);
	  x = z.xx[j+1]*z.xx[j+1]*pg->xgrid[0];
	  temp = pow(x,aa/2.0)*exp(-x/2);
	  //x = t*t;
	  t = sqrt(x);
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
      //cout << pg->xgrid[j] << endl;
      w[j] += pp*pp;
      Real p00;
      int temp, k1=0;
      Real absw;
      if (j<20) {
	t = x2[j+n];
	p0 = pp;
	pp = t*p0/sqrt(tb[1]);
	for (int i=2; i<2*n-1; i++) {
	  p00 = p0;
	  p0 = pp;
	  temp = 0;
	  pp = (t*p0 - sqrt(tb[i-1])*p00)/sqrt(tb[i]);
	  if (i%2==0) {
	    absw = abs(w[j]);
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
	tmp = -t*t;
	w[j] = 1.0/(w[j]*pow(t,2*aa));
      }
      else {
	t = pg->xgrid[j];
	p0 = pp;
	pp = (t-pg->a[0])/pg->cc[1];
	w[j] += pp*pp;
	for (int i=1; i<n-1; i++) {
	  p00 = p0;
	  p0 = pp;
	  temp = 0;
	  pp = ((t-pg->a[i])*p0-sqrt(pg->b[i])*p00)/sqrt(pg->b[i+1]);
	  absw = abs(w[j]);
	  if (absw > cM) { w[j] *= cm; ksum += 320; temp = 1;}
	  if (absw < cm) { w[j] *= cM; ksum -= 320; temp = -1;}
	  if (temp == 1) { mult_by_pow2(pp,-160); mult_by_pow2(p0,-160);}
	  if (temp == -1) { mult_by_pow2(pp,160); mult_by_pow2(p0,160);}
	  w[j] += pp*pp;
	}
	tmp = -t;
	w[j] = 1.0/(w[j]*pow(t,aa));
      }
      more = exp_(tmp,cm);
      // p[0] = 1.0/pg->cc[0];
      // p[1] = (x1[j]-pg->a[0])/pg->cc[1];
      // for (int i=1; i<n-1; i++) {
      //   p[i+1] = ((t-pg->a[i])*p[i]-sqrt(pg->b[i])*p[i-1])/sqrt(pg->b[i+1]);
      // }
      // for (int i=0; i<n; i++) {
      //   w[j] += p[i]*p[i];
      // }
      // printf("w[j]0 = %23s\n", str(w[j],0) );
      // w[j] = 1.0/(w[j]*pow(t,2*aa)); // if n is pretty big, xgrid value t would be pretty big, wo w[j] might underflow. Need to use cm, cM again. But I didn't implement here.
      for (int i=0; i<more; i++) {
	w[j] /= tmp;
	if (w[j]>cM) { w[j] *= cm; k1 += 320;}
      }
      w[j] = sqrt(w[j]);
      mult_by_pow2(w[j],(k1-ksum)/2);
      fp1 << str(pg->wgrid[j],0) << endl;
      fp1 << str(w[j],0) << endl;
      // printf("%23s\n",str(w[j],0) );
    }
#endif

     // to find the amplification factor of (recurrence matrix) in a, b and tb;
    Real a, b, c, d, a1, b1, a_temp, b_temp, c_temp, d_temp, sing1;
    Real x;
    Real sing_max;
    int sing_index;
    for (int k=0; k<30; k++) {
      x = 0.0000001*k;
      a = (x-pg->a[n-2])/sqrt(pg->b[n-1]);
      b = -sqrt(pg->b[n-2]/pg->b[n-1]);
      c = 1.0;
      d = 0.0;
      temp = a*a+b*b+c*c+d*d;
      sing1 = pow(a*a+b*b-c*c-d*d, 2)+4*pow(a*c+b*d,2);
      sing1 = sqrt((temp+sqrt(sing1))/2.0);
      sing_max = sing1;
      sing_index = n-2;
      for(int i=n-3; i>=0; i--) {
	a1 = (x-pg->a[i])/sqrt(pg->b[i+1]);
	b1 = -sqrt(pg->b[i]/pg->b[i+1]);
	c_temp = c*a1+d;
	d_temp = c*b1;
	a_temp = a*a1+b;
	b_temp = a*b1;
	a = a_temp;
	b = b_temp;
	c = c_temp;
	d = d_temp;
	temp = a*a+b*b+c*c+d*d;
	sing1 = pow(a*a+b*b-c*c-d*d, 2)+4*pow(a*c+b*d,2);
	sing1 = sqrt((temp+sqrt(sing1))/2.0);
	if (sing1 > sing_max) {
	  sing_index = i;
	  sing_max = sing1;
	}
      }
      printf("%23s, %d\n",str(sing_max,0), sing_index);
      
      t = sqrt(x);
      a = t/sqrt(pg1->b[2*n-1]);
      b = -sqrt(pg1->b[2*n-2]/pg1->b[2*n-1]);
      c = 1.0;
      d = 0.0;
      temp = a*a+b*b+c*c+d*d;
      sing1 = pow(a*a+b*b-c*c-d*d, 2)+4*pow(a*c+b*d,2);
      sing1 = sqrt((temp+sqrt(sing1))/2.0);
      sing_max = sing1;
      sing_index = 2*n-2;
      for(int i=2*n-3; i>=0; i--) {
	a1 = t/sqrt(pg1->b[i+1]);
	b1 = -sqrt(pg1->b[i]/pg1->b[i+1]);
	c_temp = c*a1+d;
	d_temp = c*b1;
	a_temp = a*a1+b;
	b_temp = a*b1;
	a = a_temp;
	b = b_temp;
	c = c_temp;
	d = d_temp;
	temp = a*a+b*b+c*c+d*d;
	sing1 = pow(a*a+b*b-c*c-d*d, 2)+4*pow(a*c+b*d,2);
	sing1 = sqrt((temp+sqrt(sing1))/2.0);
	if (sing1 > sing_max) {
	  sing_max = sing1;
	  sing_index = i;
	}
      }
      printf("%23s, %d\n",str(sing_max,0), sing_index);
    }
}
