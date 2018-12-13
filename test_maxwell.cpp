#include "poly_grid.h"
#include "real_def.h"
#include "chebyshev2.h"
#include "gamma_cx2.h"
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
    throw gen_err("\nusage:  ./test_maxwell n, nu");
  int n;
  int numIter;
  Real nu;
  int nu_int;
  Real t, x, temp;
  sscanf(argv[1], "%d", &n);
  sscanf(argv[2], "%d", &nu_int);
  str_to_real(argv[2], nu);
  Real c0 =  exp(0.5*lgamma(nu+0.5))/sqrt(2);
  int n1 = n+1;
  int ck = 0;
  int cM = 1024;
  Real cm = 1.0/cM;
  //generate a_n and b_n using algorithm from BALL
  FILE *file1 = fopen("abc","w");
  valarray<Real> g(0.0,n1+1), a(0.0,n1), b(0.0,n1);
  Real Y, C0, C1, C2, C3;
  g[0] = -nu/6.0;
  temp = exp(lgamma(nu+1)-lgamma(nu+0.5));
  g[1] = (2*nu+1)/3.0-g[0]-temp*temp;
  //g[1] = (2*nu+1)/3.0-g[0] - pow(tgamma(nu+1)/tgamma(nu+0.5), 2);
  
#if DDDD
  numIter = 16;
#else
  numIter = 8;
#endif
  for (int i=1; i<8; i++) {
    Y = 2*(i+nu);
    g[i+1] = (Y+1)/3-g[i] - pow((Y/6-g[i])*(Y/6-g[i])-nu*nu/4, 2)/((Y-1)/3-g[i]-g[i-1])/pow(Y/12+g[i], 2);
  }
  for (int i=9; i < n1+1; i++) {
    Y = 2*(i+nu);
    C0 = 1.0/36 - nu*nu/2.0;
    C1 = 23.0/432 - 11.0/12*nu*nu + 3.0/2*pow(nu,4);
    C2 = 1189.0/2592 - 409.0/48*nu*nu + 75.0/4*pow(nu,4) + 9*pow(nu,6);
    C3 = 196057.0/20736 - 153559.0/864*nu*nu + 7111.0/16*pow(nu,4) + 639.0/2*pow(nu,6) + 135.0/2*pow(nu,8);
    g[i] = C0/Y + C1/pow(Y,3) + C2/pow(Y,5) + C3/pow(Y,7);
  }
#if DDDD
  numIter = 30;
#else
  numIter = 10;
#endif
  
  for(int iter=0; iter<numIter; iter++) {
    for(int i=2; i<51; i++) {
      Y = 2*(i+nu);
      g[i] = (-81*pow(nu,4)+18*nu*nu*Y*Y+216*pow(g[i],3)*(6*g[i-1]+6*g[i+1]+Y)+9*g[i]*g[i]*(-16+72*nu*nu+24*g[i-1]*(-2+6*g[i+1]-Y)-24*g[i+1]*(-2+Y)-23*Y*Y)-Y*Y*(1+3*g[i+1]*(Y-1)+g[i-1]*(3-9*g[i+1]+3*Y))-9*g[i]*Y*(g[i+1]*(7*Y-8)+g[i-1]*(8-24*g[i+1]+7*Y)))/(Y*nu*nu/6-7*Y*Y*Y/216+Y/54)/1296;
    }
  }

  for(int i=0; i<n1; i++){
    a[i] = sqrt((2*i+2*nu+1)/3-g[i]-g[i+1]);
    b[i] = (i+nu)/6+g[i];
    if (i) {
      c0 *= sqrt(b[i]);
      if ( c0 > cM ) { c0 *= cm; ck += 320; }
      if ( c0 < cm ) { c0 *= cM; ck -= 320; }
    }
    fprintf(file1, "%s %s %s %d\n", str(a[i],0), str(b[i],0), str(c0,0), ck);
  }
  fclose(file1);


  
  // calculating tb in Fred-Lew-Quarles way
  c0 = exp(0.5*lgamma(nu+0.5))/sqrt(2);
  //c0 = sqrt(tgamma(nu+0.5)*0.5);
  ck = 0;
  Real rho = 4*nu+1;
  n1 = 2*n+100;//n1 = 2*n+1; // n1 is the right boundary.
  valarray<Real> tb(0.0, n1);
  //printf("%d,%d, %s\n", 1,1, str(tb[1],0));
  //initialize tb with asymptotic values
  valarray<Real> gamm(0.0, n1), df(0.0, n1-2);
  mp_mat<Real> ddf(n1-2,n1-2,0.0);
  for(int i=1; i<n1; i++) {
    if (i%2) {
      gamm[i] = sqrt(i+rho);
      tb[i] = gamm[i]/sqrt(3);
    } else {
      gamm[i] = sqrt(i);
      tb[i] = gamm[i]/sqrt(3);
    }
  }
  //tb[1] = tgamma((3+rho)/4)/tgamma((1+rho)/4)*2;
  tb[1] = 2*exp(lgamma((3+rho)/4)-lgamma((1+rho)/4));


  for(int i=2; i<20; i++) {
    tb[i] = gamm[i-1]*gamm[i-1]/tb[i-1]-tb[i-1]-tb[i-2];
  }
  
  int info=0;
  numIter = 6;
  valarray<int> ipiv(0, n1-2);
  FILE *file = fopen("abc_associate","w");
  fprintf(file,"%s %s %s %d\n", str(0,0), str((0.0),0), str(c0, 0), 0);
  for(int iter=0; iter<numIter; iter++) {
    for(int i=0; i<n1-2; i++) {
      if (i<n1-3) {
	ddf(i,i+1) = 1;
	ddf(i+1,i) = 1;
      }
      ddf(i,i) = 1+gamm[i+1]*gamm[i+1]/tb[i+1]/tb[i+1];
      df[i] = -gamm[i+1]*gamm[i+1]/tb[i+1]+tb[i]+tb[i+1]+tb[i+2];
    }
    dgetrf(n1-2, n1-2, ddf.p, n1-2, &ipiv[0], &info);
     cout << "haha2" << endl;
    dgetrs('N', n1-2, 1, ddf.p, n1-2, &ipiv[0] , &df[0], n1-2, &info);
     cout << "haha3" << endl;
    for(int i=1; i<n1-1; i++) {
      temp = tb[i] - df[i-1];
      //while ( (temp < 0) || (temp > gamm[i])) {
      //	df[i-1] /= 2;
      //	temp = tb[i] - df[i-1];
      //	cout << tb[i] << " " << df[i-1] << endl;
      //	cout << "hwat" << endl;
      //}
      tb[i] = temp;
      if(iter == numIter-1 && i < 2*n) {
	c0 *= sqrt(temp/2);
	if ( c0 > cM ) { c0 *= cm; ck += 320; }
	if ( c0 < cm ) { c0 *= cM; ck -= 320; }
	fprintf(file, "%s %s %s %d\n", str(0,0), str(temp/2,0), str(c0,0), ck);
      }
    } 
  }

  fclose(file);
    cout << "I'm here111 " << endl;
  poly_grid<Real> *pg = new plasma_grid<Real>(n,"abc",nu_int); // 0,1,...,n-1
  poly_grid<Real> *pg_bigger = new plasma_grid<Real>(n+1,"abc",nu_int);
  poly_grid<Real> *pg1 = new associate_plasma_grid<Real>(2*n,"abc_associate", nu_int);

  cout << "I'm here " << endl;

  
  valarray<Real> p(0.0,n+1);
  valarray<Real> tp(0.0,2*n+1);
#if DDDD
  ofstream fp3("tb");
  Real tmp1;
  cout << "tb generated for double precision usage" << endl;
  if (0) {
    for  (int i=0; i<n; i++) {
      if (i>0) {
	tb[2*i] = pg->b[i]/tb[2*i-1];
	fp3 << str(tb[2*i],0) << endl;
      }
      tb[2*i+1] = pg->a[i]-tb[2*i];
      fp3 << str(tb[2*i+1],0) << endl;
    }
  }
  if (1) {
    for (int i=1; i<2*n; i++) {
        tb[i] = pg1->b[i];
	fp3 << str(tb[i], 0) << endl;
      }
    }
#else
    if (0) {
      string line;
      ifstream in("tb");
      for (int i=1; i<2*n; i++) {
	getline(in, line);
	istringstream yan(line);
	yan >> tb[i];
	cout << str(tb[i],0) << endl;
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
    // use tb calculated directly
    if (1) {
      for (int i=1; i<2*n; i++) {
	tb[i] = pg1->b[i];
      }
    }
#endif
    
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
    // eigenvalues in ascending orders are stored in X
    dstedc('N', 2*n, &x2[0], &E2[0], QT.p, 2*n, &info);
    for (int i=n; i<400+n; i++) {
      xx[i-n] = x2[i]*x2[i];
    }
    for (int i=400; i<n; i++) {
      xx[i] = pg->xgrid[i];
    }
    alpha[n] = 1.0;
    pg_bigger->eval_der(&alpha[0], n, &xx[0], &der[0], 1.0, 0, n+1);


    
#if DDDD
    ofstream fp4("xgrid_dd");
    for (int i=0; i<n; i++) {
      fp4 << str(pg->xgrid[i],0) << endl;
      fp4 << str(pg->xgrid[i],0) << endl;
    }
#else
    ofstream fp4("xgrid");
    for (int i=n; i<n+400; i++) {
      fp4 << str(pg->xgrid[i-n],0) << endl;
      t = x2[i];
      tp[0] = 1.0/pg->cc[0];
      tp[1] = t/(pg->cc[0]*sqrt(tb[1]));
      for (int i=1; i<2*n; i++) {
	tp[i+1] = (t*tp[i]-sqrt(tb[i])*tp[i-1])/sqrt(tb[i+1]);
      }
      cout << tp[2*n] << endl;
      fp4 << str(xx[i-n]-tp[2*n]/der[i-n],0) << endl;
    }
    for(int i=400; i<n; i++) {
      t = pg->xgrid[i];
      p[0] = 1.0/pg->cc[0];
      p[1] = (t-pg->a[0])/pg->cc[1];
      cout << p[0] << " " << p[1] << endl;
      for (int i=1; i<n; i++) {
	p[i+1] = ((t-pg->a[i])*p[i]-sqrt(pg->b[i])*p[i-1])/sqrt(pg->b[i+1]);
	//cout << "the " << i+1 << " is " << p[i+1] << endl;
      }
      fp4 << str(pg->xgrid[i], 0) << endl;
      cout << xx[i] << " " << p[n] << " " << der[i] << endl;
      fp4 << str(xx[i]-p[n]/der[i], 0) << endl;
    }
#endif

    int m = 20;
#if DDDD
    ofstream fp("eval_dd");
    string line;
    ifstream in("chuan");
    for (int i=0; i<m; i++) {
      getline(in, line);
      istringstream is(line);
      is >> t >> x;
      temp = pow(x,nu)*exp(-x*x/2);
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
	temp = pow(x,nu)*exp(-x*x/2);
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
	x = z.xx[j+1]*z.xx[j+1]*pg->xgrid[0];
	temp = pow(x,nu)*exp(-x*x/2);
	t = sqrt(x);
	fprintf(fp1, "%s %s\n", str(t,0), str(x,0));
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
	  fp << str(temp*p[i],0) << endl;
	  fp << str(temp*tp[2*i],0) << endl;
	}
      }
    }
#endif
    
    
    // determine weights at xgrid
    valarray<Real> w1(0.0,n), w2(0.0,n), w(0.0,n);
    Real more = 0;
    Real p0 = 0.0;
    Real pp = 1.0;
    Real tmp = 0.0;
    int kp = 0;
    int ksum = 0;
#if DDDD
    ofstream fp2("wgrid_dd");
    for (int i=0; i<n; i++) {
      fp2 << str(pg->wgrid[i],0) << endl;
      fp2 << str(pg->wgrid[i],0) << endl;
    }
#else
    ofstream fp2("wgrid");
    for (int j=0; j<n; j++) {
      p0 = 0.0;
      pp = 1.0/pg->cc[0];
      tmp = 0.0;
      kp = 0;
      ksum = 0;
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
	    if (temp == 1) { mult_by_pow2(pp,-160); mult_by_pow2(p0,-160);}
	    if (temp == -1) { mult_by_pow2(pp,160); mult_by_pow2(p0,160);}
	    w[j] += pp*pp;
	  }
	}
	tmp = -t*t*t*t;
	w[j] = 1.0/(w[j]*pow(t,4*nu));
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
	tmp = -t*t;
	w[j] = 1.0/(w[j]*pow(t,nu*2));
      }
      more = exp_(tmp,cm);
      for (int i=0; i<more; i++) {
	w[j] /= tmp;
	if (w[j]>cM) { w[j] *= cm; k1 += 320;}
      }
      w[j] = sqrt(w[j]);
      mult_by_pow2(w[j],(k1-ksum)/2);
      fp2 << str(pg->wgrid[j],0) << endl;
      fp2 << str(w[j],0) << endl;
      // printf("%23s\n",str(w[j],0) );
    }
#endif


    
}
