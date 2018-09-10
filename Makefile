include Sysdep

headers = $(shell ls ../*.h)

testCvK : CvK.o mp_mat.o mp_mat_double.o str_double.o gammaETD.o \
	qr_class.o qr_class_cx.o poly_grid.o chebyshev.o fftjw2.o fft_template.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

kimura_green5e : kimura_green5e.o kimura_time.o nugrid3.o poly_grid.o str_double.o \
	mp_mat.o chebyshev2.o fftjw2.o hyperF.o gamma_cx2.o irk_jw2.o qr_class.o \
	fft_template.o gammaETD2.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

kimura_green5d : kimura_green5d.o kimura_time.o nugrid3.o poly_grid.o str_double.o \
	mp_mat.o chebyshev2.o fftjw2.o hyperF.o gamma_cx2.o irk_jw2.o qr_class.o \
	fft_template.o gammaETD2.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

kimura_green5c : kimura_green5c.o kimura_time.o nugrid3.o poly_grid.o str_double.o \
	mp_mat.o chebyshev2.o fftjw2.o hyperF.o gamma_cx2.o irk_jw2.o qr_class.o \
	fft_template.o gammaETD2.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

kimura_green5b : kimura_green5b.o kimura_time.o nugrid2.o poly_grid.o str_double.o \
	mp_mat.o chebyshev2.o fftjw2.o hyperF.o gamma_cx2.o irk_jw2.o qr_class.o \
	fft_template.o gammaETD2.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

kimura_green5a : kimura_green5a.o kimura_time.o nugrid2.o poly_grid.o str_double.o \
	mp_mat.o chebyshev2.o fftjw2.o hyperF.o gamma_cx2.o irk_jw2.o qr_class.o \
	fft_template.o gammaETD2.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

kimura_green4a : kimura_green4a.o nugrid2.o poly_grid.o str_double.o mp_mat.o \
	chebyshev2.o fftjw2.o hyperF.o gamma_cx2.o irk_jw2.o qr_class.o \
	fft_template.o gammaETD2.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

kimura_green3a : kimura_green3a.o nugrid2.o poly_grid.o str_double.o mp_mat.o \
	chebyshev2.o fftjw2.o hyperF.o gamma_cx2.o irk_jw2.o qr_class.o \
	fft_template.o gammaETD2.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

kimura_green3a_dd : kimura_green3a.ddo nugrid2.ddo poly_grid.ddo str_dd.ddo \
	mp_mat.ddo chebyshev2.ddo fftjw2.ddo hyperF.ddo gamma_cx2.ddo irk_jw2.ddo \
	qr_class.ddo fft_template.ddo gammaETD2.ddo lapack_dd.ddo dd_utils.ddo \
	dstedc_dd.ddo mp_mat_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_gamma_algorithm : test_gamma_algorithm.mpo str_mpfr.o fftjw2.mpo fft_template.mpo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS) -lgmpxx

kimura_green2b : kimura_green2b.o nugrid2.o poly_grid.o str_double.o mp_mat.o \
	chebyshev2.o fftjw2.o hyperF.o gamma_cx2.o irk_jw2.o qr_class.o \
	fft_template.o gammaETD2.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

kimura_green2a : kimura_green2a.o nugrid2.o poly_grid.o str_double.o mp_mat.o \
	chebyshev2.o fftjw2.o hyperF.o gamma_cx2.o irk_jw2.o qr_class.o \
	fft_template.o gammaETD2.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

kimura_green2a_dd : kimura_green2a.ddo nugrid2.ddo poly_grid.ddo str_dd.ddo \
	mp_mat.ddo chebyshev2.ddo fftjw2.ddo hyperF.ddo gamma_cx2.ddo irk_jw2.ddo \
	qr_class.ddo fft_template.ddo gammaETD2.ddo lapack_dd.ddo dd_utils.ddo \
	dstedc_dd.ddo mp_mat_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

kimura_green2 : kimura_green2.o nugrid2.o poly_grid.o str_double.o mp_mat.o \
	chebyshev2.o fftjw2.o hyperF.o gamma_cx2.o irk_jw2.o qr_class.o \
	fft_template.o gammaETD2.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

kimura_green : kimura_green.o nugrid2.o poly_grid.o str_double.o mp_mat.o \
	chebyshev2.o fftjw2.o hyperF.o gamma_cx2.o irk_jw2.o qr_class.o \
	fft_template.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

kimura_green_dd : kimura_green.ddo nugrid2.ddo poly_grid.ddo str_dd.ddo \
	mp_mat.ddo chebyshev2.ddo fftjw2.ddo qr_class.ddo fft_template.ddo \
	dd_utils.ddo lapack_dd.ddo dstedc_dd.ddo hyperF.ddo gamma_cx2.ddo \
	irk_jw2.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

# need to implement MPFR in irk_jw2 before this can work:
kimura_green_mpfr : kimura_green.mpo nugrid2.mpo poly_grid.mpo str_mpfr.mpo \
	mp_mat.mpo chebyshev2.mpo fftjw2.mpo qr_class.mpo fft_template.mpo \
	dd_utils.mpo lapack_mpfr.mpo dstedc_mpfr.mpo hyperF.mpo gamma_cx2.mpo \
	irk_jw2.mpo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_hyperF : test_hyperF.o hyperF.o str_double.o str_dd.o gamma_cx.o irk_jw2.o \
	poly_grid.o mp_mat.o qr_class.o chebyshev2.o fftjw2.o fft_template.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_hyperF_dd : test_hyperF.ddo hyperF.ddo str_dd.ddo str_mpfr.o gamma_cx.ddo \
	irk_jw2.ddo poly_grid.ddo mp_mat.ddo qr_class.ddo lapack_dd.ddo \
	dstedc_dd.ddo dd_utils.ddo chebyshev2.ddo fftjw2.ddo fft_template.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_hyperF_mpfr : test_hyperF.mpo hyperF.mpo str_mpfr.mpo gamma_cx.mpo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_gamma_cx : test_gamma_cx.o gamma_cx.o str_double.o str_dd.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_gamma_cx_dd : test_gamma_cx.ddo gamma_cx.ddo str_dd.ddo str_mpfr.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_gamma_cx_mpfr : test_gamma_cx.mpo gamma_cx.mpo str_mpfr.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_gamma_cx2 : test_gamma_cx2.o gamma_cx2.o str_double.o str_dd.o fftjw2.o fft_template.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS) -lgmpxx

test_gamma_cx2_dd : test_gamma_cx2.ddo gamma_cx2.ddo str_dd.ddo \
	str_mpfr.ddo dd_utils.ddo fftjw2.ddo fft_template.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS) -lgmpxx

test_gamma_cx2_mpfr : test_gamma_cx2.mpo gamma_cx2.mpo str_mpfr.mpo poly_grid.mpo \
	fftjw2.mpo fft_template.mpo lapack_mpfr.mpo mp_mat.mpo dstedc_mpfr.mpo \
	mp_mat_mpfr.mpo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS) -lgmpxx

kimura_setd5c : kimura_setd5c.o poly_grid.o gammaETD2.o mp_mat.o str_double.o \
	chebyshev2.o fftjw2.o qr_class.o levmar2.o mp_mat_double.o \
	fft_template.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

kimura_setd5a : kimura_setd5a.o poly_grid.o gammaETD2.o mp_mat.o str_double.o \
	chebyshev2.o fftjw2.o qr_class.o levmar2.o mp_mat_double.o \
	fft_template.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

kimura_setd5b : kimura_setd5b.o poly_grid.o gammaETD2.o mp_mat.o str_double.o \
	chebyshev2.o fftjw2.o qr_class.o levmar2.o mp_mat_double.o \
	fft_template.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

kimura_setd5 : kimura_setd5.o poly_grid.o gammaETD2.o mp_mat.o str_double.o \
	chebyshev2.o fftjw2.o qr_class.o levmar2.o mp_mat_double.o \
	fft_template.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

kimura_setd4 : kimura_setd4.o poly_grid.o gammaETD.o mp_mat.o str_double.o \
	chebyshev.o fftjw.o qr_class.o levmar2.o mp_mat_double.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

kimura_setd3 : kimura_setd3.o poly_grid.o gammaETD.o mp_mat.o str_double.o \
	chebyshev.o fftjw.o qr_class.o levmar2.o mp_mat_double.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

kimura_setd2 : kimura_setd2.o poly_grid.o gammaETD.o mp_mat.o str_double.o \
	chebyshev.o fftjw.o qr_class.o levmar2.o mp_mat_double.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

kimura_setd : kimura_setd.o poly_grid.o gammaETD.o mp_mat.o str_double.o \
	chebyshev.o fftjw.o qr_class.o levmar2.o mp_mat_double.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_setd : test_setd.o poly_grid.o gammaETD.o mp_mat.o str_double.o \
	chebyshev.o fftjw.o qr_class.o levmar2.o mp_mat_double.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_setd_dd : test_setd.ddo poly_grid.ddo gammaETD.ddo mp_mat.ddo str_dd.ddo \
	chebyshev.ddo fftjw.ddo qr_class.ddo levmar2.ddo mp_mat_dd.ddo \
	fft_template.ddo lapack_dd.ddo dd_utils.ddo str_double.ddo dstedc_dd.ddo \
	dgesvd_dd.ddo mp_mat_double.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_setd2 : test_setd2.o poly_grid.o gammaETD.o mp_mat.o str_double.o \
	chebyshev.o fftjw.o qr_class.o levmar2.o mp_mat_double.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_setd2_dd : test_setd2.ddo poly_grid.ddo gammaETD.ddo mp_mat.ddo str_dd.ddo \
	chebyshev.ddo fftjw.ddo qr_class.ddo levmar2.ddo mp_mat_dd.ddo \
	fft_template.ddo lapack_dd.ddo dd_utils.ddo str_double.ddo dstedc_dd.ddo \
	dgesvd_dd.ddo mp_mat_double.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

hmodel : hmodel.o poly_grid.o gammaETD2.o mp_mat.o str_double.o \
	chebyshev2.o fftjw2.o qr_class.o fft_template.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

hmodel_dd : hmodel.ddo poly_grid.ddo gammaETD2.ddo mp_mat.ddo str_dd.ddo \
	chebyshev2.ddo fftjw2.ddo qr_class.ddo fft_template.ddo lapack_dd.ddo \
	dd_utils.ddo dstedc_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_setd3 : test_setd3.o poly_grid.o gammaETD2.o mp_mat.o str_double.o \
	chebyshev2.o fftjw2.o qr_class.o levmar2.o mp_mat_double.o fft_template.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_setd3_burger : test_setd3_burger.o poly_grid.o gammaETD2.o mp_mat.o str_double.o \
	chebyshev2.o fftjw2.o qr_class.o levmar2.o mp_mat_double.o fft_template.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_setd3_dd : test_setd3.ddo poly_grid.ddo gammaETD2.ddo mp_mat.ddo str_dd.ddo \
	chebyshev2.ddo fftjw2.ddo qr_class.ddo levmar2.ddo mp_mat_dd.ddo \
	fft_template.ddo lapack_dd.ddo dd_utils.ddo str_double.ddo dstedc_dd.ddo \
	dgesvd_dd.ddo mp_mat_double.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_setd5 : test_setd5.o poly_grid.o gammaETD2.o mp_mat.o str_double.o \
	chebyshev2.o fftjw2.o mp_mat_double.o fft_template.o qr_class.o 
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_setd5_burger : test_setd5_burger.o poly_grid.o gammaETD2.o mp_mat.o str_double.o \
	chebyshev2.o fftjw2.o qr_class.o levmar2.o mp_mat_double.o fft_template.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_etd : test_etd.o gammaETD2.o mp_mat.o str_double.o chebyshev2.o fftjw2.o \
	fft_template.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_etd_dd : test_etd.ddo gammaETD.ddo mp_mat.ddo str_dd.ddo chebyshev.ddo \
	fftjw.ddo fft_template.ddo lapack_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_etd2 : test_etd2.o gammaETD2.o mp_mat.o str_double.o chebyshev2.o fftjw2.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_etd2_dd : test_etd2.ddo gammaETD.ddo mp_mat.ddo str_dd.ddo chebyshev.ddo \
	fftjw.ddo fft_template.ddo lapack_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_gammaETD : test_gammaETD.o gammaETD.o mp_mat.o str_double.o chebyshev.o \
	fftjw.o gammaETD.ddo str_dd.ddo lapack_dd.ddo mp_mat.ddo \
	dd_utils.ddo mp_mat_dd.ddo chebyshev.ddo fft_template.ddo fftjw.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

chebETD : chebETD.o mp_mat.o str_double.o chebyshev.o fftjw.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

chebETD_dd : chebETD.ddo mp_mat.ddo str_dd.ddo chebyshev.ddo fftjw.ddo \
	fft_template.ddo lapack_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_banded : test_banded.o banded.o str_double.o mp_mat.o chebyshev.o fftjw.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_nugrid : test_nugrid.o nugrid.o poly_grid.o str_double.o mp_mat.o \
	chebyshev.o fftjw.o qr_class.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_nugrid_dd : test_nugrid.ddo nugrid.ddo poly_grid.ddo str_dd.ddo mp_mat.ddo \
	chebyshev.ddo fftjw.ddo qr_class.ddo fft_template.ddo dd_utils.ddo \
	lapack_dd.ddo dstedc_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_nugrid2 : test_nugrid2.o nugrid2.o poly_grid.o str_double.o mp_mat.o \
	chebyshev2.o fftjw2.o qr_class.o fft_template.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_nugrid2_dd : test_nugrid2.ddo nugrid2.ddo poly_grid.ddo str_dd.ddo mp_mat.ddo \
	chebyshev2.ddo fftjw2.ddo qr_class.ddo fft_template.ddo dd_utils.ddo \
	lapack_dd.ddo dstedc_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_nugrid3 : test_nugrid3.o nugrid3.o poly_grid.o str_double.o mp_mat.o \
	chebyshev2.o fftjw2.o fft_template.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

ricci4 : ricci4.o poly_grid.o str_double.o graded_mesh.o mp_mat.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

ricci3 : ricci3.o poly_grid.o str_double.o graded_mesh.o mp_mat.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

ricci2 : ricci2.o poly_grid.o str_double.o graded_mesh.o mp_mat.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

ricci2_dd : ricci2.ddo poly_grid.ddo str_dd.ddo graded_mesh.ddo mp_mat.ddo \
	lapack_dd.ddo dstedc_dd.ddo dd_utils.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

ricci1 : ricci1.o chebyshev.o str_double.o graded_mesh.o mp_mat.o fftjw.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

ricci1_dd : ricci1.ddo chebyshev.ddo str_dd.ddo graded_mesh.ddo mp_mat.ddo fftjw.ddo \
	fft_template.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

geoP : geoP.o poly_grid.o str_double.o mp_mat.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

annulusP : annulusP.o poly_grid.o str_double.o mp_mat.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

annulusP_dd : annulusP.ddo poly_grid.ddo str_dd.ddo mp_mat.ddo lapack_dd.ddo \
	dstedc_dd.ddo dd_utils.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

diskP : diskP.o poly_grid.o str_double.o mp_mat.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

epstein1 : epstein1.o chebyshev.o fftjw.o str_double.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

epstein2 : epstein2.o chebyshev.o fftjw.o str_double.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

mchrist2 : mchrist2.o lerch.o poly_grid.o mp_mat.o str_double.o chebyshev.o fftjw.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

mchrist2_dd : mchrist2.ddo lerch.ddo poly_grid.ddo mp_mat.ddo str_dd.ddo chebyshev.ddo \
	fftjw.ddo lapack_dd.ddo fft_template.ddo dd_utils.ddo dstedc_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

mchrist : mchrist.o lerch.o poly_grid.o mp_mat.o str_double.o chebyshev.o fftjw.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

mchrist_dd : mchrist.ddo lerch.ddo poly_grid.ddo mp_mat.ddo str_dd.ddo chebyshev.ddo \
	fftjw.ddo lapack_dd.ddo fft_template.ddo dd_utils.ddo dstedc_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

mchrist1 : mchrist1.o poly_grid.o mp_mat.o str_double.o chebyshev.o fftjw.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

mchrist1_dd : mchrist.ddo poly_grid.ddo mp_mat.ddo str_dd.ddo chebyshev.ddo \
	fftjw.ddo lapack_dd.ddo fft_template.ddo dd_utils.ddo dstedc_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_lerch : test_lerch.o lerch.o poly_grid.o mp_mat.o str_double.o chebyshev.o fftjw.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_lerch_dd : test_lerch.ddo lerch.ddo poly_grid.ddo mp_mat.ddo str_dd.ddo chebyshev.ddo \
	fftjw.ddo lapack_dd.ddo fft_template.ddo dd_utils.ddo dstedc_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_gamma : test_gamma.o gammanew_eval.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_irk_jw_dd : test_irk_jw.ddo irk_jw.ddo mp_mat.ddo str_dd.ddo poly_grid.ddo \
	lapack_dd.ddo dstedc_dd.ddo dd_utils.ddo qr_class.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_irk_jw : test_irk_jw.o irk_jw.o mp_mat.o str_double.o poly_grid.o qr_class.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_irk_jw2_dd : test_irk_jw2.ddo irk_jw2.ddo mp_mat.ddo str_dd.ddo poly_grid.ddo \
	lapack_dd.ddo dstedc_dd.ddo dd_utils.ddo qr_class.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_irk_jw2 : test_irk_jw2.o irk_jw2.o mp_mat.o str_double.o poly_grid.o qr_class.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_irk_fm_dd : test_irk_fm.ddo irk_fm6.ddo mp_mat.ddo str_dd.ddo poly_grid.ddo \
	lapack_dd.ddo dstedc_dd.ddo dd_utils.ddo qr_class.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_irk_fm : test_irk_fm.o irk_fm6.o mp_mat.o str_double.o poly_grid.o qr_class.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

sturm3a : sturm3a.o irk_jw2.o str_double.o mp_mat.o chebyshev.o fftjw.o \
	qr_class.o mp_mat_double.o psi.o poly_grid.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

sturm3b : sturm3b.o str_double.o mp_mat.o chebyshev.o fftjw.o \
	qr_class.o mp_mat_double.o psi.o poly_grid.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

sturm3b_dd : sturm3b.ddo str_dd.ddo mp_mat.ddo chebyshev.ddo fftjw.ddo \
	qr_class.ddo mp_mat_dd.ddo psi_dd.ddo poly_grid.ddo fft_template.ddo \
	dd_utils.ddo lapack_dd.ddo dstedc_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

sturm3 : sturm3.o irk_jw2.o str_double.o mp_mat.o chebyshev.o fftjw.o \
	qr_class.o mp_mat_double.o psi.o poly_grid.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

sturm3_dd : sturm3.ddo irk_jw2.ddo str_dd.ddo mp_mat.ddo chebyshev.ddo fftjw.ddo \
	qr_class.ddo mp_mat_dd.ddo psi_dd.ddo poly_grid.ddo fft_template.ddo \
	dd_utils.ddo lapack_dd.ddo dstedc_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)


sturm2a : sturm2a.o irk_jw2.o chebyshev.o fftjw.o str_double.o mp_mat.o \
	qr_class.o mp_mat_double.o psi.o poly_grid.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

sturm2a_dd : sturm2a.ddo irk_jw2.ddo chebyshev.ddo fftjw.ddo str_dd.ddo mp_mat.ddo \
	qr_class.ddo mp_mat_dd.ddo psi_dd.ddo poly_grid.ddo fft_template.ddo \
	dd_utils.ddo lapack_dd.ddo dstedc_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

sturm2 : sturm2.o irk_jw2.o chebyshev.o fftjw.o str_double.o mp_mat.o \
	qr_class.o mp_mat_double.o psi.o poly_grid.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

sturm2_dd : sturm2.ddo irk_jw2.ddo chebyshev.ddo fftjw.ddo str_dd.ddo mp_mat.ddo \
	qr_class.ddo mp_mat_dd.ddo psi_dd.ddo poly_grid.ddo fft_template.ddo \
	dd_utils.ddo lapack_dd.ddo dstedc_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

sturm : sturm.o chebyshev.o fftjw.o str_double.o mp_mat.o mp_mat_double.o psi.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

sturm_dd : sturm.ddo chebyshev.ddo fftjw.ddo str_dd.ddo mp_mat.ddo \
	mp_mat_dd.ddo psi_dd.ddo fft_template.ddo dd_utils.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

sturm1 : sturm1.o chebyshev.o fftjw.o str_double.o mp_mat.o mp_mat_double.o psi.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

sturm1_dd : sturm1.ddo chebyshev.ddo fftjw.ddo str_dd.ddo mp_mat.ddo \
	mp_mat_dd.ddo psi_dd.ddo fft_template.ddo dd_utils.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_cheb_dat : test_cheb_dat.o chebyshev.o fftjw.o str_double.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_cheb_dat_dd : test_cheb_dat.ddo chebyshev.ddo fftjw.ddo str_dd.ddo \
	fft_template.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_chebyshev : test_chebyshev.o chebyshev.o fftjw.o str_double.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_chebyshev_dd : test_chebyshev.ddo chebyshev.ddo fftjw.ddo str_dd.ddo \
	fft_template.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_chebyshev2 : test_chebyshev2.o chebyshev2.o fftjw2.o str_double.o fft_template.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_chebyshev2_dd : test_chebyshev2.ddo chebyshev2.ddo fftjw2.ddo str_dd.ddo \
	fft_template.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_chebyshev2_mpfr : test_chebyshev2.mpo chebyshev2.mpo fftjw2.mpo str_mpfr.mpo \
	fft_template.mpo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

composite3 : composite_grid3.o poly_grid.o mp_mat.o str_double.o psi.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

composite3_dd : composite_grid3.ddo poly_grid.ddo mp_mat.ddo str_dd.ddo psi_dd.ddo \
	dd_utils.ddo dgesvd_dd.ddo lapack_dd.ddo dstedc_dd.ddo dgeqrf_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

composite2 : composite_grid2.o poly_grid.o mp_mat.o str_double.o psi.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

composite2_dd : composite_grid2.ddo poly_grid.ddo mp_mat.ddo str_dd.ddo psi_dd.ddo \
	dd_utils.ddo dgesvd_dd.ddo lapack_dd.ddo dstedc_dd.ddo dgeqrf_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

composite : composite_grid.o poly_grid.o mp_mat.o str_double.o psi.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

composite_dd : composite_grid.ddo poly_grid.ddo mp_mat.ddo str_dd.ddo psi_dd.ddo \
	dd_utils.ddo dgesvd_dd.ddo lapack_dd.ddo dstedc_dd.ddo dgeqrf_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

plasma_grid_segment_dd : plasma_grid_segment.ddo poly_grid.ddo mp_mat.ddo \
	str_dd.ddo dstedc_dd.ddo lapack_dd.ddo dd_utils.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

plasma_grid_segment_mpfr : plasma_grid_segment.mpo poly_grid.mpo mp_mat.mpo \
	mp_mat_mpfr.mpo str_mpfr.mpo dstedc_mpfr.mpo lapack_mpfr.mpo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

rafelski : rafelski.o poly_grid.o mp_mat.o str_double.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

rafelski_dd : rafelski.ddo poly_grid.ddo mp_mat.ddo str_dd.ddo dd_utils.ddo \
	lapack_dd.ddo dgesvd_dd.ddo dstedc_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

kimura : kimura.o poly_grid.o mp_mat.o str_double.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

kimura_dd : kimura.ddo poly_grid.ddo mp_mat.ddo str_dd.ddo dd_utils.ddo lapack_dd.ddo \
	dgesvd_dd.ddo dstedc_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

cheb_to_leg : cheb_to_leg.o poly_grid.o mp_mat.o str_double.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

cheb_to_leg_dd : cheb_to_leg.ddo poly_grid.ddo mp_mat.ddo str_dd.ddo \
	lapack_dd.ddo dgesvd_dd.ddo dstedc_dd.ddo dd_utils.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

cheb_to_leg_mpfr : cheb_to_leg.mpo poly_grid.mpo mp_mat.mpo str_mpfr.mpo \
	lapack_mpfr.mpo dgesvd_dd.mpo dstedc_mpfr.mpo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_poly_nugrid : test_poly_nugrid.o poly_grid.o mp_mat.o str_double.o nugrid2.o \
	chebyshev2.o fft_template.o fftjw2.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_poly_grid : test_poly_grid.o poly_grid.o mp_mat.o str_double.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_jacobi : test_jacobi.o poly_grid.o mp_mat.o str_double.o \
	chebyshev2.o fft_template.o fftjw2.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_jacobi_dd: test_jacobi.ddo poly_grid.ddo mp_mat.ddo str_dd.ddo \
	dstedc_dd.ddo dd_utils.ddo lapack_dd.ddo \
	chebyshev2.ddo fft_template.ddo fftjw2.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_general_laguerre : test_general_laguerre.o poly_grid.o mp_mat.o str_double.o \
	chebyshev2.o fft_template.o fftjw2.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_general_laguerre_dd: test_general_laguerre.ddo poly_grid.ddo mp_mat.ddo str_dd.ddo \
	dstedc_dd.ddo dd_utils.ddo lapack_dd.ddo \
	chebyshev2.ddo fft_template.ddo fftjw2.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_poly_grid_dd : test_poly_grid.ddo poly_grid.ddo mp_mat.ddo str_dd.ddo dstedc_dd.ddo \
	dd_utils.ddo lapack_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_poly_grid_mpfr : test_poly_grid.mpo poly_grid.mpo mp_mat.mpo str_mpfr.mpo dstedc_mpfr.mpo \
	lapack_mpfr.mpo mp_mat_mpfr.mpo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

sdc_computeS : sdc_computeS.mpo poly_grid.mpo mp_mat.mpo str_mpfr.mpo dstedc_mpfr.mpo \
	lapack_mpfr.mpo mp_mat_mpfr.mpo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

fdiff : fdiff.o expm.o poly_grid.o mp_mat.o str_double.o psi.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

fdiff_time : fdiff_time.o expm.o poly_grid.o mp_mat.o str_double.o psi.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

fdiff_dd : fdiff.ddo expm.ddo poly_grid.ddo mp_mat.ddo str_dd.ddo psi_dd.ddo dd_utils.ddo \
	lapack_dd.ddo dstedc_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

fdiff2 : fdiff2.o mp_mat.o str_double.o psi.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

fdiff2_dd : fdiff2.ddo mp_mat.ddo str_dd.ddo psi_dd.ddo dd_utils.ddo \
	lapack_dd.ddo dstedc_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_dd_read : test_dd_read.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

dd_diff_mat2 : dd_diff_mat2.ddo str_dd.ddo mp_mat.ddo mp_mat_dd.ddo dd_utils.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

dd_diff_mat : dd_diff_mat.ddo str_dd.ddo mp_mat.ddo mp_mat_dd.ddo dd_utils.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

repar : repar.ddo str_dd.ddo fft_template.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

dd_diff : dd_diff.ddo str_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

dd_diff_rel : dd_diff_rel.ddo str_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

dd_prod : dd_prod.ddo str_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

exact1 : exact1.ddo str_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

exact2 : exact2.ddo str_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

fit2_dd : fit2.ddo str_dd.ddo dd_utils.ddo mp_mat.ddo mp_mat_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

fit2 : fit2.o str_double.o mp_mat.o mp_mat_double.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

fit_dd : fit.ddo str_dd.ddo dd_utils.ddo mp_mat.ddo \
	         mp_mat_dd.ddo qr_class.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

fit : fit.o str_double.o mp_mat.o mp_mat_double.o qr_class.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho9_mpi : compute_rho9_mpi.ddo str_dd.ddo dd_utils.ddo mp_mat.ddo \
	         mp_mat_dd.ddo psi_dd.ddo fft_template.ddo integrate.ddo \
		 qr_class.ddo
	$(MPICC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rhoA_dd : compute_rhoA.ddo str_dd.ddo dd_utils.ddo mp_mat.ddo \
	         mp_mat_dd.ddo psi_dd.ddo fft_template.ddo integrate.ddo \
		 qr_class.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho9_dd : compute_rho9.ddo str_dd.ddo dd_utils.ddo mp_mat.ddo \
	         mp_mat_dd.ddo psi_dd.ddo fft_template.ddo integrate.ddo \
		 qr_class.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho9 : compute_rho9.o str_double.o mp_mat.o mp_mat_double.o \
	psi.o fft_template.o integrate.o qr_class.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho9_debug : compute_rho9.do str_double.do mp_mat.do mp_mat_double.do \
	psi.do fft_template.do integrate.do qr_class.do
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho8_dd : compute_rho8.ddo str_dd.ddo dd_utils.ddo mp_mat.ddo \
	         mp_mat_dd.ddo psi_dd.ddo extrap_dd.ddo qr_class.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho8 : compute_rho8.o str_double.o mp_mat.o mp_mat_double.o \
	psi.o extrap.o qr_class.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho8_debug : compute_rho8.do str_double.do mp_mat.do mp_mat_double.do \
	psi.do extrap.do qr_class.do
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho7_dd : compute_rho7.ddo str_dd.ddo dd_utils.ddo mp_mat.ddo \
	         mp_mat_dd.ddo psi_dd.ddo fft_template.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho7 : compute_rho7.o str_double.o mp_mat.o mp_mat_double.o \
	psi.o fft_template.o integrate.o qr_class.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho7_debug : compute_rho7.do str_double.do mp_mat.do mp_mat_double.do \
	psi.do fft_template.do integrate.do qr_class.do
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho6_dd : compute_rho6.ddo str_dd.ddo dd_utils.ddo mp_mat.ddo \
	         mp_mat_dd.ddo psi_dd.ddo extrap_dd.ddo qr_class.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho6 : compute_rho6.o str_double.o mp_mat.o mp_mat_double.o \
	psi.o extrap.o qr_class.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho6_debug : compute_rho6.do str_double.do mp_mat.do mp_mat_double.do \
	psi.do extrap.do qr_class.do
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho5_mpi : compute_rho5_mpi.ddo str_dd.ddo dd_utils.ddo mp_mat.ddo \
	         mp_mat_dd.ddo psi_dd.ddo fft_template.ddo integrate.ddo \
		 qr_class.ddo
	$(MPICC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho5_dd : compute_rho5.ddo str_dd.ddo dd_utils.ddo mp_mat.ddo \
	         mp_mat_dd.ddo psi_dd.ddo fft_template.ddo integrate.ddo \
		 qr_class.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho5 : compute_rho5.o str_double.o mp_mat.o mp_mat_double.o \
	psi.o fft_template.o integrate.o qr_class.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho5_debug : compute_rho5.do str_double.do mp_mat.do mp_mat_double.do \
	psi.do fft_template.do integrate.do qr_class.do
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho4_dd : compute_rho4.ddo str_dd.ddo dd_utils.ddo mp_mat.ddo \
	         mp_mat_dd.ddo psi_dd.ddo extrap_dd.ddo qr_class.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho4 : compute_rho4.o str_double.o mp_mat.o mp_mat_double.o \
	psi.o extrap.o qr_class.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho4_debug : compute_rho4.do str_double.do mp_mat.do mp_mat_double.do \
	psi.do extrap.do qr_class.do
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho3_dd : compute_rho3.ddo str_dd.ddo dd_utils.ddo mp_mat.ddo \
	         mp_mat_dd.ddo psi_dd.ddo fft_template.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho3 : compute_rho3.o str_double.o mp_mat.o mp_mat_double.o \
	psi.o fft_template.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho3_debug : compute_rho3.do str_double.do mp_mat.do mp_mat_double.do \
	psi.do fft_template.do
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho2_dd : compute_rho2.ddo str_dd.ddo dd_utils.ddo mp_mat.ddo \
	         mp_mat_dd.ddo psi_dd.ddo fft_template.ddo integrate.ddo \
		 qr_class.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho2 : compute_rho2.o str_double.o mp_mat.o mp_mat_double.o \
	psi.o fft_template.o integrate.o qr_class.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho2_debug : compute_rho2.do str_double.do mp_mat.do mp_mat_double.do \
	psi.do fft_template.do integrate.do qr_class.do
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho1_dd : compute_rho1.ddo str_dd.ddo dd_utils.ddo mp_mat.ddo \
	         mp_mat_dd.ddo psi_dd.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho1 : compute_rho1.o str_double.o mp_mat.o mp_mat_double.o \
	psi.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho1_debug : compute_rho1.do str_double.do mp_mat.do mp_mat_double.do \
	psi.do
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho_dd : compute_rho.ddo str_dd.ddo dd_utils.ddo mp_mat.ddo \
	         mp_mat_dd.ddo psi_dd.ddo extrap_dd.ddo qr_class.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho : compute_rho.o str_double.o mp_mat.o mp_mat_double.o \
	psi.o extrap.o qr_class.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compute_rho_debug : compute_rho.do str_double.do mp_mat.do mp_mat_double.do \
	psi.do extrap.do qr_class.do
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_integrate : test_integrate.o integrate.o mp_mat.o mp_mat_double.o \
	str_double.o qr_class.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_integrate_dd : test_integrate.ddo integrate.ddo mp_mat.ddo \
	mp_mat_dd.ddo str_dd.ddo dd_utils.ddo qr_class.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

# obsolete
pseudo : driverP.o pseudo.o mp_mat.o mp_mat_double.o str_double.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

driver2 : driver2.o orthog.o orthog_extra.o mp_mat.o mp_mat_double.o str_double.o sp_deriv.o \
	 psi.o qr_class.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

driver : driver.o orthog.o orthog_extra.o mp_mat.o mp_mat_double.o str_double.o sp_deriv.o \
	 psi.o qr_class.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

driver_time : driver.o orthog_time.o orthog_extra.o mp_mat.o mp_mat_double.o str_double.o sp_deriv.o \
	 psi.o qr_class.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

driver2_dd : driver2.ddo orthog.ddo orthog_extra.ddo mp_mat.ddo \
	mp_mat_dd.ddo str_dd.ddo \
	sp_deriv_dd.ddo dd_utils.ddo lapack_dd.ddo psi_dd.ddo \
	dgesvd_dd.ddo dgeqrf_dd.ddo dstedc_dd.ddo qr_class.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

driver_dd : driver.ddo orthog.ddo orthog_extra.ddo mp_mat.ddo mp_mat_dd.ddo \
	str_dd.ddo sp_deriv_dd.ddo dd_utils.ddo lapack_dd.ddo psi_dd.ddo \
	dgesvd_dd.ddo dgeqrf_dd.ddo dstedc_dd.ddo qr_class.ddo
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

driver_debug : driver.do orthog.do mp_mat.do mp_mat_double.do \
	       str_double.do sp_deriv.do
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_sp_deriv : mp_mat.o mp_mat_double.o sp_deriv.o test_sp_deriv.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_psi_dd : test_psi_dd.o psi_dd.o str_dd.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

compare_dd : cmp_dd_double.o str_dd.ddo 
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

%: %.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

%.o : ../%.f
	$(f77) -c $(CFAST) $<

%.o : ../%.c $(headers)
	$(cc) -c $(INC) $(cFLAGS) $<

%.o: ../%.cpp $(headers)
	$(CC) -c $(INC) $(CFLAGS) $<

%.o: ../%.cxx $(headers)
	$(CC) -c $(INC) $(CFLAGS) $<

%.o : ../%.cc $(headers)
	$(MPICC) -c $(INC) $(CFLAGS) $<

%.o : ../%.cu $(headers)
	$(NVCC) -c $(CUINC) $(CUFLAGS) $<

%.do: ../%.cpp $(headers)
	$(CC) -c -o $@ $(INC) $(CFLAGSD) $<

%.ddo: ../%.cpp $(headers)
	$(CC) -c -o $@ $(INC) $(CFLAGS) -DDDDD=1 $<

%.ddo: ../%.cc $(headers)
	$(MPICC) -c -o $@ $(INC) $(CFLAGS) -DDDDD=1 $<

%.ddo: ../%.cu $(headers)
	$(NVCC) -c -o $@ $(CUINC) $(CUFLAGS) -DDDDD=1 $<

%.mpo: ../%.cpp $(headers)
	$(CC) -c -o $@ $(INC) $(CFLAGS) -DMPFR=1 $<

irk_fm.o: ../irk_fm.cpp ../irk_fm_aux.cpp $(headers)
	$(CC) -c -o $@ $(INC) $(CFLAGS) $<

irk_fm.ddo: ../irk_fm.cpp ../irk_fm_aux.cpp $(headers)
	$(CC) -c -o $@ $(INC) $(CFLAGS) -DDDDD=1 $<

.PHONY : always
%.ps : always
	latex ../template; dvips -o $*.ps template

.PHONY : clean
clean:
	rm -f *.o *.ddo *.mpo *.do a.out
